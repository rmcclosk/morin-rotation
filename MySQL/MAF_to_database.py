#!/usr/bin/env python

# Import MAF files into the database.

import pysam
import argparse
import cancerGenome
import csv
import warnings
import subprocess
import tempfile
import sys
import re
from Bio import SeqIO

def num_mismatches(ref, query):
    """Number of matches between query and reference."""
    mismatches = 0
    for i in range(len(ref)):
        if query[i] != ref[i]:
            mismatches += 1
    return mismatches

def multi_align(ref1, ref2, queries):
    """
    Do MSA of several query sequences to two references. Return the reference
    index with the fewest mismatches for each sequence.
    """
    tfile = tempfile.NamedTemporaryFile()
    tfile.write(">ref1\n{}\n".format(ref1))
    tfile.write(">ref2\n{}\n".format(ref2))
    n_queries = 0
    for query in queries:
        tfile.write(">query{}\n{}\n".format(n_queries, query))
        n_queries += 1
    tfile.flush()
    
    outfile = tempfile.TemporaryFile()
    cmd = ["mafft", "--auto", tfile.name]
    subprocess.Popen(cmd, stdout=outfile, stderr=subprocess.PIPE).communicate()
    outfile.seek(0)
    align = SeqIO.parse(outfile, "fasta")

    ref1 = next(align)
    ref2 = next(align)
    scores = [0, 0]
    ties = 0
    for query in align:
        mismatch1 = num_mismatches(ref1, query)
        mismatch2 = num_mismatches(ref2, query)
        if mismatch1 < mismatch2:
            scores[0] += 1
        elif mismatch1 > mismatch2:
            scores[1] += 1
        else:
            ties += 1
    if ties > 0:
        msg = "Discarding {} reads aligning equally to ref and alt".format(ties)
        warnings.warn(msg)

    return scores

def count_indels(samfile, reffile, chrom, pos, ref, alt, min_mapq=20):
    """
    Count occurences of the reference and alternate indel allele at a given
    position, by alignment score.
    """
    start = max(pos-100, 0)
    end = pos+100
    ref_seq = reffile.fetch(reference=chrom, start=start, end=end).decode("utf-8")
    alt_seq = ref_seq[:pos-start-1] + alt + ref_seq[pos-start-1+len(ref):]
    
    reads = samfile.fetch(chrom, pos, pos+len(ref))
    counts = [0, 0]
    
    reads = [r.seq for r in reads if r.mapq >= min_mapq]
    scores = multi_align(ref_seq, alt_seq, reads)
    
    return scores

def count_bases_pileup(pileup, position, min_baseq=15, min_mapq=20):
    """
    Count the number of times each nucleotide occurs at a given position in a
    pileup object, subject to minimum base quality and map quality constraints.
    Return a dictionary keyed by nucleotides.
    """
    counts = dict.fromkeys(["A", "T", "C", "G"], 0)
    for x in pileup:
        if position == x.pos + 1:
            for read in x.pileups:
                dup = read.alignment.is_duplicate
                qc_fail = read.alignment.is_qcfail
                low_mapq = read.alignment.mapq < min_mapq
                if not (dup or qc_fail or low_mapq):
                    base_qual = ord(read.alignment.qual[read.qpos])-33
                    if base_qual >= min_baseq:
                        base = read.alignment.seq[read.qpos]
                        counts[base] += 1

    return counts

def count_bases(samfile, reffile, chrom, pos):
    """
    Count the number of times each nucleotide occurs at a given position in a
    bam file. Return a dictionary keyed by nucleotides.
    """
    start = max(pos-200, 0)
    end = pos+200
    reffile.fetch(reference=chrom, start=pos-1, end=pos)
    ref_base = reffile.fetch(reference=chrom, start=pos-1, end=pos).decode("utf-8")
    pileup = samfile.pileup(chrom, start, end, fastafile=reffile)
    return count_bases_pileup(pileup, pos)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="upload MAF file to database")
    parser.add_argument("maf", metavar="M", help="MAF file to add to database")
    parser.add_argument("library_id", metavar="L", type=int, help="library ID")
    parser.add_argument("db_host", metavar="H", help="database host")
    parser.add_argument("db_name", metavar="D", help="database name")
    parser.add_argument("db_user", metavar="U", help="database user")
    parser.add_argument("db_pass", metavar="P", help="database password")
    parser.add_argument("--normal-bam", dest="normal_bam", help="normal BAM file")
    parser.add_argument("--tumour-bam", dest="tumour_bam", help="tumour BAM file")
    parser.add_argument("--reference", dest="reference", help="reference FASTA file")
    args = parser.parse_args()

    db = cancerGenome.cancerGenomeDB(
            database_name=args.db_name,
            database_host=args.db_host,
            database_user=args.db_user,
            database_password=args.db_pass
    )

    if args.normal_bam or args.tumour_bam:
        if args.reference:
            reffile = pysam.Fastafile(args.reference)
        else:
            sys.exit("You must supply a reference when using BAM files")
    
    samfiles = {"tumour": None, "normal": None}
    if args.normal_bam:
        samfiles["normal"] = pysam.Samfile(args.normal_bam)
    if args.tumour_bam:
        samfiles["tumour"] = pysam.Samfile(args.tumour_bam)

    cursor = db.db.cursor()
    reader = csv.DictReader((line for line in open(args.maf) if not line.startswith("#")), delimiter="\t")
    num_ptn = re.compile("\d+")
    cdna_ptn = re.compile("c[.](?P<pos>\d+)(?P<ref>[ATCG])>(?P<alt>[ATCG])")
    
    for row in reader:
        params = {"library_id": args.library_id}
        params["chromosome"] = row["Chromosome"]
        ref_base = row["Reference_Allele"]
        nref_base = row["Tumor_Seq_Allele2"]
        protein_change = row["HGVSp_Short"].replace("p.", "").replace("=", "")
    
        # SNV
        if len(ref_base) == len(nref_base) and ref_base != "-" and nref_base != "-":
            protein_change = row["HGVSp_Short"].replace("p.", "").replace("=", "")
            params["annotation"] = protein_change
            params["position"] = row["Start_Position"]
            params["ensembl_gene_id"] = row["Entrez_Gene_Id"]
            params["base_change"] = "{}>{}".format(ref_base, nref_base)
            params["protein_altering"] = protein_change != ""
            match = cdna_ptn.search(row["HGVSc"])
            if match:
                params["cdna_change"] = "{}{}{}".format(match.group("ref"), match.group("pos"), match.group("alt"))
            else:
                params["cdna_change"] = ""
            print(row["HGVSc"], params["cdna_change"])
            params["identifiers"] = row["dbSNP_RS"].replace("&", ",")
            params["splice_site"] = row["Variant_Classification"] == "Splice_Site"
            params["mutation_seq_probability"] = 0
            params["sift_score"] = 0
            params["polyphen_score"] = 0
            params["mutation_ass_score"] = 0
    
            params["tumour_ref"] = row["t_ref_count"]
            params["tumour_nref"] = row["t_alt_count"]
            params["normal_ref"] = row["n_ref_count"]
            params["normal_nref"] = row["n_alt_count"]
            for sample, samfile in samfiles.items():
                ref_key = "{}_ref".format(sample)
                alt_key = "{}_nref".format(sample)
                if samfile:
                    chrom = row["Chromosome"]
                    pos = row["Start_Position"]
                    counts = count_bases(samfile, reffile, chrom, int(pos))
                    params[ref_key] = counts[row["Reference_Allele"]]
                    params[alt_key] = counts[row["Tumor_Seq_Allele2"]]
            db.addMutation(**params)
    
        # indel
        else:
            if protein_change == "":
                params["annotation"] = ""
            else:
                params["annotation"] = num_ptn.search(protein_change).group()

            if row["Entrez_Gene_Id"] != "0":
                query = "SELECT id FROM gene WHERE ensembl_id = '%s'"
                cursor.execute(query, row["Entrez_Gene_Id"])
                params["ensembl_id"] = cursor.fetchone()[0]
            else:
                params["ensembl_id"] = None

            params["start"] = row["Start_Position"]
            params["end"] = row["End_Position"]
            params["ref"] = ref_base
            params["alt"] = nref_base
            # TODO: update the allowed values of indel.effect in the database
            if row["Consequence"] in ('frameshift_variant','inframe_deletion','inframe_insertion','splice_donor_variant','splice_acceptor_variant','coding_sequence_variant'):
                params["effect"] = row["Consequence"]
            else:
                params["effect"] = ""
            db.addIndel(**params)
