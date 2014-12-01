#!/usr/bin/env python

# Import MAF files into the database. Optionally, obtain counts of reference
# and non-reference alleles from a tumor and/or normal BAM file.
#
# This script requires the pysam and BioPython libraries to be installed. 

# There is code in here to count allele frequencies for small indels, which is
# currently unused, since the indel table in the database does not have fields
# for allele counts. Please leave the code in, since it will probably be needed
# in the future.

import pysam
import argparse
import cancer_db as cancerGenome
import csv
import warnings
import subprocess
import tempfile
import sys
import re
from Bio import SeqIO

###############################################################################
# Functions for counting alleles                                              #
###############################################################################

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

################################################################################
# Functions for parsing the MAF file                                           #
################################################################################

def get_protein_change(maf_row):
    """Parse the protein change from a MAF row"""
    try:
        prot_change = maf_row["Protein_Change"]
    except KeyError:
        prot_change = maf_row["HGVSp_Short"]
    return prot_change.replace("p.", "").replace("=", "")

def get_nref_allele(maf_row):
    """Get the nonreference allele from a MAF row"""
    ref_allele = maf_row["Reference_Allele"]
    if maf_row["Tumor_Seq_Allele1"] != ref_allele:
        return maf_row["Tumor_Seq_Allele1"]
    return maf_row["Tumor_Seq_Allele2"]
    
def is_snv(maf_row):
    """True if this row describes a SNV, False for an indel."""
    ref_allele = maf_row["Reference_Allele"]
    nref_allele = get_nref_allele(maf_row)
    return (len(ref_allele) == len(nref_allele) and 
            "-" not in [ref_allele, nref_allele])

def get_ensembl_id(maf_row):
    """Get the Ensembl gene ID for a MAF row."""
    if row["Entrez_Gene_Id"].startswith("ENSG"):
        return row["Entrez_Gene_Id"]
    return row["Gene"]

def get_base_change(maf_row):
    """Get the base change from a MAF row describing a SNV"""
    ref_allele = maf_row["Reference_Allele"]
    nref_allele = get_nref_allele(maf_row)
    return "{}>{}".format(ref_allele, nref_allele)

def get_cdna_change(maf_row):
    """Get the cDNA change from a MAF row describing a SNV"""
    try:
        return row["cDNA_Change"]
    except KeyError:
        cdna_ptn = "c[.](\d+)([ATCG])>([ATCG])"
        match = re.search(cdna_ptn, maf_row["HGVSc"])
        if match:
            return "{}{}{}".format(*match.group(2,1,3))
        return ""

def get_triplet(maf_row):
    """Get the triplet context from a MAF row describing a SNV"""
    try:
        return maf_row["Codons"].split("/")[0]
    except KeyError:
        return ""
    
def get_transcript(maf_row):
    """Get the transcript from a MAF row"""
    try:
        return row["Annotation_Transcript"]
    except KeyError:
        return row["Transcript_ID"]

def get_identifiers(maf_row):
    """Get all known identifiers from a MAF row"""
    try:
        ids = row["Existing_variation"].replace(",", "&")
    except KeyError:
        ids = row["dbSNP_RS"]
    ids = ids.replace(",", " ").replace(";", " ")
    ids = ids.split()
    return "&".join(id for id in ids if id != "novel")

def get_codon_pos(maf_row):
    """Get the codon position of a MAF row representing an indel"""
    try:
        positions = re.findall("\d+", maf_row["Transcript_Position"])
        positions = [int((int(i)-1)/3)+1 for i in positions]
    except KeyError:
        positions = [int(i) for i in re.findall("\d+", maf_row["HGVSc"])]

    if len(positions) == 1:
        return positions[0]
    else:
        return "{}-{}".format(min(positions), max(positions))
        
def get_effect(maf_row):
    """Get the effect from a MAF row representing an indel"""
    try:
        return maf_row["Consequence"]
    except KeyError:
        indel_class = maf_row["Variant_Classification"]
        if indel_class.startswith("Frame_Shift"):
            return "frameshift_variant"
        elif indel_class == "In_Frame_Del":
            return "inframe_deletion"
        elif indel_class == "In_Frame_Ins":
            return "inframe_insertion"
        else:
            msg = "Inserting NULL effect for an indel of type {}"
            warnings.warn(msg.format(indel_class))
            return ""

def get_allele_counts(maf_row):
    """Get allele counts from a row of a MAF file


    Returns a 4-tuple (normal_ref, normal_nref, tumour_ref, tumour_nref)"""
    counts = []
    for key in ["n_ref_count", "n_alt_count", "t_ref_count", "t_alt_count"]:
        try:
            counts.append(int(maf_row[key]))
        except (KeyError, ValueError):
            counts.append(0)
    return counts
         

################################################################################
# Main                                                                         #
################################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="upload MAF file to database")
    parser.add_argument("maf", help="MAF file to add to database")
    parser.add_argument("library_id", type=int, help="library ID")
    parser.add_argument("db_host", help="database host")
    parser.add_argument("db_name", help="database name")
    parser.add_argument("db_user", help="database user")
    parser.add_argument("db_pass", help="database password")
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
    cursor = db.db.cursor()

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

    reader = csv.DictReader((line for line in open(args.maf) if not line.startswith("#")), delimiter="\t")
    num_ptn = re.compile("\d+")
    cdna_ptn = re.compile("c[.](?P<pos>\d+)(?P<ref>[ATCG])>(?P<alt>[ATCG])")
    
    allowed_classes = ["Frame_Shift_Del",
                       "Frame_Shift_Ins",
                       "In_Frame_Del", 
                       "In_Frame_Ins", 
                       "Missense_Mutation", 
                       "Nonsense_Mutation", 
                       "Silent", 
                       "Splice_Site", 
                       "Translation_Start_Site", 
                       "Nonstop_Mutation", 
                       "RNA", 
                       "Targeted_Region", 
                       "De_novo_Start_InFrame", 
                       "De_novo_Start_OutOfFrame"]
    for row in reader:
        variant_class = row["Variant_Classification"]
        if not variant_class in allowed_classes:
            warnings.warn("Skipping mutation of type {}".format(variant_class))
            continue

        params = {"library_id": args.library_id, 
                  "chromosome": row["Chromosome"]}
    
        # SNV
        if is_snv(row):
            ref_base = row["Reference_Allele"]
            nref_base = get_nref_allele(row)
            protein_change = get_protein_change(row)
            params.update({"position": row["Start_Position"],
                           "ensembl_gene_id": get_ensembl_id(row),
                           "base_change": get_base_change(row),
                           "protein_altering": protein_change != "",
                           "annotation": protein_change,
                           "cdna_change": get_cdna_change(row),
                           "identifiers": get_identifiers(row),
                           "splice_site": variant_class == "Splice_Site",
                           "triplet": get_triplet(row),
                           "sift_score": 0,
                           "polyphen_score": 0,
                           "mutation_ass_score": 0,
                           "transcript": get_transcript(row)})
        
            keys = ["normal_ref", "normal_nref", "tumour_ref", "tumour_nref"]
            params.update(dict(zip(keys, get_allele_counts(row))))

            for sample, samfile in samfiles.items():
                ref_key = "{}_ref".format(sample)
                alt_key = "{}_nref".format(sample)
                if samfile:
                    chrom = row["Chromosome"]
                    pos = row["Start_Position"]
                    counts = count_bases(samfile, reffile, chrom, int(pos))
                    params[ref_key] = counts[ref_base]
                    params[alt_key] = counts[nref_base]
            db.addMutation(**params)
    
        # indel
        else:
            params.update({"annotation": get_codon_pos(row),
                           "ensembl_id": get_ensembl_id(row),
                           "start": row["Start_Position"],
                           "end": row["End_Position"],
                           "ref": row["Reference_Allele"],
                           "alt": get_nref_allele(row),
                           "effect": get_effect(row)})
            db.addIndel(**params)
