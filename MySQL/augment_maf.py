#!/usr/bin/env python

# Update a MAF file with variant and reference allele frequencies from a BAM
# file.

import pysam
import csv
import argparse
import warnings
import subprocess
import tempfile
from Bio import SeqIO

NUCLEOTIDES = ["A", "T", "C", "G"]
NA_STRINGS = ["NA", "None", ""]

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

def count_bases_pileup(pileup, position, min_baseq=15, min_mapq=20):
    """
    Count the number of times each nucleotide occurs at a given position in a
    pileup object, subject to minimum base quality and map quality constraints.
    Return a dictionary keyed by nucleotides.
    """
    counts = dict.fromkeys(NUCLEOTIDES, 0)
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

    desc = "Update allele counts in a MAF file"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("bam_type", metavar="T", choices=["normal", "tumour"])
    parser.add_argument("maf", metavar="M")
    parser.add_argument("bam", metavar="B")
    parser.add_argument("reference", metavar="R")
    parser.add_argument("outfile", metavar="O")
    args = parser.parse_args()

    reffile = pysam.Fastafile(args.reference)
    samfile = pysam.Samfile(args.bam)
    reader = csv.DictReader(open(args.maf), delimiter="\t")
    writer = csv.DictWriter(open(args.outfile, "w"), delimiter="\t", fieldnames=reader.fieldnames)
    writer.writeheader()

    for row in reader:
        counts = count_bases(samfile, reffile, row["Chromosome"], int(row["Start_position"]))
        ref_key = "{}_ref_count".format(args.bam_type[0])
        alt_key = "{}_alt_count".format(args.bam_type[0])
        try:
            row[ref_key] = int(row[ref_key]) + counts[row["Reference_Allele"]]
        except ValueError:
            row[ref_key] = counts[row["Reference_Allele"]]

        try:
            row[alt_key] = int(row[alt_key]) + counts[row["Tumor_Seq_Allele1"]]
        except ValueError:
            row[alt_key] = counts[row["Tumor_Seq_Allele1"]]
        writer.writerow(row)
