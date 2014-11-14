#!/usr/bin/env python

import pysam
import sys
import csv
import imp
import os
import re
import argparse
import warnings
import subprocess
import tempfile
import itertools

NUCLEOTIDES = ["A", "T", "C", "G"]
NA_STRINGS = ["NA", "None", ""]

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

    desc = "Augment counts in a MAF file with another BAM file"
    parser = argparse.ArgumentParser(description=desc)
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
        row["t_ref_count"] = int(row["t_ref_count"]) + counts[row["Reference_Allele"]]
        row["t_alt_count"] = int(row["t_alt_count"]) + counts[row["Tumor_Seq_Allele1"]]
        writer.writerow(row)
