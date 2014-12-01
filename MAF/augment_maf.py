#!/usr/bin/env python

import pysam
import csv
import argparse
import mafUtils
import bamUtils
import warnings

if __name__ == "__main__":
    desc = "Add normal and/or tumour allele counts to a MAF file"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--replace", "-r", action="store_true", default=False)
    parser.add_argument("--normal-bam", action="append", default=[])
    parser.add_argument("--tumour-bam", action="append", default=[])
    parser.add_argument("maf")
    parser.add_argument("reference")
    parser.add_argument("outfile")
    args = parser.parse_args()

    reffile = pysam.Fastafile(args.reference)
    normal_sams = [pysam.Samfile(bam) for bam in args.normal_bam]
    tumour_sams = [pysam.Samfile(bam) for bam in args.tumour_bam]
    samfiles = {"normal": normal_sams, "tumour": tumour_sams}
    lines = (line for line in open(args.maf) if not line.startswith("#"))
    reader = csv.DictReader(lines, delimiter="\t")

    count_keys = ["n_ref_count", "n_alt_count", "t_ref_count", "t_alt_count"]
    fields = reader.fieldnames
    for key in count_keys:
        if not key in fields:
            fields.append(key)

    outfile = open(args.outfile, "w")
    writer = csv.DictWriter(outfile, delimiter="\t", fieldnames=fields)
    writer.writeheader()

    for row in reader:
        chrom = row["Chromosome"]
        pos = int(row["Start_Position"])
        old_counts = mafUtils.get_allele_counts(row)
        if sum(old_counts) > 0:
            if args.replace:
                warnings.warn("Replacing existing allele counts")
                old_counts = (0, 0, 0, 0)
            else:
                warnings.warn("Adding to existing allele counts")
        
        row.update(dict(zip(count_keys, old_counts)))
        ref = row["Reference_Allele"]
        alt = mafUtils.get_nref_allele(row)

        for sample, sams in samfiles.items():
            ref_key = "{}_ref_count".format(sample[0])
            alt_key = "{}_alt_count".format(sample[0])
            for samfile in sams:
                if mafUtils.is_snv(row):
                    counts = bamUtils.count_bases(samfile, reffile, chrom, pos)
                else:
                    counts = bamUtils.count_indels(samfile, reffile, chrom, pos, ref, alt)
                row[ref_key] += counts[ref]
                row[alt_key] += counts[alt]

        writer.writerow(row)
