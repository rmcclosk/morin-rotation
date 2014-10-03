#!/home/rmccloskey/bin/python3

import pysam
import sys
import csv
import imp
import os
import itertools

settings = imp.load_source("settings", "/home/rmccloskey/morin-rotation/settings.conf")
NUCLEOTIDES = ["A", "T", "C", "G"]

def countBasesAtPileupPosition(pileup_object, 
                               position, 
                               ref_allele, 
                               alt_allele=None,
                               minimum_base_qual=15,
                               minimum_mapping_qual=20):

    counts = dict.fromkeys(NUCLEOTIDES, 0)
    for x in pileup_object:
        if position == x.pos + 1:
            for read in x.pileups:
                dup = read.alignment.is_duplicate
                qc_fail = read.alignment.is_qcfail
                low_mapq = read.alignment.mapq < minimum_mapping_qual
                if not (dup or qc_fail or low_mapq):
                    base_qual = read.alignment.qual[read.qpos]-33
                    if base_qual >= minimum_base_qual:
                        base = chr(read.alignment.seq[read.qpos])
                        counts[base] += 1

    ref_count = sum(counts[a] for a in ref_allele)
    alt_count = sum(counts[a] for a in alt_allele)
    return sum(counts.values()), ref_count, alt_count

def read_vcf_file(path):
    """Get positions from VCF file."""
    vcf_fields = ["chrom", "pos", "id", "ref", "alt", "qual", "filter", "info", 
                  "format", "normal", "tumor"]
    reader = csv.DictReader(filter(lambda row: row[0]!='#', open(path)), 
                            delimiter="\t", fieldnames=vcf_fields)
    positions = set([])
    for row in reader:
        positions.add((row["chrom"], int(row["pos"]), row["ref"], row["alt"]))
    return positions

def get_vcf(samfile, reffile, chrom, pos, ref_base, alt_base):
    """Get the variant allele fraction at a particular position."""
    start = max(pos-200, 0)
    end = pos+200
    reffile.fetch(reference=chrom, start=pos-1, end=pos)
    ref_base = reffile.fetch(reference=chrom, start=pos-1, end=pos).decode("utf-8")
    pileup = samfile.pileup(chrom, start, end, fastafile=reffile)
    return countBasesAtPileupPosition(pileup, pos, ref_base, alt_base)

if __name__ == "__main__":

    reader = csv.DictReader(open(settings.METADATA), delimiter="\t")
    writer = csv.writer(sys.stdout, delimiter="\t")
    samples = {}
    for row in reader:
        if row["normal.sample"]:
            try:
                samples[row["patient"]].append(row["normal.sample"])
            except KeyError:
                samples[row["patient"]] = [row["normal.sample"]]

        try:
            samples[row["patient"]].append(row["tumor.sample"])
        except KeyError:
            samples[row["patient"]] = [row["tumor.sample"]]

    vcf_dir = os.path.join(settings.WORK_DIR, "11_strelka")
    bam_dir = os.path.join(settings.WORK_DIR, "01_fixbams")
    reffile = pysam.Fastafile(settings.HUMAN_REF)

    row = ["patient", "sample", "chrom", "pos", "ref", "alt", "ref.count", "alt.count", "depth"]
    writer.writerow(row)
    for patient in samples:
        positions = set([])
        for sample in samples[patient]:
            vcf_file = os.path.join(vcf_dir, "{}.vcf".format(sample))
            if os.path.exists(vcf_file):
                positions |= read_vcf_file(vcf_file)
        
        for sample in samples[patient]:
            bam_file = os.path.join(bam_dir, "{}.bam".format(sample))
            samfile = pysam.Samfile(bam_file)
            for chrom, pos, ref_base, alt_base in positions:
                depth, ref_count, alt_count = get_vcf(samfile, reffile, chrom, pos, ref_base.split(","), alt_base.split(","))
                row = [patient, sample, chrom, pos, ref_base, alt_base, ref_count, alt_count, depth]
                writer.writerow(row)
