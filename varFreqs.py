#!/home/rmccloskey/bin/python

import pysam
import sys
import csv
import imp
import os
import re

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
                    base_qual = ord(read.alignment.qual[read.qpos])-33
                    if base_qual >= minimum_base_qual:
                        base = read.alignment.seq[read.qpos]
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

def read_maf_file(path):
    """Read lines from MAF file, index by position and allele."""
    if not os.path.exists(path):
        return {}
    rows = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = (row["Chromosome"], 
                   row["Start_Position"], 
                   row["Reference_Allele"], 
                   row["Tumor_Seq_Allele1"])
            rows[key] = row
    return rows

def get_counts(samfile, reffile, chrom, pos, ref_base, alt_base):
    """Get the variant allele fraction at a particular position."""
    start = max(pos-200, 0)
    end = pos+200
    reffile.fetch(reference=chrom, start=pos-1, end=pos)
    ref_base = reffile.fetch(reference=chrom, start=pos-1, end=pos).decode("utf-8")
    pileup = samfile.pileup(chrom, start, end, fastafile=reffile)
    return countBasesAtPileupPosition(pileup, pos, ref_base, alt_base)

def try_append(d, key, value):
    try:
        d[key].append(value)
    except KeyError:
        d[key] = [value]
    return d

if __name__ == "__main__":

    reader = csv.DictReader(open(settings.METADATA), delimiter="\t")

    fields = ["patient", "sample", "chrom", "pos", "ref", "alt", "ref.count", 
              "alt.count", "depth", "gene", "class", "rs", "cosmic", "esp",
              "prot.change"]
    writer = csv.DictWriter(sys.stdout, delimiter="\t", fieldnames=fields)
    writer.writeheader()

    samples = {}
    for row in reader:
        samples = try_append(samples, row["patient"], row["tumor.sample"])

    maf_dir = os.path.join(settings.WORK_DIR, "11_snv")
    bam_dir = os.path.join(settings.WORK_DIR, "01_fixbams")
    reffile = pysam.Fastafile(settings.HUMAN_REF)

    for patient in samples:
        all_rows = {}
        for sample in samples[patient]:
            maf_file = os.path.join(maf_dir, "{}.maf".format(sample))
            all_rows.update(read_maf_file(maf_file))
        
        for sample in samples[patient]:
            bam_file = os.path.join(bam_dir, "{}.bam".format(sample))
            if not os.path.exists(bam_file):
                continue
            samfile = pysam.Samfile(bam_file)
            for key, in_row in all_rows.items():
                row = dict.fromkeys(fields)
                row["chrom"], row["pos"], row["ref"], row["alt"] = key
                row["pos"] = int(row["pos"])
                res = get_counts(samfile, reffile, row["chrom"], row["pos"], 
                                 in_row["Reference_Allele"], 
                                 in_row["Tumor_Seq_Allele1"])
                row["depth"], row["ref.count"], row["alt.count"] = res
                row["patient"] = patient
                row["sample"] = sample
                row["gene"] = in_row["Hugo_Symbol"]
                row["class"] = in_row["Variant_Classification"]
                
                row["cosmic"] = ",".join(re.findall("COSM\d+", in_row["dbSNP_RS"]))
                row["rs"] = ",".join(re.findall("rs\d+", in_row["dbSNP_RS"]))
                row["esp"] = ",".join(re.findall("TMP_ESP_[\dXY]{1,2}_\d+", in_row["dbSNP_RS"]))
                row["prot.change"] = in_row["Protein_Change"].replace("p.", "")
                writer.writerow(row)
