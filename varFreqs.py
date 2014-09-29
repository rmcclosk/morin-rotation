#!/home/rmccloskey/bin/python3

import pysam
import sys
import csv
import imp
import os

settings = imp.load_source("settings", "/home/rmccloskey/morin-rotation/settings.conf")

def countBasesAtPileupPosition(pileup_object, 
                               position, 
                               ref_allele, 
                               nonref_allele=None,
                               minimum_base_qual=10,
                               minimum_mapping_qual=1):

    counts = {'A':0,'T':0,'C':0,'G':0}
    for x in pileup_object:
        if position == x.pos + 1:
            for read in x.pileups:
                if read.alignment.mapq >= minimum_mapping_qual:
                    base = chr(read.alignment.seq[read.qpos])
                    #weird that this is already done for mapping quality 
                    #base_qual = ord(read.alignment.qual[read.qpos])-33 
                    base_qual = read.alignment.qual[read.qpos]-33
                    if base_qual >= minimum_base_qual:
                        counts[base] += 1

    if nonref_allele:
        total = sum(counts.values())
        other_counts = total - counts[ref_allele] - counts[nonref_allele]
        return (counts[ref_allele], counts[nonref_allele], other_counts)

    else:
        best_base, best_n = max(counts.items(), key=lambda x: (x[0] != ref_allele, x[1]))
        #return the nonref allele since it was unknown
        return (counts[ref_allele], best_n, best_base) 

def read_vaf_file(path):
    """Read variant allele fractions."""
    reader = csv.DictReader(open(path), delimiter=" ")
    positions = set([])
    for row in reader:
        positions.add((row["chr"], int(row["pos"]), float(row["vaf"])))
    return positions

if __name__ == "__main__":

    reader = csv.DictReader(open(settings.EXOME_METADATA))
    samples = {}
    for row in reader:
        try:
            samples[row["patient_id"]].append(row["sample_id"])
        except KeyError:
            samples[row["patient_id"]] = [row["sample_id"]]

    vaf_dir = os.path.join(settings.WORK_DIR, "11_vaf")
    for patient in samples:
        positions = set([])
        for sample in samples[patient]:
            vaf_file = os.path.join(vaf_dir, "{}.dat".format(sample))
            if os.path.exists(vaf_file):
                positions |= read_vaf_file(vaf_file)

        for pos in positions:

    #if len(sys.argv) != 5:
    #    print("Usage: varFreqs.py [bamfile] [reference fasta] [chrom] [positions]")
    #    sys.exit()
    #bam, genome_ref, chrom, positions = sys.argv[1:]
    #positions = [int(i) for i in positions.split(",")]
    
    dir = "/genesis/extscratch/morinlab/shared/rmccloskey"
    bam = dir + "/colorectal/01_fixbams/01-050_A45028_exome.GRCh37-lite.aln.bam"
    genome_ref = "/genesis/extscratch/morinlab/shared/rmccloskey/GRCh37-lite.fa" 
    chrom = "14" 
    positions = [90863402]

    for pos in positions:
        start = pos-200
        end = pos+200
        samfile = pysam.Samfile(bam)
        reffile = pysam.Fastafile(genome_ref)
        reffile.fetch(reference=chrom, start=pos-1, end=pos)
        ref_base = reffile.fetch(reference=chrom, start=pos-1, end=pos).decode("utf-8")
        try:
            pileup = samfile.pileup(chrom, start, end, fastafile=reffile)
        except ValueError:
            #continue
            raise
        (ref, nonref, other_base) = countBasesAtPileupPosition(pileup, 
                                                               pos,
                                                               ref_base,
                                                               minimum_base_qual=15,
                                                               minimum_mapping_qual=20)
        if nonref > 0:
            #VAF is ref / (nonref + ref)
            print("{} {},{} {}>{}".
                  format(pos, ref, nonref, ref_base, other_base))

    sys.exit()

    reader = csv.DictReader(open(settings.EXOME_METADATA))
    samples = {}
    for row in reader:
        try:
            samples[row["patient_id"]].append(row["sample_id"])
        except KeyError:
            samples[row["patient_id"]] = [row["sample_id"]]

    vaf_dir = os.path.join(settings.WORK_DIR, "11_vaf")
    for patient in samples:
        positions = []
        for sample in samples[patient]:
            vaf_file = os.path.join(vaf_dir, "{}.dat".format(sample))
            if os.path.exists(vaf_file):
                vaf_reader = csv.DictReader(open(vaf_file))
