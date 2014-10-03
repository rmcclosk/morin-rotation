#!/home/rmccloskey/bin/python3

import pysam
import sys
import csv
import imp
import os

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
                    base = chr(read.alignment.seq[read.qpos])
                    base_qual = read.alignment.qual[read.qpos]-33
                    if base_qual >= minimum_base_qual:
                        counts[base] += 1

    ref_count = sum(counts[a] for a in ref_allele)
    if alt_allele is None:
        best_base, best_n = max(counts.items(), key=lambda x: (x[0] != ref_allele, x[1]))
        return ref_count, best_n, best_base, sum(counts.values()) - ref_count - best_n
    else:
        alt_count = sum(counts[a] for a in alt_allele)
        return ref_count, alt_count, alt_allele, sum(counts.values()) - ref_count - alt_count

def read_vaf_file(path):
    """Read variant allele fractions."""
    reader = csv.DictReader(open(path), delimiter=" ")
    positions = set([])
    for row in reader:
        positions.add((row["chr"], int(row["pos"]), row["ref"], row["alt"]))
    return positions

def get_vaf(samfile, reffile, chrom, pos, ref_base, alt_base):
    """Get the variant allele fraction at a particular position."""
    start = max(pos-200, 0)
    end = pos+200
    reffile.fetch(reference=chrom, start=pos-1, end=pos)
    ref_base = reffile.fetch(reference=chrom, start=pos-1, end=pos).decode("utf-8")
    try:
        pileup = samfile.pileup(chrom, start, end, fastafile=reffile)
    except ValueError:
        #continue
        raise
    ref, nonref, other, _ = countBasesAtPileupPosition(pileup, pos, ref_base, alt_base)
    if (ref + nonref == 0):
        sys.stderr.write("{} {} {} {} {} {}\n".format(samfile.filename, chrom, pos, ref_base, alt_base, other))
    return ref, nonref

if __name__ == "__main__":

    reader = csv.DictReader(open(settings.EXOME_METADATA))
    writer = csv.writer(sys.stdout)
    samples = {}
    for row in reader:
        filename = "{}_{}_exome.GRCh37-lite.aln.bam".format(row["patient_id"], row["gsc_exome_library_id"])
        try:
            samples[row["patient_id"]].append((row["sample_id"], filename))
        except KeyError:
            samples[row["patient_id"]] = [(row["sample_id"], filename)]

    vaf_dir = os.path.join(settings.WORK_DIR, "11_vaf")
    bam_dir = os.path.join(settings.WORK_DIR, "01_fixbams")
    reffile = pysam.Fastafile(settings.HUMAN_REF)

    row = ["patient", "sample", "chrom", "pos", "ref", "alt", "ref.count", "alt.count"]
    writer.writerow(row)
    for patient in samples:
        positions = set([])
        samples[patient] = sorted(samples[patient], key=lambda x:x[1].split("_")[1])
        for sample, _ in samples[patient]:
            vaf_file = os.path.join(vaf_dir, "{}.dat".format(sample))
            if os.path.exists(vaf_file):
                positions |= read_vaf_file(vaf_file)
        
        for sample, filename in samples[patient]:
            bam_file = os.path.join(bam_dir, filename)
            samfile = pysam.Samfile(bam_file)
            for chrom, pos, ref_base, alt_base in positions:
                ref_count, all_count = get_vaf(samfile, reffile, chrom, pos, ref_base.split(","), alt_base.split(","))
                row = [patient, sample, chrom, pos, ref_base, alt_base, ref_count, all_count]
                writer.writerow(row)
