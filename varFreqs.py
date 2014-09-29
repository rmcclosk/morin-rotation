#!/usr/bin/env python3

import pysam
import sys

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
                    base = read.alignment.seq[read.qpos]
                    #weird that this is already done for mapping quality 
                    base_qual = ord(read.alignment.qual[read.qpos])-33 
                    if base_qual >= minimum_base_qual:
                        counts[base] += 1

    if nonref_allele:
        total = sum(counts.values())
        other_counts = total - counts[ref_allele] - counts[nonref_allele]
        return (counts[ref_allele], counts[nonref_allele], other_counts)

    else:
        best_base, best_n = max(counts.items(), key=lambda x: x[1])
        #return the nonref allele since it was unknown
        return (counts[ref_allele], best_n, best_base) 

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: varFreqs.py [bamfile] [reference fasta] [positions]")
        sys.exit()

    bam, genome_ref, positions = sys.argv[1:]
    positions = [int(i) for i in positions.split(",")]

    for pos in positions:
        start = pos-200
        end = pos+200
        samfile = pysam.Samfile(bam, "rb")
        reffile = pysam.Fastafile(genome_ref)
        ref_base = reffile.fetch(reference="12",start=pos-1,end=pos)
        try:
            pileup = samfile.pileup(chrom, start, end, fastafile=reffile)
        except ValueError:
            continue
        (ref, nonref, other_base) = countBasesAtPileupPosition(pileup, 
                                                               pos, 
                                                               ref_base)
        if nonref > 0:
            #VAF is ref / (nonref + ref)
            print("{} {},{} {}>{}".
                  format(pos, ref, nonref, ref_base, other_base))
