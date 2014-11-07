#!/usr/bin/env Rscript

library(Rsamtools)
library(parallel)

source(file="/home/rmccloskey/morin-rotation/settings.conf")

metadata <- read.table(METADATA, header=T, stringsAsFactors=F)
genome.dir <- file.path(sub("exome", "genome", WORK_DIR), "01_fixbams")
exome.dir <- file.path(sub("genome", "exome", WORK_DIR), "01_fixbams")
pooled.dir <- file.path(sub("(ex|gen)ome", "pooled", WORK_DIR), "01_fixbams")
dir.create(pooled.dir, showWarnings=F)

samples <- unique(c(metadata$tumor.sample, metadata$normal.sample))
samples <- samples[!is.na(samples)]
genome.bams <- file.path(genome.dir, paste0(samples, ".bam"))
exome.bams <- file.path(exome.dir, paste0(samples, ".bam"))
keeps <- file.exists(genome.bams) & file.exists(exome.bams)
genome.bams <- genome.bams[keeps]
names(genome.bams) <- samples[keeps]
exome.bams <- exome.bams[keeps]
names(exome.bams) <- samples[keeps]

intervals <- read.table(INTERVAL_LIST, stringsAsFactors=F)
split <- strsplit(intervals$V1, "[:-]")
intervals$chrom <- factor(sapply(split, "[[", 1), levels=c(1:22, "X", "Y"))
intervals$start <- as.integer(sapply(split, "[[", 2))
intervals$end <- as.integer(sapply(split, "[[", 3))

ranges <- as.list(by(intervals, intervals$chrom, function (x) {
    IRanges(start=x$start, end=x$end)
}))
ranges <- RangedData(do.call(RangesList, ranges))

samples <- samples[1]
mclapply(samples, function (sample) {
    in.bams <- c(exome.bams[sample], genome.bams[sample])
    out.bam <- file.path(pooled.dir, paste0(sample, ".bam"))
    print(out.bam)
    if (!file.exists(out.bam))
        mergeBam(in.bams, out.bam, region=ranges)
}, mc.cores=1)
