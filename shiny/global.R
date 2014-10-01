#!/usr/bin/env Rscript

classes <- c(rep("factor", 3), "numeric", rep("factor", 2), "numeric")
d <- read.csv("freqs.dat", header=T, colClasses=classes, as.is=T)
d <- d[d$chrom %in% c(1:22, "X", "Y"),]
d$chrom <- factor(d$chrom, levels=c(1:22, "X", "Y"))

# keep patients and samples in the same order
d$patient <- factor(d$patient, levels=unique(d$patient))
d$sample <- factor(d$sample, levels=unique(d$sample))

nunique <- function (x) length(unique(x))
upto <- function (n) 1:n
sample.counts <- aggregate(sample~patient, d, nunique)
sample.idx <- unlist(sapply(sample.counts$sample, upto))
sample.counts <- data.frame(sample=levels(d$sample), sample.num=sample.idx)

d <- merge(sample.counts, d)
