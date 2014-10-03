#!/usr/bin/env Rscript

classes <- c(rep("factor", 3), "numeric", rep("factor", 2), rep("numeric", 3))
d <- read.table("freqs.dat", header=T, colClasses=classes)
d <- d[d$chrom %in% c(1:22, "X", "Y"),]
d$chrom <- factor(d$chrom, levels=c(1:22, "X", "Y"))
d <- d[d$depth > 0,]
d$vaf <- d$alt.count/d$depth

# keep patients and samples in the same order
d$patient <- factor(d$patient, levels=unique(d$patient))
d$sample <- factor(d$sample, levels=unique(d$sample))

# mark sample 1, 2, ... for each patient
nunique <- function (x) length(unique(x))
upto <- function (n) 1:n
sample.counts <- aggregate(sample~patient, d, nunique)
sample.idx <- unlist(sapply(sample.counts$sample, upto))
keep.patients <- sample.counts[sample.counts$sample > 2,"patient"]
sample.counts <- data.frame(patient=rep(sample.counts$patient, sample.counts$sample), 
                            sample=levels(d$sample), 
                            sample.num=sample.idx)
sample.counts <- sample.counts[sample.counts$patient %in% keep.patients,]
d <- merge(sample.counts, d)

# take out normal samples
d <- d[d$sample.num > 1,]
d$patient <- factor(d$patient, levels=unique(d$patient))
d$sample <- factor(d$sample, levels=unique(d$sample))

# give a unique key so we can make a tooltip
d$key <- 1:nrow(d)
