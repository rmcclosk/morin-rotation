#!/usr/bin/env Rscript

d <- read.table("annotations.dat", header=T)
d$length <- d$end-d$start+1
d$exome.cnv <- d$exome.copy != 2
d$genome.cnv <- d$genome.copy != 2
d$exome.cnv.length <- d$exome.cnv*d$length
d$genome.cnv.length <- d$genome.cnv*d$length

total.length <- aggregate(length~sample+purity, d, sum)
exome.cnv.length <- aggregate(exome.cnv.length~sample+purity, d, sum)

agg.d <- merge(total.length, exome.cnv.length)
agg.d$prop.exome.cnv <- agg.d$exome.cnv.length/agg.d$length

agg.d <- agg.d[agg.d$sample == "01-2-049-9",]
