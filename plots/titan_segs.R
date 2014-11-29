#!/usr/bin/env Rscript

options(warn=1)

all.sample.groups <- function (x) {
    if (length(unique(x$time.point)) == 1)
        return(matrix(x$sample, ncol=1))
    min.time.point <- min(x$time.point)
    first <- subset(x, time.point==min.time.point)
    res <- all.sample.groups(subset(x, time.point>min.time.point))
    do.call(rbind, lapply(first$sample, function (s) {
        cbind(rep(s, nrow(res)), res)
    }))
}

metadata <- read.table("../metadata.tsv", header=T, stringsAsFactors=F)
segs <- read.table("../TITAN/titan.seg", header=T)
chrs <- read.table("../data/chr-lengths.tsv", header=T)
chrs$chr.start <- c(0, tail(cumsum(as.numeric(chrs$length)), -1))
segs <- merge(merge(segs, metadata), chrs)

# select dominant subclones only
segs <- merge(segs, aggregate(prevalence~sample, segs, max))
