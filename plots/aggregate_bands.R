#!/usr/bin/env Rscript

source(file="bedUtils.R")

chrs <- c(1:22, "X")

# read Diep data
diep <- read.table("../data/diep.tsv", header=T, stringsAsFactors=F)
diep$chr <- factor(sub("[pq].*", "", diep$locus), levels=chrs)
diep$band <- sub("^[0-9X]+", "", diep$locus)

# read band locations
bands <- read.table("../data/bands.bed", 
                    col.names=c("chr", "start", "end", "band"))
bands$band <- paste0(bands$chr, bands$band)

# read segments
segs <- read.delim("../TITAN/titan.seg")
colnames(segs)[colnames(segs) == "chrom"] <- "chr"
segs <- subset(segs, copy.number != 2 & end-start > 0)

# read metadata
metadata <- read.delim("../metadata_colorectal.csv")
metadata <- metadata[,c("patient", "tumor.sample")]
metadata <- setNames(metadata, c("patient", "sample"))
segs <- merge(segs, metadata)

# select clonal events only
clonal.prev <- aggregate(prevalence~sample, segs, max, na.rm=T)
segs <- merge(segs, clonal.prev, by=c("sample"), suffixes=c("", ".clonal"))
segs <- subset(segs, prevalence == prevalence.clonal)
seg.cols <- match(c("chr", "start", "end"), colnames(segs))
segs <- segs[,c(seg.cols, setdiff(1:ncol(segs), seg.cols))]

# aggregate loss and gain counts by sample
loss.segs <- subset(segs, copy.number > 2)
gain.segs <- subset(segs, copy.number < 2)

loss.bands.sample <- bedtools("intersect", bands, loss.segs)
gain.bands.sample <- bedtools("intersect", bands, gain.segs)

# then by patient
