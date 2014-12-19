#!/usr/bin/env Rscript

library(ggplot2)
library(gtable)
library(gridExtra)

source(file="bedUtils.R")
chrs <- c(1:22, "X")

# read Diep data
diep <- read.table("../data/diep.tsv", header=T, stringsAsFactors=F)
diep$chr <- factor(sub("[pq].*", "", diep$band), levels=chrs)

# read band locations
bands <- read.table("../data/bands.bed", 
                    col.names=c("chr", "start", "end", "band"))
bands$band <- paste0(bands$chr, bands$band)
bands <- subset(bands, chr %in% c(1:22, "X"))
bands$chr <- factor(bands$chr, levels=c(1:22, "X"))
chr.data <- aggregate(end~chr, bands, max)
chr.data$chr.start <- c(0, head(cumsum(as.numeric(chr.data$end)), -1))
chr.data$chr.midpoint <- chr.data$chr.start + chr.data$end/2
genome.end <- sum(as.numeric(chr.data$end))
chr.data <- chr.data[,c("chr", "chr.start", "chr.midpoint")]
bands <- merge(bands, chr.data)

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
do.aggregate <- function (segs, group.col) {
    agg <- do.call(rbind, as.list(by(segs, segs[,group.col], function (x) {
        bedtools("intersect", bands, x, args="-u")
    })))
    agg <- setNames(aggregate(chr~band, agg, length), c("band", "count"))
    agg$count <- agg$count * 100 / length(unique(segs[,group.col]))
    agg
}

loss.segs <- subset(segs, copy.number < 2)
gain.segs <- subset(segs, copy.number > 2)
loss.sample <- do.aggregate(loss.segs, "sample")
loss.sample$data <- "QCROC samples"
gain.sample <- do.aggregate(gain.segs, "sample")
gain.sample$data <- "QCROC samples"

# then by patient
loss.patient <- do.aggregate(loss.segs, "patient")
loss.patient$data <- "QCROC patients"
gain.patient <- do.aggregate(gain.segs, "patient")
gain.patient$data <- "QCROC patients"

# merge together with Diep data
diep.gain.primary <- setNames(diep[,c("band", "gains")], c("band", "count"))
diep.gain.primary$data <- "Diep et al. primary"
diep.loss.primary <- setNames(diep[,c("band", "losses")], c("band", "count"))
diep.loss.primary$data <- "Diep et al. primary"
diep.gain.meta <- setNames(diep[,c("band", "gains.liver")], c("band", "count"))
diep.gain.meta$data <- "Diep et al. metastases"
diep.loss.meta <- setNames(diep[,c("band", "losses.liver")], c("band", "count"))
diep.loss.meta$data <- "Diep et al. metastases"

do.plot <- function (plot.data, title) {
    plot.data <- merge(plot.data, bands, by=c("band"))
    plot.data$start <- plot.data$start + plot.data$chr.start
    plot.data$end <- plot.data$end + plot.data$chr.start

    ggplot(plot.data, aes(x=start, y=count, color=data)) + 
        geom_step() +
        theme_bw() +
        ylab("copy number") +
        xlab("chromosome") +
        theme(axis.ticks.x=element_blank(), legend.position="bottom") +
        geom_vline(data=chr.data, aes(xintercept=chr.start), color="grey", linetype="dashed") +
        scale_x_continuous(breaks=chr.data$chr.midpoint, labels=c(1:22, "X"), limits=c(0, genome.end), expand=c(0, 0)) +
        scale_y_continuous(limits=c(0, 80)) +
        ggtitle(title)
}

pdf("diep/patient-primary.pdf", width=12, height=5)
print(do.plot(rbind(gain.patient, diep.gain.primary), "Gains by patient"))
print(do.plot(rbind(loss.patient, diep.loss.primary), "Losses by patient"))
dev.off()

pdf("diep/patient-metastasis.pdf", width=12, height=5)
print(do.plot(rbind(gain.patient, diep.gain.meta), "Gains by patient"))
print(do.plot(rbind(loss.patient, diep.loss.meta), "Losses by patient"))
dev.off()

pdf("diep/sample-primary.pdf", width=12, height=5)
print(do.plot(rbind(gain.sample, diep.gain.primary), "Gains by sample"))
print(do.plot(rbind(loss.sample, diep.loss.primary), "Losses by sample"))
dev.off()

pdf("diep/sample-metastasis.pdf", width=12, height=5)
print(do.plot(rbind(gain.sample, diep.gain.meta), "Gains by sample"))
print(do.plot(rbind(loss.sample, diep.loss.meta), "Losses by sample"))
dev.off()
