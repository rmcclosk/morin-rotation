#!/usr/bin/env Rscript

options(warn=1)
library(ggplot2)
source(file="bedUtils.R")

plot.segs <- function (segs, title, width.col) {
    segs <- merge(segs, chrs)
    segs$start <- segs$start + segs$chr.start
    segs$end <- segs$end + segs$chr.start
    segs$width <- segs[,width.col]
    
    ggplot(segs, aes(x=start, y=copy.number)) +
        theme_bw() +
        ylab("copy number") +
        xlab("chromosome") +
        theme(axis.ticks.x=element_blank()) + #, legend.position="top", legend.box="horizontal") +
        geom_vline(data=chrs, aes(xintercept=chr.start), color="grey", linetype="dashed") +
        scale_x_continuous(breaks=chrs$midpoint, labels=c(1:22, "X"), limits=c(0, genome.end), expand=c(0, 0)) +
        geom_segment(aes(xend=end, yend=copy.number, size=width)) +
        ggtitle(title)
}

metadata <- read.table("../metadata.tsv", header=T, stringsAsFactors=F)

chrs <- read.table("../data/chr-lengths.tsv", header=T)
chrs <- subset(chrs, chr %in% c(1:22, "X"))
chrs$chr.start <- c(0, head(cumsum(as.numeric(chrs$length)), -1))
chrs$midpoint <- chrs$chr.start+chrs$length/2
genome.end <- max(chrs$chr.start+chrs$length)

cat("Reading segments... ")
segs <- read.table("../TITAN/titan.seg", header=T)
segs <- subset(segs, copy.number != 2)
segs <- merge(segs, metadata)
colnames(segs)[colnames(segs) == "chrom"] <- "chr"
n.unique <- function (x) length(unique(x))
keeps <- subset(aggregate(sample~patient, segs, n.unique), sample > 1)
segs <- droplevels(subset(segs, patient %in% keeps$patient))

bed.cols <- match(c("chr", "start", "end"), colnames(segs))
data.cols <- setdiff(1:ncol(segs), bed.cols)
segs <- segs[,c(bed.cols, data.cols)]
cat("done\n")

primary <- aggregate(time.point~patient, segs, min)

do.combine <- function (s1, s2) {
    cat("Merging", nrow(s1), "rows with", nrow(s2), "... ")
    if (is.null(s1)) return (s2)
    if (is.null(s2)) return (s1)
    sub1 <- bedtools("subtract", s1, s2)
    sub2 <- bedtools("subtract", s2, s1)
    int1 <- bedtools("intersect", s1, s2)
    int2 <- bedtools("intersect", s2, s1)
    intersect <- rbind(int1, int2)
    intersect <- intersect[!duplicated(intersect),]
    count.cols <- grepl("count", colnames(intersect))
    intersect[,count.cols] <- intersect[,count.cols] + 1
    cat("done\n")
    rbind(sub1, sub2, intersect)
}

all.segs <- lapply(unique(segs$patient), function (by.patient) {
    pat.segs <- subset(segs, patient == by.patient)
    primary.tpt <- subset(primary, patient == by.patient)$time.point
    metastases <- subset(metadata, patient == by.patient & time.point != primary.tpt)$sample
    
    primary.segs <- subset(pat.segs, time.point == primary.tpt)
    sub.segs <- lapply(metastases, function (m) {
        cat("Subtracting", m, "from primary... ")
        meta.segs <- subset(pat.segs, sample == m)
        res <- bedtools("subtract", meta.segs, primary.segs)
        res$count <- 1
        cat("done\n")
        res
    })
    comb.segs <- Reduce(do.combine, sub.segs)
    comb.segs <- comb.segs[order(comb.segs$chr, comb.segs$start),]
    pdf(file.path("metastatic-cnv", paste0(by.patient, ".pdf")), width=12, height=5)
    print(plot.segs(comb.segs, by.patient, "count"))
    dev.off()
    comb.segs$count.patient <- 1
    comb.segs
})

sample.segs <- Reduce(do.combine, all.segs)
sample.segs <- sample.segs[order(sample.segs$chr, sample.segs$start),]
pdf(file.path("metastatic-cnv", "all-samples.pdf"), width=12, height=5)
print(plot.segs(sample.segs, "all by sample", "count"))
dev.off()

pdf(file.path("metastatic-cnv", "all-patients.pdf"), width=12, height=5)
print(plot.segs(sample.segs, "all by patient", "count.patient"))
dev.off()
