#!/usr/bin/env Rscript

options(warn=1)
library(ggplot2)
source(file="bedUtils.R")
bed.colnames <- c("chr", "start", "end")

plot.segs <- function (segs, title, width.col) {
    segs <- merge(segs, chrs)
    segs$start <- segs$start + segs$chr.start
    segs$end <- segs$end + segs$chr.start
    segs$count <- segs[,width.col]
    
    ggplot(segs, aes(x=start, y=copy.number)) +
        theme_bw() +
        ylab("copy number") +
        xlab("chromosome") +
        theme(axis.ticks.x=element_blank()) + 
        geom_vline(data=chrs, aes(xintercept=chr.start), color="grey", linetype="dashed") +
        scale_x_continuous(breaks=chrs$midpoint, labels=c(1:22, "X"), limits=c(0, genome.end), expand=c(0, 0)) +
        geom_segment(aes(xend=end, yend=copy.number, size=count)) +
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

bed.cols <- match(bed.colnames, colnames(segs))
data.cols <- setdiff(1:ncol(segs), bed.cols)
segs <- segs[,c(bed.cols, data.cols)]
cat("done\n")

primary <- aggregate(time.point~patient, segs, min)

do.combine <- function (s1, s2) {
    cat("Merging", nrow(s1), "rows with", nrow(s2), "... ")
    if (is.null(s1)) return (s2)
    if (is.null(s2)) return (s1)

    count.cols <- grepl("count", colnames(s1))
    res <- do.call(rbind, lapply(unique(c(s1$copy.number, s2$copy.number)), function (cn) {
        ss1 <- subset(s1, copy.number == cn)
        ss2 <- subset(s2, copy.number == cn)
        if (nrow(ss1) == 0) return (ss2)
        if (nrow(ss2) == 0) return (ss1)
        sub1 <- bedtools("subtract", ss1, ss2)
        sub2 <- bedtools("subtract", ss2, ss1)
        intersect <- bedtools("intersect", ss1, ss2)
        intersect[,count.cols] <- intersect[,count.cols] + 1
        rbind(sub1, sub2, intersect)
    }))
    cat("done\n")
    res
}

do.subtract <- function (primary.segs, meta.segs) {
    sub.segs <- do.call(rbind, lapply(unique(primary.segs$copy.number), function (cn) {
        s1 <- subset(meta.segs, copy.number == cn)
        s2 <- subset(primary.segs, copy.number == cn)
        if (nrow(s1) == 0) return (s2)
        bedtools("subtract", s1, s2)
    }))
    sub.segs$count <- 1
    sub.segs
}

all.segs <- lapply(unique(segs$patient), function (by.patient) {
    pat.segs <- subset(segs, patient == by.patient)
    primary.tpt <- subset(primary, patient == by.patient)$time.point
    metastases <- subset(metadata, patient == by.patient & time.point != primary.tpt)$sample
    
    primary.segs <- subset(pat.segs, time.point == primary.tpt)
    sub.segs <- lapply(metastases, function (m) {
        cat("Subtracting", m, "from primary... ")
        meta.segs <- subset(pat.segs, sample == m)
        res <- do.subtract(primary.segs, meta.segs)
        cat("done\n")
        write.table(res[,c(bed.colnames, "copy.number")], paste0("metastatic-cnv/", m, ".bed"), row.names=F, quote=F, sep="\t", col.names=F)
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

# write data
# TODO: this should not be in this script
sample.segs <- Reduce(do.combine, all.segs)
sample.segs <- sample.segs[order(sample.segs$chr, sample.segs$start),]
sample.segs <- sample.segs[,c(bed.colnames, "copy.number", "count", "count.patient")]
write.table(sample.segs, file="combined-segs.bed", row.names=F, quote=F, sep="\t", col.names=F)

pdf(file.path("metastatic-cnv", "all-samples.pdf"), width=12, height=5)
print(plot.segs(sample.segs, "all by sample", "count"))
dev.off()

pdf(file.path("metastatic-cnv", "all-patients.pdf"), width=12, height=5)
print(plot.segs(sample.segs, "all by patient", "count.patient"))
dev.off()
