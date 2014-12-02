#!/usr/bin/env Rscript

options(warn=1)
#library(ggplot2)

all.sample.groups <- function (samples, time.points) {
    if (length(unique(time.points)) == 1)
        return(matrix(samples, ncol=1))
    min.time.point <- min(time.points)
    first.idx <- which(time.points == min.time.point)
    first <- samples[first.idx]
    res <- all.sample.groups(samples[-first.idx], time.points[-first.idx])
    do.call(rbind, lapply(first, function (s) {
        cbind(rep(s, nrow(res)), res)
    }))
}

plot.segs.trappings <- function (seg.plot) {
    seg.plot +
        theme_bw() +
        ylab("copy number") +
        xlab("chromosome") +
        theme(axis.ticks.x=element_blank()) + #, legend.position="top", legend.box="horizontal") +
        geom_vline(data=chrs, aes(xintercept=chr.start), color="grey", linetype="dashed") +
        scale_x_continuous(breaks=chrs$midpoint, labels=c(1:22, "X"), limits=c(0, genome.end), expand=c(0, 0))
}

plot.segs <- function (segs) {
    p <- ggplot(segs, aes(x=start, y=copy.number, color=`time point`))
    p <- plot.segs.trappings(p)
    p + geom_segment(aes(xend=end, yend=copy.number, size=prevalence)) +
        ggtitle(segs[1, "patient"])
}

# read in segments
metadata <- read.table("../metadata.tsv", header=T, stringsAsFactors=F)
segs <- read.table("../TITAN/titan.seg", header=T)
chrs <- read.table("../data/chr-lengths.tsv", header=T)
chrs <- subset(chrs, chr %in% c(1:22, "X"))
chrs$chr.start <- c(0, head(cumsum(as.numeric(chrs$length)), -1))
chrs$midpoint <- chrs$chr.start+chrs$length/2
genome.end <- max(chrs$chr.start+chrs$length)
segs <- merge(merge(segs, metadata), chrs)
segs$start <- segs$start + segs$chr.start
segs$end <- segs$end + segs$chr.start

# read in genes
gene.files <- list.files("../genes/by-sample", full.names=T)
gene.data <- lapply(gene.files, read.table)
names(gene.data) <- sapply(strsplit(basename(gene.files), ".", fixed=T), "[[", 1)

# fill in copy number 2 segments for rest of prevalence
het <- subset(segs, copy.number != 2)
het$copy.number <- 2
het$prevalence <- 1-het$prevalence
segs <- rbind(segs, het)
segs <- segs[order(segs$patient, segs$sample, segs$start),]

# make plots
by(segs, segs$patient, function (pat.segs) {
    by.patient <- pat.segs[1, "patient"]
    samples <- as.character(unique(pat.segs$sample))
    time.points <- unique(pat.segs$time.point)
    n.tpt <- length(time.points)
    dodge <- 0.25
    adj <- data.frame(time.point=time.points, adj=seq(-dodge*(n.tpt-1)/2, dodge*n.tpt/2, dodge))
    pat.segs <- merge(pat.segs, adj)
    pat.segs$copy.number <- pat.segs$copy.number + pat.segs$adj
    pat.segs$`time point` <- factor(pat.segs$time.point)
    groups <- all.sample.groups(samples, time.points)
    pdf(file.path("titan-segs", paste0(by.patient, ".pdf")), width=12, height=5)
    apply(groups, 1, function (sample.group) {
        plot.segs(subset(pat.segs, sample %in% sample.group))
    })
    dev.off()
})
