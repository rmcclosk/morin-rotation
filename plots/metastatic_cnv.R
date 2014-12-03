#!/usr/bin/env Rscript

options(warn=1)
library(ggplot2)

write.bed <- function (bed) {
    f <- tempfile()
    write.table(bed, f, col.names=F, row.names=F, quote=F, sep="\t")
    f
}

bedtools <- function (cmd, b1, b2, args=NULL) {
    f1 <- write.bed(b1)
    f2 <- write.bed(b2)
    cmd <- paste("bedtools", cmd, args, "-a", f1, "-b", f2)
    res <- read.table(textConnection(system(cmd, intern=T)))
    setNames(res, colnames(b1))
}

plot.segs <- function (segs) {
    segs <- merge(segs, chrs)
    segs$start <- segs$start + segs$chr.start
    segs$end <- segs$end + segs$chr.start
    
    ggplot(segs, aes(x=start, y=copy.number)) +
        theme_bw() +
        ylab("copy number") +
        xlab("chromosome") +
        theme(axis.ticks.x=element_blank()) + #, legend.position="top", legend.box="horizontal") +
        geom_vline(data=chrs, aes(xintercept=chr.start), color="grey", linetype="dashed") +
        scale_x_continuous(breaks=chrs$midpoint, labels=c(1:22, "X"), limits=c(0, genome.end), expand=c(0, 0)) +
        geom_segment(aes(xend=end, yend=copy.number, size=count)) +
        ggtitle(segs[1, "patient"])
}

metadata <- read.table("../metadata.tsv", header=T, stringsAsFactors=F)
counts <- subset(aggregate(sample~patient, metadata, length), sample > 1)
metadata <- droplevels(subset(metadata, patient %in% counts$patient))

chrs <- read.table("../data/chr-lengths.tsv", header=T)
chrs <- subset(chrs, chr %in% c(1:22, "X"))
chrs$chr.start <- c(0, head(cumsum(as.numeric(chrs$length)), -1))
chrs$midpoint <- chrs$chr.start+chrs$length/2
genome.end <- max(chrs$chr.start+chrs$length)

cat("Reading segments... ")
segs <- read.table("../TITAN/titan.seg", header=T)
colnames(segs)[colnames(segs) == "chrom"] <- "chr"

segs <- merge(segs, metadata)
segs[is.na(segs$prevalence),"prevalence"] <- 1
segs <- subset(segs, copy.number != 2)
bed.cols <- match(c("chr", "start", "end"), colnames(segs))
data.cols <- setdiff(1:ncol(segs), bed.cols)
segs <- segs[,c(bed.cols, data.cols)]
cat("done\n")

primary <- aggregate(time.point~patient, metadata, min)

do.combine <- function (s1, s2) {
    cat("Merging... ")
    if (nrow(s1) == 0) return (s2)
    if (nrow(s2) == 0) return (s1)
    sub1 <- bedtools("subtract", s1, s2)
    sub2 <- bedtools("subtract", s2, s1)
    intersect <- bedtools("intersect", s1, s2)
    intersect$count <- intersect$count + 1
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
    print(head(comb.segs))
    comb.segs <- comb.segs[order(comb.segs$chr, comb.segs$start),]
    pdf(file.path("metastatic-cnv", paste0(by.patient, ".pdf")), width=12, height=5)
    plot.segs(comb.segs)
    dev.off()
})

all.segs <- Reduce(do.combine, all.segs)
pdf(file.path("metastatic-cnv", "all.pdf"), width=12, height=5)
plot.segs(all.segs)
dev.off()
