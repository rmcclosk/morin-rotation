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
    output <- system(cmd, intern=T)
    if (length(output) == 0) return (NULL)
    res <- setNames(read.table(textConnection(output)), colnames(b1))
    res$chr <- factor(res$chr, levels=c(1:22, "X"))
    res
}

plot.segs <- function (segs, title) {
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
head(segs)

primary <- aggregate(time.point~patient, segs, min)

do.combine <- function (s1, s2) {
    cat("Merging", nrow(s1), "rows with", nrow(s2), "... ")
    if (is.null(s1)) return (s2)
    if (is.null(s2)) return (s1)
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
    comb.segs <- comb.segs[order(comb.segs$chr, comb.segs$start),]
    pdf(file.path("metastatic-cnv", paste0(by.patient, ".pdf")), width=12, height=5)
    print(plot.segs(comb.segs, by.patient))
    dev.off()
    comb.segs
})

all.segs <- Reduce(do.combine, all.segs)
pdf(file.path("metastatic-cnv", "all.pdf"), width=12, height=5)
print(plot.segs(all.segs, "all"))
dev.off()
