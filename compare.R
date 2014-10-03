#!/usr/bin/env Rscript

source(file="settings.conf")

vaf.dir <- file.path(WORK_DIR, "11_vaf")

rep.ends <- function (x) {
    c(head(x, 1), x, tail(x, 1))
}

seg.data <- read.table("annotations.dat", header=T)

# remove segments likely to cause error
seg.data$arm <- substr(seg.data$locus, 1, 1)
seg.data <- subset(seg.data, (! chr %in% c(13, 14, 15, 21, 22)) | arm != "p")
seg.data <- subset(seg.data, chr != "Y")

# use only the purity closest to the pathologist's estimate
purity.to.use <- do.call(rbind, by(seg.data, seg.data$sample, function (x) {
    x[which.min(abs(x$purity-x$report.purity)), c("sample", "purity")]
}, simplify=F))
seg.data <- merge(seg.data, purity.to.use)

seg.data$genome.copy <- as.integer(seg.data$genome.copy)
seg.data[seg.data$genome.copy > 3,]$genome.copy <- 3
seg.data[seg.data$genome.copy < 1,]$genome.copy <- 1

by(seg.data, seg.data$sample, function (by.seg.data) {
    vaf.path <- file.path(vaf.dir, paste0(by.seg.data[1,"sample"], ".dat"))
    vaf.data <- read.table(vaf.path, header=T)

    pdf(file.path("plots", paste0(by.seg.data[1,"sample"], ".pdf")), width=9, height=9)
    par(mfrow=c(2, 2), oma=c(0,0,0,1))

    sapply(c(1:22, "X"), function (plot.chr) {
        by.vaf <- subset(vaf.data, chr==plot.chr)
        plot.seg <- subset(by.seg.data, chr==plot.chr)
        plot.seg <- plot.seg[order(plot.seg$start),]
        bounds <- c(plot.seg$start, tail(plot.seg$end, 1))
        genome.fun <- stepfun(bounds, rep.ends(plot.seg$genome.copy)+0.02)
        exome.fun <- stepfun(bounds, rep.ends(plot.seg$exome.copy)-0.02)
    
        plot(genome.fun, col="red", do.points=F,
             main="",
             xlab=paste("Position on chromosome", plot.chr),
             ylab="Copy number", lwd=2, xlim=c(min(bounds), max(bounds)),
             xaxs="i", ylim=c(0.5,3.5), yaxs="i", yaxt="n")
        plot(exome.fun, col="blue", do.points=F, add=T, lwd=2, yaxt="n")
        axis(2, at=1:3, labels=c("loss", "neutral", "gain"))
        points(by.vaf$pos, 2*by.vaf$vaf+1)
        abline(h=mean(2*by.vaf$vaf+1), lty=2, lwd=2)
        axis(4, at=1:3, labels=c(0, 0.5, 1))
        mtext("variant allelic fraction", side=4)

        xmax <- par("usr")[2]
        ymax <- par("usr")[4]
        par(xpd=T)
        legend(xmax, ymax, legend=c("HMMcopy", "ExomeCNV", "SNP"), 
               col=c("red", "blue", "black"), lty=c(1, 1, NA), 
               pch=c(NA, NA, 1), xjust=1, yjust=0)
        par(xpd=F)
    })
    dev.off()
})
