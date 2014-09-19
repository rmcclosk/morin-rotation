#!/usr/bin/env Rscript

source(file="settings.conf")

sample.data <- read.csv("BC Genome Profiling-HQC.csv", header=T, fill=T, stringsAsFactors=F)
sample.data$purity <- (sample.data$X..tumor/100) * 
                      ((sample.data$X..viable.neoplastic.cells/100) +
                       (sample.data$X..necrosis/100))
sample.data$purity <- ifelse(is.na(sample.data$purity), 
                             sample.data$X..Tumor/100, 
                             sample.data$purity)
sample.data <- sample.data[,c("Sample.ID", "purity")]
colnames(sample.data) <- c("sample", "report.purity")
sample.data <- sample.data[grepl("^01", sample.data$sample),]
sample.data$patient.id <- paste0("01-", sapply(strsplit(sample.data$sample, "-"), "[[", 3))
sample.data$report.purity <- sub(".0$", "", floor((sample.data$report.purity)/0.05)*0.05)

seg.data <- read.table("annotations.dat", header=T)

apply(sample.data, 1, function (row) {
    pdf(file.path("plots"), paste0(sample, ".pdf")), width=9, height=9)
    par(mfrow=c(2, 2), xaxs="i")

    sapply(c(1:22, "X", "Y"), function (plot.chr) {
        plot.seg <- subset(seg.data, sample==row["sample"] & 
                                     purity==row["report.purity"] &
                                     chr==plot.chr)
        genome.fun <- stepfun(plot.seg$start, c(0, plot.seg$genome.copy))
        exome.fun <- stepfun(plot.seg$start, c(0, plot.seg$exome.copy))
    
        plot(genome.fun, col="red", do.points=F,
             main="",
             xlab=paste("Position on chromosome", plot.chr),
             ylab="Copy number", lwd=2, xlim=c(min(segments$start), max(segments$end)))
        plot(exome.fun, col="blue", do.points=F, add=T, lwd=2)

        xmax <- par("usr")[2]
        ymax <- par("usr")[4]
        par(xpd=T)
        legend(xmax, ymax, legend=c("HMMcopy", "ExomeCNV"), col=c("red", "blue"), lty=1,
               xjust=1, yjust=0)
        par(xpd=F)
    })
    dev.off()
})
