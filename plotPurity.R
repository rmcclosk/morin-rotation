#!/usr/bin/env Rscript

sample.data <- read.csv("BC Genome Profiling-HQC.csv", header=T, fill=T)
sample.data$purity <- (sample.data$X..tumor/100) * 
                      ((sample.data$X..viable.neoplastic.cells/100) +
                       (sample.data$X..necrosis/100))
sample.data$purity <- ifelse(is.na(sample.data$purity), 
                             sample.data$X..Tumor/100, 
                             sample.data$purity)
sample.data <- sample.data[,c("Sample.ID", "purity")]
colnames(sample.data) <- c("sample", "report.purity")

d <- read.table("annotations.dat", header=T)
d$length <- d$end-d$start+1
d$exome.cnv <- d$exome.copy != 2
d$genome.cnv <- d$genome.copy != 2
d$exome.cnv.length <- d$exome.cnv*d$length
d$genome.cnv.length <- d$genome.cnv*d$length

total.length <- aggregate(length~sample+purity, d, sum)
exome.cnv.length <- aggregate(exome.cnv.length~sample+purity, d, sum)
genome.cnv.length <- aggregate(genome.cnv.length~sample+purity, d, sum)

agg.d <- merge(total.length, exome.cnv.length)
agg.d <- merge(agg.d, genome.cnv.length)
agg.d <- merge(agg.d, sample.data)
agg.d$prop.exome.cnv <- agg.d$exome.cnv.length/agg.d$length
agg.d$prop.genome.cnv <- agg.d$genome.cnv.length/agg.d$length

pdf("cnv.pdf")
par(mfrow=c(2,2))
. <- by(agg.d, agg.d$sample, function (x) {
    plot(x$purity*100, x$prop.exome.cnv*100, type="l",
         main=paste("Sample", x[1,"sample"]),
         xlab="% estimated tumor purity",
         ylab="% exome with CNV")
    abline(h=mean(x$prop.genome.cnv)*100, lty=2, col="blue")
    abline(v=x$report.purity*100, lty=2, col="red")
    legend("topright", col=c("blue", "red"), lty=2, 
           legend=c("HMMcopy", "pathologist"), bg="white")
})
dev.off()
