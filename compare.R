#!/home/rmccloskey/bin/Rscript

source(file="settings.conf")

patient.id <- "01-007"
sample <- "01-2-007-4"

exome.dir <- file.path(WORK_DIR, "03_cnv", patient.id, sample)
genome.dir <- file.path(WORK_DIR, "05_hmmcopy", patient.id, sample)

admixture.rate <- 0.3
exome.file <- file.path(exome.dir, 
                        paste0(admixture.rate, ".segment.copynumber.txt"))
exome.seg <- read.table(exome.file)
colnames(exome.seg) <- c("chr", "start", "end", "copy.number")
exome.seg$chr <- factor(sub("chr", "", exome.seg$chr))

genome.file <- file.path(genome.dir, "1409141605_segments.dat")
genome.seg <- read.table(genome.file, header=T)
genome.seg$copy.number <- genome.seg$state - 1

quit()

pdf(paste0(patient.id, "_", sample, ".pdf"), width=9, height=9)
par(mfrow=c(2, 2), xaxs="i")

sapply(c(1:22, "X", "Y"), function (plot.chr) {

    genome.chr.seg <- subset(genome.seg, chr==plot.chr)
    exome.chr.seg <- subset(exome.seg, chr==plot.chr)
    
    if (plot.chr == 20) {
        print(genome.chr.seg)
        print(exome.chr.seg)
    }

    starts <- sort(unique(c(genome.chr.seg$start, genome.chr.seg$end+1, exome.chr.seg$start, exome.chr.seg$end+1)))
    starts <- starts[starts >= min(exome.chr.seg$start) & starts < max(exome.chr.seg$end)]
    segments <- data.frame(start=head(starts, -1), end=tail(starts, -1)-1)
    
    segments <- cbind(segments, t(apply(segments, 1, function (s) {
        c(genome.copy=subset(genome.chr.seg, start <= s[1] & end >= s[2], select=copy.number)[1,1],
          exome.copy=subset(exome.chr.seg, start <= s[1] & end >= s[2], select=copy.number)[1,1])
    })))
    segments <- segments[!is.na(segments$exome.copy),]
    segments$length <- segments$end-segments$start+1
    
    genome.fun <- stepfun(segments$start, c(0, segments$genome.copy))
    exome.fun <- stepfun(segments$start, c(0, segments$exome.copy))
    
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
