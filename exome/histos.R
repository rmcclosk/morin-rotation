#!/usr/bin/env Rscript

source(file="../settings.conf")
coverage.dir <- file.path(WORK_DIR, "02_coverage")

read.cnv <- function (cnv.file) {
    cat("Reading ", cnv.file, "... ")
    cnv <- read.table(cnv.file, col.names=c("chrom", "start", "end", "logR"))
    cnv$ratio <- 2^cnv$logR
    cnv$chrom <- sub("chr", "", cnv$chrom)
    cnv$sample <- basename(dirname(cnv.file))
    cat("done\n")
    cnv
}

find.start <- function (chrom, pos, intervals) {
    subset(intervals[[chrom]], start <= pos & end >= pos)[1,"start"]
}

parse.segments <- function (coverage.file) {
    cat("Reading", coverage.file, "... ")
    coverage <- read.table(coverage.file, header=T, sep="\t", stringsAsFactors=F)[,1:3]
    target <- strsplit(coverage$Target, ":|-")
    coverage$chr <- sapply(target, "[[", 1)
    coverage$start <- as.numeric(sapply(target, "[[", 2))
    coverage$end <- as.numeric(sapply(target, "[[", 3))
    coverage$length <- coverage$end-coverage$start+1
    coverage$total_coverage <- as.numeric(coverage$total_coverage)
    avg <- sum(coverage$total_coverage)/sum(coverage$length)
    coverage$proportion <- coverage$total_coverage/avg
    cat("done\n")
    coverage
}

s <- read.table(METADATA, header=T, stringsAsFactors=F)
s <- s[!is.na(s$normal.sample),]

cat("Reading SNVs... ")
d <- read.table("snv.tsv", header=T, sep="\t", fill=T)
d <- droplevels(d[d$chrom %in% c(1:22, "X", "Y"),])
d$chrom <- factor(d$chrom, levels=c(1:22, "X", "Y"))
d$vaf <- d$alt.count/d$depth
cat("done\n")

a <- read.table("absolute-exomecnv-0.txt", header=T, sep="\t", fill=T)
colnames(a)[colnames(a)=="purity"] <- "purity.absolute.1"
a <- a[,c("sample", "purity.absolute.1")]
d <- merge(d, a)

a <- read.table("absolute-exomecnv-bootstrap.txt", header=T, sep="\t", fill=T)
colnames(a)[colnames(a)=="purity"] <- "purity.absolute.2"
a <- a[,c("sample", "purity.absolute.2")]
d <- merge(d, a)

a <- read.table("titan.txt", col.names=c("sample", "purity.titan"))
d <- merge(d, a)

d <- d[!is.na(d$vaf),]
d <- merge(d, s, by.x=c("sample"), by.y=c("tumor.sample"))

pdf("histos.pdf")
. <- by(d, d$sample, function (ss) {
    mu <- mean(ss$vaf)
    sigma <- sd(ss$vaf)
    ss <- subset(ss, abs((vaf-mu)/sigma) <= 3) # trim outliers
    ss1 <- subset(ss, vaf >= 0.1)
    if (nrow(ss) == 0 | nrow(ss1) == 0) return(NULL)
    
    dens <- density(ss$vaf)
    dens1 <- density(ss1$vaf)
    ymax <- max(c(dens$y, dens1$y))
    cat(levels(d$sample)[ss[1, "sample"]], 2*dens$x[which.max(dens1$y)], "\n")
    
    plot(dens, col="blue", main=ss[1,"sample"], xlab="variant allelic fraction",
         ylim=c(0, ymax))
    lines(dens1, col="purple")

    vlines <- c(ss[1,"purity"]/200, ss[1,"purity.absolute.1"]/2,
                ss[1,"purity.absolute.2"]/2, ss[1,"purity.titan"]/2)
    vlines[duplicated(vlines)] <- vlines[duplicated(vlines)] + 0.01
    abline(v=vlines[1], lwd=2, lty=2, col="red")
    abline(v=vlines[2], lwd=2, lty=2, col="forestgreen")
    abline(v=vlines[3], lwd=2, lty=2, col="cyan")
    abline(v=vlines[4], lwd=2, lty=2, col="darkgoldenrod")
    legend("topright", 
           col=c("blue", "purple", "red", "forestgreen", "cyan", "darkgoldenrod"),
           lty=c(1,1,2,2,2,2),
           lwd=2,
           legend=c("all", "VAF >= 0.1", "pathologist", "absolute basic", "absolute bootstrap", "titan"),
           bg="white")
})
dev.off()
