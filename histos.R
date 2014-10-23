#!/usr/bin/env Rscript

source(file="settings.conf")
WORK_DIR="/extscratch/morinlab/shared/rmccloskey/colorectal_exome"
cnv.dir <- file.path(WORK_DIR, "03_cnv")
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
    coverage <- read.table(coverage.file, header=T, stringsAsFactors=F)
    target <- strsplit(coverage$Target, ":|-")
    coverage$chr <- sapply(target, "[[", 1)
    coverage$start <- as.numeric(sapply(target, "[[", 2))
    coverage$end <- as.numeric(sapply(target, "[[", 3))
    coverage$length <- coverage$end-coverage$start+1
    coverage$total_coverage <- as.numeric(coverage$total_coverage)
    avg <- sum(coverage$total_coverage*coverage$length)/sum(coverage$length)
    coverage$proportion <- coverage$total_coverage/avg
    cat("done\n")
    coverage
}

s <- read.table(METADATA, header=T, stringsAsFactors=F)
# DEBUG
s <- head(s, 1)
samples <- unique(c(s$tumor.sample, s$normal.sample))
get.coverage.file <- function (s) {
    file.path(coverage.dir, paste0(s, ".sample_interval_summary"))
}
coverage.files <- sapply(samples, get.coverage.file)
coverage <- lapply(coverage.files, parse.segments)
names(coverage) <- samples

apply(s, 1, function (row) {
    tumor.prop <- coverage[[row["tumor.sample"]]]$proportion
    normal.prop <- coverage[[row["normal.sample"]]]$proportion
    tumor.prop/normal.prop
})
quit()
                         
# read coverage data
norm.coverage <- parse.segments(read.table(opt$normcoverage, header=T))
tum.coverage <- parse.segments(read.table(opt$tumcoverage, header=T))

# scale so average coverage is the same
norm.ave <- sum(norm.coverage$total_coverage) / sum(norm.coverage$length)
tum.ave <- sum(tum.coverage$total_coverage) / sum(tum.coverage$length)
tum.coverage$average_coverage <- (norm.ave/tum.ave) * tum.coverage$average_coverage
logR <- log(tum.coverage$average_coverage) - log(norm.coverage$average_coverage)

quit()

cat("Reading intervals list... ")
intervals <- read.table(INTERVAL_LIST, stringsAsFactors=F)
split <- strsplit(intervals$V1, "[:-]")
intervals$chrom <- factor(sapply(split, "[[", 1), levels=c(1:22, "X", "Y"))
intervals$start <- as.integer(sapply(split, "[[", 2))
intervals$end <- as.integer(sapply(split, "[[", 3))
intervals <- intervals[,c("chrom", "start", "end")]
intervals <- as.list(by(intervals, intervals$chrom, identity))
cat("done\n")
names(intervals) <- c(1:22, "X", "Y")

cat("Reading SNVs... ")
d <- read.table("freqs.dat", header=T, sep="\t", fill=T)
cnv.files <- sapply(levels(d$sample), function (s) {
    file.path(cnv.dir, s, "0.exon.lrr.txt")
})
d <- droplevels(d[file.exists(cnv.files),])
d <- droplevels(d[d$chrom %in% c(1:22, "X", "Y"),])
d$chrom <- factor(d$chrom, levels=c(1:22, "X", "Y"))
cat("done\n")

cat("Finding exon start positions... ")
d$start <- mapply(find.start, as.character(d$chrom), d$pos, MoreArgs=list(intervals))
d$vaf <- d$alt.count/d$depth
cat("done\n")

cat("Reading metadata... ")
s <- read.table(METADATA, header=T)
d <- merge(d, s, by.x=c("patient", "sample"), by.y=c("patient", "tumor.sample"))
cat("done\n")

cnv <- do.call(rbind, lapply(cnv.files[file.exists(cnv.files)], read.cnv))
d <- merge(d, cnv)
d <- d[!is.na(d$vaf) & !is.na(d$ratio),]

pdf("test.pdf")
by(d, d$sample, function (ss) {
    mu <- mean(ss$ratio*ss$vaf)
    sigma <- sd(ss$ratio*ss$vaf)
    ss <- subset(ss, abs((ratio*vaf-mu)/sigma) <= 3) # trim outliers
    ss1 <- subset(ss, depth >= 10)
    ss2 <- subset(ss1, vaf >= 0.1)
    
    dens <- density(ss$vaf*ss$ratio)
    dens1 <- density(ss1$vaf*ss1$ratio)
    dens2 <- density(ss2$vaf*ss2$ratio)
    
    ymax <- max(c(dens$y, dens1$y, dens2$y))
    
    plot(density(2*ss$vaf*ss$ratio), col="blue",
         main=ss[1,"sample"],
         xlab="2*copy ratio*VAF",
         ylim=c(0, ymax))
    lines(density(2*ss1$vaf*ss1$ratio), col="forestgreen")
    lines(density(2*ss2$vaf*ss2$ratio), col="purple")
    abline(v=ss[1,"purity"]/100, lwd=2, lty=2, col="red")
    legend("topright", 
           col=c("blue", "forestgreen", "purple", "red"),
           lty=c(1,1,1,2),
           legend=c("all", "read depth >= 10", "read depth >= 10 and VAF >= 0.1", "pathologist"),
           bg="white")
})
dev.off()
