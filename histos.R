#!/usr/bin/env Rscript

source(file="settings.conf")
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
#samples <- unique(c(s$tumor.sample, s$normal.sample))
#get.coverage.file <- function (s) {
#    file.path(coverage.dir, paste0(s, ".sample_interval_summary"))
#}
#coverage.files <- sapply(samples, get.coverage.file)
#coverage <- lapply(coverage.files, parse.segments)
#names(coverage) <- samples

#cnv <- do.call(rbind, lapply(1:nrow(s), function (i) {
#    tumor.coverage <- coverage[[s[i,"tumor.sample"]]]
#    normal.prop <- coverage[[s[i,"normal.sample"]]]$proportion
#    tumor.coverage$sample <- s[i,"tumor.sample"]
#    tumor.coverage$ratio <- tumor.coverage$proportion/normal.prop
#    tumor.coverage
#}))
#                         
#cat("Reading intervals list... ")
#intervals <- read.table(INTERVAL_LIST, stringsAsFactors=F)
#split <- strsplit(intervals$V1, "[:-]")
#intervals$chrom <- factor(sapply(split, "[[", 1), levels=c(1:22, "X", "Y"))
#intervals$start <- as.integer(sapply(split, "[[", 2))
#intervals$end <- as.integer(sapply(split, "[[", 3))
#intervals <- intervals[,c("chrom", "start", "end")]
#intervals <- as.list(by(intervals, intervals$chrom, identity))
#cat("done\n")
#names(intervals) <- c(1:22, "X", "Y")

cat("Reading SNVs... ")
d <- read.table("freqs.dat", header=T, sep="\t", fill=T)
d <- droplevels(d[d$chrom %in% c(1:22, "X", "Y"),])
d$chrom <- factor(d$chrom, levels=c(1:22, "X", "Y"))
d$vaf <- d$alt.count/d$depth
cat("done\n")

#cat("Finding exon start positions... ")
#d$start <- mapply(find.start, as.character(d$chrom), d$pos, MoreArgs=list(intervals))
#cat("done\n")

d <- merge(d, s, by.x=c("patient", "sample"), by.y=c("patient", "tumor.sample"))

#d <- merge(d, cnv)
#d <- d[!is.na(d$vaf) & !is.na(d$ratio),]
d <- d[!is.na(d$vaf),]

pdf("test.pdf")
by(d, d$sample, function (ss) {
    mu <- mean(ss$vaf)
    sigma <- sd(ss$vaf)
    ss <- subset(ss, abs((vaf-mu)/sigma) <= 3) # trim outliers
    ss1 <- subset(ss, depth >= 10)
    ss2 <- subset(ss1, vaf >= 0.1)
    if (nrow(ss) == 0 | nrow(ss1) == 0 | nrow(ss2) == 0) return(NULL)
    
    dens <- density(ss$vaf)
    dens1 <- density(ss1$vaf)
    dens2 <- density(ss2$vaf)
    
    ymax <- max(c(dens$y, dens1$y, dens2$y))
    
    plot(dens, col="blue", main=ss[1,"sample"], xlab="variant allelic fraction",
         ylim=c(0, ymax))
    lines(dens2, col="purple")
    abline(v=ss[1,"purity"]/200, lwd=2, lty=2, col="red")
    legend("topright", 
           col=c("blue", "purple", "red"),
           lty=c(1,1,2),
           legend=c("all", "VAF >= 0.1", "pathologist"),
           bg="white")
})
dev.off()
