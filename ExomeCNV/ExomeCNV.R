#!/usr/bin/env Rscript

library(ExomeCNV)
library(getopt)

spec <- matrix(c(
    "tumor.coverage", "a", 1, "character",
    "normal.coverage", "b", 1, "character",
    "read.length", "c", 1, "integer",
    "admixture.rate", "d", 1, "double",
    "chr.list", "e", 1, "character",
    "sens.exon", "f", 1, "double",
    "spec.exon", "r", 1, "double",
    "opt.exon", "h", 1, "character",
    "sens.segment", "i", 1, "double",
    "spec.segment", "j", 1, "double",
    "opt.segment", "k", 1, "character",
    "coverage.cutoff", "l", 1, "integer",
    "cnv.png", "m", 1, "character",
    "cnv.txt", "n", 1, "character",
    "exon.lrr.txt", "o", 1, "character",
    "segment.copynumber.txt", "p", 1, "character",
    "segment.lrr.txt", "q", 1, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec)

# Print usage.
if (!is.null(opt$help)) {
    cat(getopt(spec, usage=T))
    quit(status=1)
}

# Enforce required arguments.
if (is.null(opt$cnv.png) & is.null(opt$cnv.txt) & is.null(opt$exon.lrr.txt) &
        is.null(opt$segment.copynumber.txt) & is.null(opt$segment.lrr.txt)) {
    cat("You must specify names for at least one output file\n", file=stderr())
    quit(status=1)
}
if (is.null(opt$tumor.coverage) | is.null(opt$normal.coverage)) {
    cat("You must specify coverage files for the tumor and normal samples\n", file=stderr())
    quit(status=1)
}
if (is.null(opt$read.length)) {
    cat("You must specify the average read length\n", file=stderr())
    quit(status=1)
}

# Set default arguments.
opt$chr.list <- if (is.null(opt$chr.list)) {
    paste0("chr", c(1:22, "X", "Y"))
} else {
    paste0("chr", strsplit(opt$chr.list, ",")[[1]])
}
if (is.null(opt$sens.exon)) { opt$sens.exon <- 0.9999 }
if (is.null(opt$spec.exon)) { opt$spec.exon <- 0.9999 }
if (is.null(opt$opt.exon)) { opt$opt.exon <- "spec" }
if (is.null(opt$sens.segment)) { opt$sens.segment <- 0.99 }
if (is.null(opt$spec.segment)) { opt$spec.segment <- 0.99 }
if (is.null(opt$opt.segment)) { opt$opt.segment <- "auc" }
if (is.null(opt$admixture.rate)) { opt$admixture.rate <- 0.3 }
if (is.null(opt$coverage.cutoff)) { opt$coverage.cutoff <- 5 }

# Warn when in violation of the user guide.
if (opt$spec.exon < 0.9999 | opt$opt.exon != "spec") {
    warning(paste("The ExomeCNV user guide recommends a high specificity (0.9999) and",
                  "option \"spec\" to be conservative against false positives when",
                  "calling CNV for each exon."))
}
if ((opt$opt.exon == "auc" & opt$sens.exon != opt$spec.exon) | 
    (opt$opt.segment == "auc" & opt$sens.segment != opt$spec.segment)) {
    warning(paste("The ExomeCNV user guide recommends setting sensitivity equal to",
                  "specificity when the option \"auc\" is used."))
}

# Read in coverage.
read.coverage <- function (file.name) {
    if (grepl(".bz2$", file.name))
        read.coverage.gatk(bzfile(opt$tumor.coverage))
    else if (grepl(".gz$", file.name))
        read.coverage.gatk(gzfile(opt$tumor.coverage))
    else
        read.coverage.gatk(opt$tumor.coverage)
}
tumor.coverage <- read.coverage(opt$tumor.coverage)
normal.coverage <- read.coverage(opt$normal.coverage)

# Run ExomeCNV.
options(bitmapType="cairo")
logR <- calculate.logR(normal.coverage, tumor.coverage)
eCNV <- do.call(rbind, lapply(opt$chr.list, function (chr) {
    idx <- normal.coverage$chr == chr
    classify.eCNV(normal=normal.coverage[idx,], tumor=tumor.coverage[idx,],
                  logR=logR[idx], min.spec=opt$spec.exon, min.sens=opt$sens.exon,
                  option=opt$opt.exon, c=opt$admixture.rate, l=opt$read.length)
}))
cnv <-  multi.CNV.analyze(normal.coverage, tumor.coverage, logR=logR,
                          all.cnv.ls=list(eCNV),
                          coverage.cutoff=opt$coverage.cutoff,
                          min.spec=opt$spec.segment, min.sens=opt$sens.segment,
                          option=opt$opt.segment, c=opt$admixture.rate)

# Copy outputs to their proper locations.
stem <- tempfile()
write.output(eCNV, cnv, stem)
copy.compress <- function (src, dest) {
    if (!is.null(dest)) {
        if (grepl(".bz2$", dest) | grepl(".gz$", dest)) {
            ext <- regmatches(dest, regexpr("[.][^.]+$", dest))
            base <- sub(ext, "", dest)
            compress <- ifelse(ext == ".bz2", "bzip2", "gzip")
            file.rename(src, base)
            system2(compress, args=c(base))
        } else {
            file.rename(src, dest)
        }
    }
}
copy.compress(paste0(stem, ".cnv.png"), opt$cnv.png)
copy.compress(paste0(stem, ".cnv.txt"), opt$cnv.txt)
copy.compress(paste0(stem, ".exon.lrr.txt"), opt$exon.lrr.txt)
copy.compress(paste0(stem, ".segment.copynumber.txt"), opt$segment.copynumber.txt)
copy.compress(paste0(stem, ".segment.lrr.txt"), opt$segment.lrr.txt)
