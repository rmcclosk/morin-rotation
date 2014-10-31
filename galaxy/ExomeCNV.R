#!/usr/bin/env Rscript

#library(ExomeCNV)
library(getopt)

spec <- matrix(c(
    "tumor.coverage", "a", 1, "character",
    "normal.coverage", "b", 1, "character",
    "read.length", "c", 1, "integer",
    "admixture.rate", "d", 1, "double",
    "chr.list", "e", 1, "character",
    "sens.exon", "f", 1, "double",
    "spec.exon", "g", 1, "double",
    "opt.exon", "h", 1, "character",
    "sens.segment", "i", 1, "double",
    "spec.segment", "j", 1, "double",
    "opt.segment", "k", 1, "character",
    "coverage.cutoff", "l", 1, "integer",
    "output.stem", "m", 1, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec)

# Print usage.
if (!is.null(opt$help)) {
    cat(getopt(spec, usage=T))
    quit(status=1)
}

# Enforce required arguments.
if (is.null(opt$output.stem)) {
    cat("You must specify a stem for output files\n", file=stderr())
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
    strsplit(opt$chr.list, ",")[[1]]
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

# Create directory for outputs, if needed.
dir.create(dirname(opt$output.stem), showWarnings=F, recursive=T)

# Read in coverage.
tumor.coverage <- read.coverage.gatk(opt$tumor.coverage)
normal.coverage <- read.coverage.gatk(opt$normal.coverage)

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
write.output(eCNV, cnv, opt$output.stem)
