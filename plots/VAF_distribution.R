#!/usr/bin/env Rscript

# Plot a kernel density of variant allele fractions from a MAF file.
# Optionally, also add vertical lines at anticipated VAF peaks for various copy
# numbers.

library(getopt)

spec <- matrix(c(
    "help", "h", 0, "logical",
    "maf", "m", 1, "character",
    "purity", "p", 1, "double",
    "outfile", "o", 1, "character",
    "title", "t", 1, "character"
), byrow=T, ncol=4)
opt <- getopt(spec)

usage <- function (message=NA) {
    if (!is.na(message))
        cat(message, "\n")
    cat(getopt(spec, usage=TRUE))
    quit(status=1)
}

if (!is.null(opt$help))
    usage()
if (is.null(opt$maf))
    usage("You must specify a MAF file with -m")
if (!file.exists(opt$maf))
    usage("Provided MAF file does not exist")
if (is.null(opt$outfile))
    usage("You must specify an output PDF file with -o")

if (!grepl(".pdf$", opt$outfile))
    opt$outfile <- paste0(opt$outfile, ".pdf")

maf <- read.delim(opt$maf)
maf$vaf <- maf$t_alt_count/(maf$t_alt_count+maf$t_ref_count)

pdf(opt$outfile)
plot(density(maf$vaf), main="", xlab="variant allele fraction")
if (!is.null(opt$title))
    title(main=opt$title)
if (!is.null(opt$purity)) {
    t <- opt$purity
    abline(v=t, lty=2, col=2)
    abline(v=t/2, lty=2, col=4)
    abline(v=t/(3*t+2*(1-t)), lty=2, col=3)
    abline(v=2*t/(3*t+2*(1-t)), lty=2, col=3)
    legend("topright", col=c(2,4,3), lty=2, bg="white",
           legend=c("copy number 1", "copy number 2", "copy number 3"),
           title="Clonal VAF")
}
dev.off()
