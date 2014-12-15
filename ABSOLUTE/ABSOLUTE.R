#!/usr/bin/env Rscript

# Run ABSOLUTE on segments, and optionally a MAF file, for an individual
# sample. This is only the first step of ABSOLUTE. The second step is to
# summarize the results from all the samples collectively (see summarize.R).
# 
# Arguments
# ---------
# See the [official documentation](http://www.broadinstitute.org/cancer/software/genepattern/modules/docs/ABSOLUTE/1?print=true) 
# for the meaning of the parameters not documented here.
#
#   segments: file containing CNV segments, optionally compressed with gzip or
#             bzip2. Has 5 columns: Sample, Chromosome, Start, End, Num_Probes,
#             Segment_Mean.
#
#   maf-file: optional. A MAF file containing SNVs, with the columns
#             t_ref_count and t_alt_count filled in.
#
#     outdir: directory to write output files.

library(ABSOLUTE)
library(getopt)

options("scipen"=100)

spec = matrix(c(
    "help" , "h", 0, "logical",
    "segments", "s", 1, "character",
    "platform", "n", 1, "character",
    "outdir", "o", 1, "character",
    "primary-disease", "d", 1, "character",
    "sigma-p", "p", 1, "numeric",
    "max-sigma-h", "j", 1, "numeric",
    "min-ploidy", "i", 1, "numeric",
    "max-ploidy", "a", 1, "numeric",
    "max-as-seg-count", "c", 1, "numeric",
    "max-non-clonal", "l", 1, "numeric",
    "max-neg-genome", "x", 1, "numeric",
    "maf-file", "m", 1, "character",
    "min-mut-af", "f", 1, "numeric"
), byrow=TRUE, ncol=4)
opt <- getopt(spec)

usage <- function (error.msg="") {
    cat(error.msg, "\n")
    cat(getopt(spec, usage=T))
    quit(status=1)
}

if (!is.null(opt$help)) usage()
if (is.null(opt$segments)) usage("Must specify a segments file")

if (is.null(opt$platform)) opt$platform <- "Illumina_WES"
if (is.null(opt$outdir)) opt$outdir <- "."
if (is.null(opt$`primary-disease`)) opt$`primary-disease` <- "cancer"
if (is.null(opt$`sigma-p`)) opt$`sigma-p` <- 0.01
if (is.null(opt$`max-sigma-h`)) opt$`max-sigma-h` <- 0.02
if (is.null(opt$`min-ploidy`)) opt$`min-ploidy` <- 1
if (is.null(opt$`max-ploidy`)) opt$`max-ploidy` <- 6
if (is.null(opt$`max-as-seg-count`)) opt$`max-as-seg-count` <- 1500
if (is.null(opt$`max-non-clonal`)) opt$`max-non-clonal` <- 0.1
if (is.null(opt$`max-neg-genome`)) opt$`max-neg-genome` <- 0.01
if (is.null(opt$`min-mut-af`)) opt$`min-mut-af` <- 0.05

archive.or.file <- function (fn) {
    if (grepl(".bz2$", fn))
        bzfile(fn)
    if (grepl(".gz$", fn))
        gzfile(fn)
    else
        fn
}

seg.file <- archive.or.file(opt$segments)
if (!is.null(opt$`maf-file`)) {
    maf.file <- archive.or.file(opt$`maf-file`)
} else { 
    maf.file <- NULL
}

sample <- sub(".seg(.gz|.bz2)?", "", basename(opt$segments))
dir.create(opt$outdir, showWarnings=F)
stem <- file.path(opt$outdir, sample)

sink(file=paste0(stem, ".log"))
RunAbsolute(seg.file,
            opt$`sigma-p`, 
            opt$`max-sigma-h`,
            opt$`min-ploidy`, 
            opt$`max-ploidy`,
            opt$`primary-disease`, 
            opt$platform, 
            sample, 
            opt$outdir, 
            opt$`max-as-seg-count`,
            opt$`max-non-clonal`, 
            opt$`max-neg-genome`, 
            "total", 
            verbose=TRUE,
            maf.fn=maf.file,
            min.mut.af=opt$`min-mut-af`)
sink()
