#!/usr/bin/env Rscript

# This script summarizes the results of running ABSOLUTE on multiple samples.
# 
# Arguments
# ---------
#    -p: path containing RData files produced by ABSOLUTE.R
#    -n: name of summary files to produce (default = "summary")
#    -d: directory to put summary files in (default = ".")

library(ABSOLUTE)
library(getopt)

spec = matrix(c(
    "help" , "h", 0, "logical",
    "path", "p", 1, "character",
    "outfile-name", "n", 1, "character",
    "out-dir", "d", 1, "character"
), byrow=TRUE, ncol=4)
opt <- getopt(spec)

usage <- function (error.msg="") {
    cat(error.msg, "\n")
    cat(getopt(spec, usage=T))
    quit(status=1)
}

if (!is.null(opt$help)) usage()
if (is.null(opt$path)) usage ("Must specify a path containing RData files")
if (is.null(opt$`outfile-name`)) opt$`outfile-name` <- "summary"
if (is.null(opt$`out-dir`)) opt$`out-dir` <- "."

absolute.files <- Sys.glob(file.path(opt$path, "*.ABSOLUTE.RData"))
log.file <- paste0(file.path(opt$`out-dir`, opt$`outfile-name`), ".log")
sink(log.file)
CreateReviewObject(opt$`outfile-name`, absolute.files, opt$`out-dir`, "total", verbose=TRUE)
sink()
