#!/usr/bin/env Rscript

# Copy number analysis on EXOME data.
# TODO: reduce duplication between here and 05_hmmcopy (parsing metadata).

library(ExomeCNV)
library(parallel)

source(file="/home/rmccloskey/morin-rotation/settings.conf")

bam.dir <- file.path(WORK_DIR, "01_fixbams")
coverage.dir <- file.path(WORK_DIR, "02_coverage")
out.dir <- file.path(WORK_DIR, "03_cnv")
chr.list <- paste0("chr", c(1:22, "X", "Y"))
min.sens <- 0.9999
min.spec <- 0.9999
cnv.option <- "spec"
cnv.length <- 100
cnv.contam <- as.numeric(commandArgs(trailingOnly=T)[1])
ncpus <- 1

dir.create(out.dir, showWarnings=F)

sample.data <- read.table(METADATA, header=T)
sample.data$tumor.bam <- file.path(bam.dir, 
                                   paste0(sample.data$tumor.sample, ".bam"))
sample.data$normal.bam <- file.path(bam.dir, 
                                    paste0(sample.data$normal.sample, ".bam"))
sample.data$tumor.coverage.file <- file.path(coverage.dir,
                                             paste0(sample.data$tumor.sample, 
                                                    ".sample_summary"))
sample.data$normal.coverage.file <- file.path(coverage.dir,
                                              paste0(sample.data$normal.sample,
                                                     ".sample_summary"))

sample.data <- sample.data[file.exists(sample.data$normal.coverage.file) &
                           file.exists(sample.data$tumor.coverage.file),]

sample.data$filename.stem <- file.path(out.dir, sample.data$tumor.sample, cnv.contam)
sample.data <- sample.data[!file.exists(paste0(sample.data$filename.stem, ".cnv.txt")),]

if (nrow(sample.data) == 0) {
	cat("Nothing to do\n", stderr())
	quit()
}

# Read mean coverage for each sample.
read.coverage <- function (f) {
    read.table(f, header=T, fill=T, sep="\t")[1,"mean"]
}

sample.data$tumor.coverage <- sapply(sample.data$tumor.coverage.file, 
                                     read.coverage)
sample.data$normal.coverage <- sapply(sample.data$normal.coverage.file, 
                                      read.coverage)

# Scale coverage by the lower of the two samples.
worst.coverage <- min(sample.data$tumor.coverage, sample.data$normal.coverage)
sample.data$tumor.scale <- worst.coverage/sample.data$tumor.coverage
sample.data$normal.scale <- worst.coverage/sample.data$normal.coverage

. <- sapply(dirname(sample.data$filename.stem), dir.create, showWarnings=F, recursive=T)

# Read in GATK coverage and scale down.
sample.data$tumor.interval.file <- gsub("sample_summary$", 
                                        "sample_interval_summary", 
                                        sample.data$tumor.coverage.file)
sample.data$normal.interval.file <- gsub("sample_summary$", 
                                         "sample_interval_summary", 
                                         sample.data$normal.coverage.file)

gatk.scale <- function (f, scale.factor) {
    d <- read.coverage.gatk(f)
    d$coverage <- d$coverage * scale.factor
    d$average.coverage <- d$average.coverage * scale.factor
    d$base.with..10.coverage <- d$base.with..10.coverage * scale.factor
    d
}

# Run ExomeCNV.
options(bitmapType="cairo")
mclapply(1:nrow(sample.data), function (i) {
    norm.coverage <- gatk.scale(sample.data[i, "normal.interval.file"],
                                sample.data[i, "normal.scale"])
    tum.coverage <- gatk.scale(sample.data[i, "tumor.interval.file"],
                               sample.data[i, "tumor.scale"])
    logR <- calculate.logR(norm.coverage, tum.coverage)
    eCNV <- do.call(rbind, lapply(chr.list, function (chr) {
        idx <- norm.coverage$chr == chr
        classify.eCNV(normal=norm.coverage[idx,], tumor=tum.coverage[idx,], 
                      logR=logR[idx], min.spec=min.spec, 
                      min.sens=min.sens, option=cnv.option, 
                      c=cnv.contam, l=cnv.length)
    }))
    cnv <-  multi.CNV.analyze(norm.coverage, tum.coverage, logR=logR, 
                              all.cnv.ls=list(eCNV), coverage.cutoff=5, 
                              min.spec=0.99, min.sens=0.99, option="auc", 
                              c=cnv.contam)
    write.output(eCNV, cnv, sample.data[i, "filename.stem"])
}, mc.cores = ncpus)
