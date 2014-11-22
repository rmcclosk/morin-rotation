#!/usr/bin/env Rscript

exome.dir <- "/extscratch/morinlab/shared/rmccloskey/colorectal_exome/02_coverage"
genome.dir <- "/extscratch/morinlab/shared/rmccloskey/colorectal_genome/02_coverage"
pooled.dir <- "/extscratch/morinlab/shared/rmccloskey/colorectal_pooled/02_coverage"

dir.create(pooled.dir, showWarnings=F)

options(stringsAsFactors = FALSE)

files <- Sys.glob(file.path(exome.dir, "*.sample_interval_summary"))
files <- files[file.exists(file.path(genome.dir, basename(files)))]

sink("/dev/null")
sapply(files, function (exome.file) {
    genome.file <- file.path(genome.dir, basename(exome.file))
    pooled.file <- file.path(pooled.dir, basename(exome.file))
    
    classes = c("character", rep("numeric", 2), rep("NULL", 6))
    genome.coverage <- read.table(genome.file, header=T, colClasses=classes, sep="\t")
    exome.coverage <- read.table(exome.file, header=T, colClasses=classes, sep="\t")
    
    genome.coverage[,2] <- genome.coverage[,2] + exome.coverage[,2]
    genome.coverage[,3] <- genome.coverage[,3] + exome.coverage[,3]
    write.table(genome.coverage, pooled.file, quote=F, row.names=F, sep="\t")
})
sink()
