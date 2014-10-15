#!/usr/bin/env Rscript

source(file="settings.conf")

abs.dir <- file.path(WORK_DIR, "13_absolute")
d <- read.table(METADATA, header=T)

read.absolute <- function (sample) {
    abs.file <- file.path(abs.dir, paste0(sample, ".PP-calls_tab.txt"))
    if (file.exists(abs.file)) {
        abs.data <- read.table(abs.file, header=T, sep="\t", fill=T)
        abs.data$Cancer.DNA.fraction
    } else {
        NA
    }
}

cor.test(sapply(d$tumor.sample, read.absolute), d$purity, na.rm=T)
