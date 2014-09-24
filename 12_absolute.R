#!/usr/bin/env Rscript

#library(ABSOLUTE)

source(file="/home/rmccloskey/morin-rotation/settings.conf")

coverage.dir <- file.path(WORK_DIR, "02_coverage")
cnv.dir <- file.path(WORK_DIR, "02_coverage")

d <- read.csv(EXOME_METADATA, header=T)
d$file.name <- paste(d$patient_id, d$gsc_exome_library_id, 
                     "exome.GRCh37-lite.aln.sample_interval_summary", sep="_")
d <- d[,c("sample_id", "patient_id", "gsc_exome_library_id")]
normal <- d[grepl("(BF|WB)", d$sample_id),]
tumor <- d[!grepl("(BF|WB)", d$sample_id),]
d <- merge(normal, tumor, by=c("patient_id"), suffixes=c(".normal", ".tumor"))

coverage.file <- file.path(coverage.dir, d$file.name)
file.exists(coverage.file)
