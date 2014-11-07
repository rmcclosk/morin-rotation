#!/usr/bin/env Rscript

library(ABSOLUTE)
library(parallel)

source(file="/home/rmccloskey/morin-rotation/settings.conf")
options("scipen"=100, "warn"=1)
ncpus <- 1

d <- read.table(METADATA, header=T, fill=T)
d <- d[!is.na(d$normal.sample),]

coverage.dir <- file.path(WORK_DIR, "02_coverage")
cnv.dir <- file.path(WORK_DIR, "03_cnv")
out.dir <- file.path(WORK_DIR, "13_absolute")
dir.create(out.dir, showWarnings=F)

d$stem <- file.path(out.dir, d$tumor.sample)
d$rdata.path <- paste0(d$stem, ".ABSOLUTE.RData")

d$ltt.file <- file.path(cnv.dir, d$tumor.sample, "0.segment.lrr.txt")
d <- d[file.exists(d$ltt.file),]

stopifnot(nrow(d) > 0)

all.seg <- do.call(rbind, lapply(1:nrow(d), function (i) {
    # read segments 
    seg <- read.table(d[i, "ltt.file"], header=T,
                      colClasses=c("character", rep("numeric", 2), "numeric"),
                      col.names=c("Chromosome", "Start", "End", "Segment_Mean"))
    seg$Chromosome <- sub("chr", "", seg$Chromosome)
    seg$Num_Probes <- floor((seg$End-seg$Start)/1000)
    
    seg$Chromosome <- factor(seg$Chromosome, levels=c(1:22, "X", "Y"))
    seg$Sample <- d[i, "tumor.sample"]
    seg[order(seg$Chromosome, seg$Start),]
}))
write.table(all.seg, file="absolute_input.seg", quote=F, row.names=F, col.names=T)

cat("Running absolute...")
sink(out.dir, "absolute.log")
RunAbsolute(paste0(d[i, "stem"], ".seg"),
            min.ploidy=0.95, 
            max.ploidy=10, 
            max.sigma.h=0.02, 
            sigma.p=0, 
            platform="Illumina_WES", 
            copy_num_type="total",
            results.dir=out.dir,
            primary.disease="cancer", 
            sample.name=d[i, "tumor.sample"], 
            max.as.seg.count=10E10,
            max.non.clonal=0,
            max.neg.genome=0,
            verbose=T)
sink()
cat(" done\n")

CreateReviewObject(d[i, "tumor.sample"], d[i, "rdata.path"], 
                   out.dir, "total", verbose=T)

calls.path = file.path(out.dir,
                       paste0(d[i, "tumor.sample"], ".PP-calls_tab.txt"))
modes.path = file.path(out.dir,
                       paste0(d[i, "tumor.sample"], ".PP-modes.data.RData"))
ExtractReviewedResults(calls.path, d[i, "tumor.sample"], modes.path, out.dir, "absolute", "total")
file.rename(file.path(out.dir, "reviewed"), d[i, "stem"])
