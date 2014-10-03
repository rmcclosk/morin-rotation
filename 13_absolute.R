#!/usr/bin/env Rscript

library(ABSOLUTE)

source(file="/home/rmccloskey/morin-rotation/settings.conf")
options("scipen"=100)

sample.data <- read.table("sample_info.dat", header=T)

coverage.dir <- file.path(WORK_DIR, "02_coverage")
#cnv.dir <- file.path(WORK_DIR, "05_hmmcopy", "10000")
cnv.dir <- file.path(WORK_DIR, "03_exomecnv")
out.dir <- file.path(WORK_DIR, "13_absolute")
dir.create(out.dir, showWarnings=F)

# read metadata about samples
d <- read.csv(EXOME_METADATA, header=T, stringsAsFactors=F)
d <- d[,c("sample_id", "patient_id", "gsc_exome_library_id")]
d$file.name <- paste(d$patient_id, d$gsc_exome_library_id, 
                     "exome.GRCh37-lite.aln.sample_interval_summary", sep="_")
normal <- d[grepl("(BF|WB)", d$sample_id),]
tumor <- d[!grepl("(BF|WB)", d$sample_id),]
d <- merge(normal, tumor, by=c("patient_id"), suffixes=c(".normal", ".tumor"))
d$out.dir <- file.path(out.dir, d$sample_id.tumor)
d$rdata.path <- file.path(out.dir, paste0(d$sample_id.tumor, ".ABSOLUTE.RData"))

d$ltt.file <- file.path(cnv.dir, d$patient_id, paste0(d$sample_id.tumor, "_default_segments.dat"))
d <- d[file.exists(d$ltt.file),]

sapply(1:nrow(d), function (i) {
    # read segments from HMMcopy
    seg <- read.table(d[i, "ltt.file"], header=T,
                      colClasses=c("character", rep("numeric", 2), "NULL", "numeric"),
                      col.names=c("Chromosome", "Start", "End", "state", "Segment_Mean"))
    seg$Chromosome <- sub("chr", "", seg$Chromosome)
    seg$Num_Probes <- round((seg$End-seg$Start)/10000)
    
    seg$Chromosome <- factor(seg$Chromosome, levels=c(1:22, "X", "Y"))
    seg <- seg[order(seg$Chromosome, seg$Start),]
    write.table(seg, file=paste0(d[i, "out.dir"], ".seg"), row.names=F, quote=F, sep="\t")
    
    cat("Running absolute on sample", d[i, "sample_id.tumor"], "...")
    sink(paste0(d[i, "out.dir"], ".log"))
    RunAbsolute(paste0(d[i, "out.dir"], ".seg"),
                min.ploidy=0.95, 
                max.ploidy=10, 
                max.sigma.h=0.02, 
                sigma.p=0, 
                platform="Illumina_WES", 
                copy_num_type="total",
                results.dir=out.dir,
                primary.disease="cancer", 
                sample.name=d[i, "sample_id.tumor"], 
                max.as.seg.count=10E10,
                max.non.clonal=0,
                max.neg.genome=0,
                verbose=T)
    sink()
    cat(" done\n")
    
    CreateReviewObject(d[i, "sample_id.tumor"], d[i, "rdata.path"], 
                       out.dir, "total", verbose=T)
    
    calls.path = file.path(out.dir,
                           paste0(d[i, "sample_id.tumor"], ".PP-calls_tab.txt"))
    modes.path = file.path(out.dir,
                           paste0(d[i, "sample_id.tumor"], ".PP-modes.data.RData"))
    ExtractReviewedResults(calls.path, d[i, "sample_id.tumor"], modes.path, out.dir, "absolute", "total")
})
