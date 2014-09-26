#!/usr/bin/env Rscript

library(ABSOLUTE)

source(file="/home/rmccloskey/morin-rotation/settings.conf")

coverage.dir <- file.path(WORK_DIR, "02_coverage")
cnv.dir <- file.path(WORK_DIR, "03_cnv")
out.dir <- file.path(WORK_DIR, "12_absolute")
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
d$rdata.path <- file.path(d$out.dir, paste0(d$sample_id, ".ABSOLUTE.RData"))

coverage.file <- file.path(coverage.dir, d$file.name.tumor)
ltt.file <- file.path(cnv.dir, d$patient_id, d$sample_id.tumor, "0.exon.lrr.txt")

# read LTT segments from ExomeCNV
ltt <- read.table(ltt.file[1], 
                  col.names=c("Chromosome", "Start", "End", "Segment_Mean"))
ltt$Chromosome <- sub("chr", "", ltt$Chromosome)

# read coverage
coverage <- read.table(coverage.file[1], header=T,
                       colClasses=c("character", "numeric", rep("NULL", 7)))
colnames(coverage)[2] <- "Num_Probes"
coverage$Chromosome <- sapply(strsplit(coverage$Target, ":"), "[[", 1)
range <- sapply(strsplit(coverage$Target, ":"), "[[", 2)
coverage$Start <- sapply(strsplit(range, "-"), "[[", 1)
coverage$End <- sapply(strsplit(range, "-"), "[[", 2)
coverage <- coverage[,c("Chromosome", "Start", "End", "Num_Probes")]
coverage <- coverage[coverage$Chromosome %in% c(1:22, "X"),]

# we lose a few here, not quite sure why
seg <- merge(ltt, coverage)
seg$Chromosome <- factor(seg$Chromosome, levels=c(1:22, "X", "Y"))
seg <- seg[order(seg$Chromosome, seg$Start),]
write.table(seg, file=paste0(d[1, "out.dir"], ".seg"), row.names=F, quote=F, sep="\t")

cat("Running absolute on sample", d[1, "sample_id.tumor"], "...")
RunAbsolute(paste0(d[1, "out.dir"], ".seg"),
            min.ploidy=0.95, 
            max.ploidy=10, 
            max.sigma.h=0.02, 
            sigma.p=0, 
            platform="Illumina_WES", 
            copy_num_type="total",
            results.dir=d[1, "out.dir"],
            primary.disease="cancer", 
            sample.name=d[1, "sample_id.tumor"], 
            max.as.seg.count=10E10,
            max.non.clonal=0,
            max.neg.genome=0)
cat(" done\n")

CreateReviewObject(d[1, "sample_id.tumor"], d[1, "rdata.path"], d[1, "out.dir"],
                   "total", verbose=T)
calls.path = file.path("test-out", "test.PP-calls_tab.txt")
modes.path = file.path("test-out", "test.PP-modes.data.RData")
output.path = "test-extract"
ExtractReviewedResults(calls.path, "test", modes.path, output.path, "absolute", "total")
