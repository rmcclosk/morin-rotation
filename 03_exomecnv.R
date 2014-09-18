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
ncpus <- 4

dir.create(out.dir, showWarnings=F)

# Get BAM files for each sample.
if (grepl("colorectal", WORK_DIR)) {
    sample.data <- read.csv(EXOME_METADATA, header=T)
    prefix <- paste(sample.data$patient_id, sample.data$gsc_exome_library_id, "*.bam", sep="_")
    sample.data$bam.file <- sapply(prefix, function (x) {
        file <- Sys.glob(file.path(bam.dir, x))
        if (length(file) == 0) NA else file
    })
    sample.data <- sample.data[!is.na(sample.data$bam.file),]
    sample.data$normal <- grepl("(BF|WB)", sample.data$sample_id)
} else if (grepl("dlbcl", WORK_DIR)) {
    sample.data <- data.frame(bam.file=Sys.glob(file.path(bam.dir, "*.bam")))
    sample.data$bam.file <- as.character(sample.data$bam.file)

    match <- regexpr("^(PT)?[0-9]+", basename(sample.data$bam.file))
    sample.data$patient_id <- regmatches(basename(sample.data$bam.file), match)

    match <- regexpr("(T2|N|GL|T)", sub("PT", "", basename(sample.data$bam.file)))
    sample.data$sample <- regmatches(sub("PT", "", basename(sample.data$bam.file)), match)

    sample.data$normal <- sample.data$sample %in% c("GL", "N")
}

sample.data$coverage.file <- file.path(coverage.dir, sub(".bam$", ".sample_summary", basename(sample.data$bam.file)))
sample.data <- sample.data[file.exists(sample.data$coverage.file),]

# Delete patients with no normal sample.
# TODO: also annotate with pre- and post-biopsy, from the clinical csv.
sample.data <- sample.data[order(sample.data$patient_id, -sample.data$normal),]
normal.counts <- aggregate(normal~patient_id, sample.data, sum)
normal.counts <- normal.counts[normal.counts$normal == 1,]
sample.data <- merge(sample.data, normal.counts, by=c("patient_id"), suffixes=c("", ".count"))

sample.data$patient_id <- factor(sample.data$patient_id, levels=unique(sample.data$patient_id))
sample.data$filename.stem <- file.path(out.dir, sample.data$patient_id, sample.data$sample_id, cnv.contam)

# Remove late time points which have already been analysed.
sample.data <- sample.data[!file.exists(paste0(sample.data$filename.stem, ".cnv.txt")),]

if (all(sample.data$normal)) {
	cat("Nothing to do\n", stderr())
	quit()
}

# Delete patients with only one sample.
bam.counts <- aggregate(bam.file~patient_id, sample.data, length)
bam.counts <- bam.counts[bam.counts$bam.file > 1,]
sample.data <- merge(sample.data, bam.counts, by=c("patient_id"), suffixes=c("", ".count"))

# Read mean coverage for each sample.
sample.data$coverage <- sapply(sample.data$coverage.file, function (f) {
    read.table(f, header=T, fill=T, sep="\t")[1,"mean"]
})

# Remove patients with only one BAM file.
counts <- aggregate(bam.file~patient_id, sample.data, length)
keep.patients <- subset(counts, bam.file > 1, select=c("patient_id"))
sample.data <- merge(sample.data, keep.patients)
sample.data$patient_id <- factor(sample.data$patient_id, 
                                 levels=keep.patients$patient_id)

# Find the minimum coverage by patient and scaling factors.
worst.coverage <- aggregate(coverage~patient_id, sample.data, min)
sample.data <- merge(sample.data, worst.coverage, by=c("patient_id"), 
                     suffixes=c("", ".worst"))
sample.data$scale.factor <- sample.data$coverage.worst/sample.data$coverage

sample.data$patient_id <- factor(sample.data$patient_id, levels=unique(sample.data$patient_id))
. <- sapply(dirname(sample.data$filename.stem), dir.create, showWarnings=F, recursive=T)

# Read in GATK coverage and scale down.
sample.data <- sample.data[order(sample.data$patient_id, -sample.data$normal),]

sample.data$coverage.interval.file <- gsub("sample_summary$", "sample_interval_summary", sample.data$coverage.file)
coverage <- mclapply(1:nrow(sample.data), function (i) {
    d <- read.coverage.gatk(sample.data[i, "coverage.interval.file"])
    scale.factor <- sample.data[i, "scale.factor"]
    d$coverage <- d$coverage*scale.factor
    d$average.coverage <- d$average.coverage*scale.factor
    d$base.with..10.coverage <- d$base.with..10.coverage*scale.factor
    d
}, mc.cores=ncpus)
names(coverage) <- sample.data$bam.file

# Run ExomeCNV.
options(bitmapType="cairo")
. <- by(sample.data, sample.data$patient_id, function (x) {
    pairs <- expand.grid(1, 2:nrow(x))
    apply(pairs, 1, function (pair) {
        bam.files <- x[pair, "bam.file"]
        coverage1 <- coverage[[bam.files[1]]]
        coverage2 <- coverage[[bam.files[2]]]
        logR <- calculate.logR(coverage1, coverage2)
        eCNV <- do.call(rbind, mclapply(chr.list, function (chr) {
            idx <- coverage1$chr == chr
            classify.eCNV(normal=coverage1[idx,], tumor=coverage2[idx,], 
                          logR=logR[idx], min.spec=min.spec, 
                          min.sens=min.sens, option=cnv.option, 
                          c=cnv.contam, l=cnv.length)
        }, mc.cores=ncpus))
        cnv <-  multi.CNV.analyze(coverage1, coverage2, logR=logR, 
                                  all.cnv.ls=list(eCNV), coverage.cutoff=5, 
                                  min.spec=0.99, min.sens=0.99, option="auc", 
                                  c=cnv.contam)
        write.output(eCNV, cnv, x[pair[2], "filename.stem"])
    })
})
