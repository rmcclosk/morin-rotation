#!/home/rmccloskey/bin/Rscript

library(parallel)
library(HMMcopy)

source(file="settings.conf")

out.dir <- file.path(WORK_DIR, "05_hmmcopy")
wig.dir <- file.path(WORK_DIR, "04_readcounts")
ncpus <- 1

dir.create(out.dir, showWarnings=F)

# Get BAM files for each sample.
sample.data <- read.csv(GENOME_METADATA, header=T)
prefix <- paste(sample.data$patient_id, sample.data$gsc_genome_library_id, "*.bam", sep="_")
sample.data$bam.file <- sapply(prefix, function (x) {
    file <- Sys.glob(file.path(GENOME_BAM_DIR, x))
    if (length(file) == 0) NA else file
})
sample.data <- sample.data[!is.na(sample.data$bam.file),]

sample.data$wig.file <- file.path(wig.dir, sub(".bam$", ".wig", basename(sample.data$bam.file)))
sample.data <- sample.data[file.exists(sample.data$wig.file),]

# Annotate which samples are normal and which are tumor.
# Delete patients with no normal sample.
# TODO: also annotate with pre- and post-biopsy, from the clinical csv.
sample.data$normal <- grepl("BF", sample.data$sample_id)
sample.data <- sample.data[order(sample.data$patient_id, -sample.data$normal),]
normal.counts <- aggregate(normal~patient_id, sample.data, sum)
normal.counts <- normal.counts[normal.counts$normal == 1,]
sample.data <- merge(sample.data, normal.counts, by=c("patient_id"), suffixes=c("", ".count"))

# Delete patients with only one sample.
bam.counts <- aggregate(bam.file~patient_id, sample.data, length)
bam.counts <- bam.counts[bam.counts$bam.file > 1,]
sample.data <- merge(sample.data, bam.counts, by=c("patient_id"), suffixes=c("", ".count"))

sample.data$patient_id <- factor(sample.data$patient_id, levels=unique(sample.data$patient_id))

Sys.setlocale('LC_ALL','C')
wigsToRangedData(sample.data$wig.file[1], GC_FILE, MAPPABILITY_FILE)
quit()
# Apply HMMCopy to each normal/tumor pair.
by(sample.data, sample.data$patient_id, function (d) {
    normal <- d[1,]
    sapply(2:nrow(d), function (i) {
        tum.uncorrected.reads <- wigsToRangedData(d[i,"wig.file"], GC_FILE, MAPPABILITY_FILE)
    })
})
