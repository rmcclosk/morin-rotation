#!/gsc/software/linux-x86_64-centos5/R-3.0.2/bin/Rscript

library(parallel)
library(HMMcopy)

source(file="/extscratch/morinlab/shared/rmccloskey/settings.conf")

out.dir <- file.path(WORK_DIR, "05_hmmcopy")
wig.dir <- file.path(WORK_DIR, "04_readcounts")
ncpus <- 1

dir.create(out.dir, showWarnings=F)

# Get BAM files for each sample.
sample.data <- read.csv(GENOME_METADATA, header=T)
prefix <- paste(sample.data$patient_id, sample.data$gsc_genome_library_id, "*.bam", sep="_")
sample.data$bam.file <- sapply(prefix, function (x) {
    file <- Sys.glob(file.path(WORK_DIR, "00_bams", x))
    if (length(file) == 0) NA else file
})
sample.data <- sample.data[!is.na(sample.data$bam.file),]

sample.data$wig.file <- file.path(wig.dir, sub(".bam$", ".wig", basename(sample.data$bam.file)))
sample.data <- sample.data[file.exists(sample.data$wig.file),]

# Annotate which samples are normal and which are tumor.
# Delete patients with no normal sample.
# TODO: also annotate with pre- and post-biopsy, from the clinical csv.
sample.data$normal <- grepl("BF", sample.data$sample_id)
normal.counts <- aggregate(normal~patient_id, sample.data, sum)
normal.counts <- normal.counts[normal.counts$normal == 1,]
sample.data <- merge(sample.data, normal.counts, by=c("patient_id"), suffixes=c("", ".count"))

# Delete patients with only one sample.
bam.counts <- aggregate(bam.file~patient_id, sample.data, length)
bam.counts <- bam.counts[bam.counts$bam.file > 1,]
sample.data <- merge(sample.data, bam.counts, by=c("patient_id"), suffixes=c("", ".count"))

sample.data$patient_id <- factor(sample.data$patient_id, levels=unique(sample.data$patient_id))
timestamp <- format(Sys.time(), "%y%m%d%H%M")
sample.data$filename.stem <- file.path(out.dir, sample.data$patient_id, sample.data$sample_id, timestamp)
. <- sapply(dirname(sample.data$filename.stem), dir.create, showWarnings=F, recursive=T)

. <- Sys.setlocale('LC_ALL','C')

write.output <- function (stem, copy, segments=NULL, params=NULL) {
    png(paste0(stem, "_bias.png"))
    par(cex.main = 0.7, cex.lab = 0.7, cex.axis = 0.7, mar = c(4, 4, 2, 0.5))
    plotBias(copy, pch = 20, cex = 0.5)
    dev.off()

    png(paste0(stem, "_correction.png"))
    par(cex.main = 0.7, cex.lab = 0.7, cex.axis = 0.7, mar = c(4, 4, 2, 0.5))
    plotCorrection(copy, pch=".")
    dev.off()

    if (is.null(segments)) return()
    dput(segments, file=paste0(stem, "_segments.dat"))

    if (is.null(params)) params <- HMMsegment(copy, getparam=T)
    dput(params, file=paste0(stem, "_params.dat"))

    rangedDataToSeg(copy, file=paste0(stem, ".seg"))

    png(paste0(stem, "_params.png"))
    plotParam(segments, params)
    dev.off()

    png(paste0(stem, "_segments.png"))
    par(cex.main = 0.5, cex.lab = 0.5, cex.axis = 0.5, mar = c(2, 1.5, 0, 0),
        mgp = c(1, 0.5, 0))
    cols <- stateCols()
    plotSegments(copy, segments, pch = ".", ylab = "Tumour Copy Number", 
                 xlab = "Chromosome Position")
    legend("topleft", c("HOMD", "HETD", "NEUT", "GAIN", "AMPL", "HLAMP"),
           fill = cols, horiz = TRUE, bty = "n", cex = 0.5)
    dev.off()
}

# Apply HMMCopy to each normal/tumor pair.
options(bitmapType="cairo")

sample.data <- sample.data[order(sample.data$patient_id, -sample.data$normal),]
mclapply(sample.data$patient_id, function (p) {
    d <- subset(sample.data, patient_id == p)
    norm.uncorrected.reads <- wigsToRangedData(d[1,"wig.file"], GC_FILE, MAPPABILITY_FILE)
    norm.corrected.copy <- correctReadcount(norm.uncorrected.reads)
    write.output(d[1, "filename.stem"], norm.corrected.copy)

    sapply(2:nrow(d), function (i) {
        tum.uncorrected.reads <- wigsToRangedData(d[i,"wig.file"], GC_FILE, MAPPABILITY_FILE)
        tum.corrected.copy <- correctReadcount(tum.uncorrected.reads)
        tum.corrected.copy$copy <- tum.corrected.copy$copy - norm.corrected.copy$copy
        segments <- HMMsegment(tum.corrected.copy)
        write.output(d[i,"filename.stem"], tum.corrected.copy, segments=segments)
    })
}, mc.cores=ncpus)
