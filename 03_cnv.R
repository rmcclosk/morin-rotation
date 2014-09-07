#!/home/rmccloskey/bin/Rscript

library(ExomeCNV)
library(parallel)

source(file="settings.conf")

coverage.dir <- file.path(EXOME_WORK_DIR, "01_coverage")
out.dir <- file.path(EXOME_WORK_DIR, "03_cnv")
chr.list <- paste0("chr", c(1:22, "X", "Y"))
chr.list <- c("chr1")
min.sens <- 0.9999
min.spec <- 0.9999
cnv.option <- "spec"
cnv.length <- 100
cnv.contam <- 0.5
ncpus <- 1

dir.create(out.dir, showWarnings=F)

# Keep only rows where at least one comparison can be made.
matching.ncols <- max(count.fields(EXOME_MATCHING, sep="\t"))
matching <- read.table(EXOME_MATCHING, sep="\t", fill=T, header=F,
                       na.strings=c("", "Insufficient DNA"), 
                       stringsAsFactors=F, col.names=1:matching.ncols)
matching <- matching[rowSums(!is.na(matching)) > 1,]

# Extract patient ID and file names.
patient.id <- sapply(strsplit(matching[,1], " "), "[[", 1)
matching <- data.frame(apply(matching, 2, function (col) {
    ifelse(is.na(col), NA, sapply(strsplit(col[!is.na(col)], " "), "[[", 2))
}), stringsAsFactors=F)

# Make a data frame with information about each bam file.
sample.data <- data.frame(
    patient.id=rep(patient.id, rowSums(!is.na(matching))),
    bam.file=t(matching)[!is.na(t(matching))],
    timepoint=((which(t(!is.na(matching)))-1) %% ncol(matching)) + 1
)

# Read mean coverage for each sample.
coverage.files <- Sys.glob(file.path(coverage.dir, "*.sample_summary"))
tmp <- data.frame(
    coverage.file=sub(".sample_summary$", ".sample_interval_summary", basename(coverage.files)),
    bam.file=sub(".sample_summary$", ".bam", basename(coverage.files)),
    coverage=sapply(coverage.files, function (f) {
        read.table(f, header=T, fill=T, sep="\t")[1,"mean"]
    })
)
sample.data <- merge(sample.data, tmp)

# Remove patients with only one BAM file.
# TODO: delete this when BAM file corruption issues are resolved.
counts <- aggregate(bam.file~patient.id, sample.data, length)
keep.patients <- subset(counts, bam.file > 1, select=c("patient.id"))
sample.data <- merge(sample.data, keep.patients)
sample.data$patient.id <- factor(sample.data$patient.id, 
                                 levels=keep.patients$patient.id)

# Find the minimum coverage by patient and scaling factors.
worst.coverage <- aggregate(coverage~patient.id, sample.data, min)
sample.data <- merge(sample.data, worst.coverage, by=c("patient.id"), 
                     suffixes=c("", ".worst"))
sample.data$scale.factor <- sample.data$coverage.worst/sample.data$coverage

sample.data <- sample.data[1:2,]
sample.data$patient.id <- factor(sample.data$patient.id, levels=unique(sample.data$patient.id))

coverage <- mclapply(1:nrow(sample.data), function (i) {
    d <- read.coverage.gatk(file.path(coverage.dir, sample.data[i, "coverage.file"]))
    scale.factor <- sample.data[i, "scale.factor"]
    d$coverage <- d$coverage*scale.factor
    d$average.coverage <- d$average.coverage*scale.factor
    d$base.with..10.coverage <- d$base.with..10.coverage*scale.factor
    d
}, mc.cores=ncpus)
names(coverage) <- sample.data$bam.file

sample.data <- sample.data[order(sample.data$patient.id, sample.data$timepoint),]
by(sample.data, sample.data$patient.id, function (x) {
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
        patient.id <- x[pair[1], "patient.id"]
        timepoint <- x[pair[2], "timepoint"]
        filename.stem <- paste0("patient", patient.id, "_timepoint", timepoint)
        write.output(eCNV, cnv, file.path(out.dir, filename.stem))
    })
})
