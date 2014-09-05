#!/home/rmccloskey/bin/Rscript

#library(ExomeCNV)

source(file="settings.conf")

coverage.dir <- file.path(WORK_DIR, "01_coverage")

# Keep only rows where at least one comparison can be made.
matching.ncols <- max(count.fields(MATCHING, sep="\t"))
matching <- read.table(MATCHING, sep="\t", fill=T, header=F,
                       na.strings=c("", "Insufficient DNA"), 
                       stringsAsFactors=F, col.names=1:matching.ncols)
matching <- matching[rowSums(!is.na(matching)) > 1,]

# Extract patient ID and file names.
patient.id <- sapply(strsplit(matching[,1], " "), "[[", 1)
matching <- data.frame(apply(matching, 2, function (col) {
    ifelse(is.na(col), NA, sapply(strsplit(col[!is.na(col)], " "), "[[", 2))
}), stringsAsFactors=F)
rownames(matching) <- patient.id

# Make a data frame with information about each bam file.
sample.data <- data.frame(do.call(rbind, lapply(rownames(matching), function (patient.id) {
    row <- matching[patient.id,]
    do.call(rbind, lapply(which(!is.na(row)), function (i) {
        c(patient.id, row[i], i)
    }))
})))
colnames(sample.data) <- c("patient.id", "bam.file", "timepoint")

# Read mean coverage for each sample.
coverage.files <- Sys.glob(file.path(coverage.dir, "*.sample_summary"))
bam.files <- sub(".sample_summary$", ".bam", basename(coverage.files))
mean.coverage <- sapply(coverage.files, function (f) {
    bam.file <- sub(".sample_summary$", ".bam", basename(f))
    coverage <- read.table(f, header=T, fill=T, sep="\t", na.strings=c("", "N/A"))
    coverage[coverage$sample_id == "Total", "mean"]
}, USE.NAMES=F)
coverage <- data.frame(coverage=mean.coverage, bam.file=basename(bam.files))
sample.data <- merge(coverage, sample.data

sample.data
aggregate(mean.coverage~patient.id, sample.data, min)
quit()

# Find factor to down-sample by, to equate means.
scale.factors <- do.call(rbind, apply(matching, 1, function (row) {
    row <- row[!is.na(row)]
    data.frame(scale.factor=min(coverage[row,], na.rm=T)/coverage[row,],
           bam.file=row)
}))

row <- matching[1,]

normal.file <- "01-009_A44993_exome.GRCh37-lite.aln.sample_interval_summary"
diagnosis.file <- "01-009_A44994_exome.GRCh37-lite.aln.sample_interval_summary"
relapse.file <- NA

coverage.files <- c(normal=normal.file, diagnosis=diagnosis.file, relapse=relapse.file)
coverage.paths <- file.path(coverage.dir, coverage.files[!is.na(coverage.files)])
coverage.paths <- ifelse(is.na(coverage.files), NA, coverage.paths)
coverage <- sapply(coverage.paths[!is.na(coverage.paths)], read.coverage.gatk, simplify=F, USE.NAMES=T)


#scale.factor <- coverage[sub(".sample_interval_summary$", ".bam", normal.file),]
#normal$average.coverage <- normal$average.coverage/scale.factor
#normal$coverage <- normal$coverage/scale.factor
