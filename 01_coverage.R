#!/home/rmccloskey/bin/Rscript

source(file="settings.conf")

out.dir <- file.path(WORK_DIR, "01_coverage")
dir.create(out.dir, showWarnings=F)
matching <- read.table(MATCHING, stringsAsFactors=F, fill=T, sep="\t")
files <- c(matching$V1, matching$V2, matching$V3)
files <- files[!files %in% c("", "Insufficient DNA")]
files <- sapply(files, function (f) strsplit(f, " ")[[1]][2])

main <- sapply(files, function (bam.file) {
    bam.path <- file.path(BAM_DIR, bam.file)
    coverage.path <- file.path(out.dir, sub(".bam$", "", bam.file))
    if (! file.exists(paste0(coverage.path, ".sample_summary"))) {
        command <- paste(JAVA_BIN, "-jar", GATK_JAR, "-T DepthOfCoverage",
                         "-omitBaseOutput -omitLocusTable", "-R", HUMAN_REF,
                         "-I", bam.path, "-L", INTERVAL_LIST, "-o", coverage.path)
        cat(command, "\n")
    } else {
        cat("Output file", coverage.path, "already exists\n", file=stderr())
    }
})
