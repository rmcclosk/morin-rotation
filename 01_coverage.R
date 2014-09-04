#!/home/rmccloskey/bin/Rscript

source(file="settings.R")

dir.create(coverage.dir, showWarnings=F)
matching <- read.table(matching.file, stringsAsFactors=F, fill=T, sep="\t")
files <- c(matching$V1, matching$V2, matching$V3)
files <- files[!files %in% c("", "Insufficient DNA")]
files <- sapply(files, function (f) strsplit(f, " ")[[1]][2])

main <- sapply(files, function (bam.file) {
    bam.path <- file.path(bam.dir, bam.file)
    coverage.path <- file.path(coverage.dir, sub(".bam$", "", bam.file))
    if (! file.exists(paste0(coverage.path, ".sample_summary"))) {
        command <- sprintf("%s -jar %s -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R %s -I %s -L %s -o %s",
                           java.bin, gatk.jar, human.fasta, bam.path, exome.region.list, coverage.path)
        system(command)
    } else {
        cat("Output file", coverage.path, "already exists")
    }
})
