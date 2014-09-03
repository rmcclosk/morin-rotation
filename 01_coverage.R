#!/usr/bin/env Rscript

source(file="settings.R")

dir.create(coverage.dir, showWarnings=F)
matching <- read.table(matching.file, stringsAsFactors=F)
main <- sapply(matching[,c(2,4,6)], function (bam.file) {
    bam.path <- file.path(bam.dir, bam.file)
    coverage.path <- file.path(coverage.dir, sub(".bam$", "", bam.file))
    if (! file.exists(paste0(coverage.path, ".sample_summary"))) {
        command <- sprintf("java -Xmx1024m -jar %s -T DepthOfCoverage -omitBaseOutput -omitLocusTable -R %s -I %s -L %s -o %s",
                           gatk.jar, human.fasta, bam.path, exome.region.list, coverage.path)
        system(command)
    } else {
        cat("Output file", coverage.path, "already exists")
    }
})
