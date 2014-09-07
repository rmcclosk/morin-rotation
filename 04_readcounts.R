source(file="settings.conf")

out.dir <- file.path(GENOME_WORK_DIR, "04_readcounts")
dir.create(out.dir, showWarnings=F)
matching <- read.table(GENOME_MATCHING, stringsAsFactors=F, fill=T, sep="\t")
files <- c(matching$V1, matching$V2, matching$V3)
files <- files[!files %in% c("", "Insufficient DNA")]
files <- sapply(files, function (f) strsplit(f, " ")[[1]][2])

main <- sapply(files, function (bam.file) {
    bam.path <- file.path(GENOME_BAM_DIR, bam.file)
    out.path <- file.path(out.dir, sub(".bam$", ".wig", basename(bam.file)))
    if (! file.exists(out.path)) {
        command <- paste(file.path(HMMCOPY_DIR, "readCounter"), bam.path,
                         ">", out.path)
        cat(command, "\n")
    } else {
        cat("Output file", out.path, "already exists\n", file=stderr())
    }
})
