#!/usr/bin/env Rscript

source(file="/home/rmccloskey/morin-rotation/settings.conf")

out.dir <- file.path(WORK_DIR, "11_vaf")
dir.create(out.dir, showWarnings=F)

d <- read.csv(EXOME_METADATA, header=T)
d <- d[,c("patient_id", "sample_id", "gsc_exome_library_id")]
d$strelka.id <- paste(d$sample_id, d$gsc_exome_library_id, sep="-")
d$out.file <- paste0(file.path(out.dir, d$sample_id), ".dat")

norm <- d[grepl("(BF|WB)", d$sample_id),]
tum <- d[!grepl("(BF|WB)", d$sample_id),]
d <- merge(norm, tum, by=c("patient_id"), suffixes=c(".normal", ".tumor"))

d$strelka.dir <- paste0(d$strelka.id.normal, "_vs_", d$strelka.id.tumor)
vcf.file <- file.path("results", "passed.somatic.snvs.vcf")
d$strelka.file <- file.path(STRELKA_PATH, d$strelka.dir, vcf.file)

vcf.cols <- c("chr", "pos", "id", "ref", "alt", "qual", "filter", "info", 
              "format", "normal", "tumor")
nts <- c("A", "C", "G", "T")

. <- apply(d, 1, function (row) {
    cat("Processing sample", row["sample_id.tumor"], "...")
    vcf <- read.table(row["strelka.file"], col.names=vcf.cols, stringsAsFactors=F)
    counts <- strsplit(vcf$tumor, ":")
    read.depth <- as.integer(sapply(counts, "[[", 1))

    # sum counts for all alternate alleles
    alt.idx <- sapply(sapply(strsplit(vcf$alt, ","), match, nts), "+", 4)
    var.counts <- mapply("[", counts, alt.idx, SIMPLIFY=F)
    var.counts <- rapply(var.counts, strsplit, classes="character", 
                         how="replace", split=",")
    var.counts <- sapply(var.counts, function (l) sapply(l, "[[", 1))
    var.counts <- sapply(var.counts, as.integer)
    var.counts <- sapply(var.counts, sum)

    vcf$vaf <- var.counts/read.depth
    vcf <- vcf[,c("chr", "pos", "ref", "alt", "vaf")]
    write.table(vcf, file=row["out.file.tumor"], col.names=T, row.names=F)
    cat(" done\n")
})
