#!/usr/bin/env Rscript

chrs <- c(1:22, "X")

diep <- read.table("../data/diep.tsv", header=T, stringsAsFactors=F)
diep$chr <- factor(sub("[pq].*", "", diep$locus), levels=chrs)
diep$band <- sub("^[0-9X]+", "", diep$locus)

bed.files <- file.path("beds", list.files("beds"))

do.intersect <- function (bed.file.1, bed.file.2, data.cols) {
    cmd <- paste("bedtools intersect -a", bed.file.1, "-b", bed.file.2)
    cat(paste0(cmd, "\n"))
    res <- read.table(textConnection(system(cmd, intern=T)))
    colnames(res)[1:3] <- c("chr", "start", "end")
    colnames(res)[4:ncol(res)] <- data.cols
    res
}

bands <- do.call(rbind, lapply(bed.files, function (f) {
    res <- merge(do.intersect(f, "bands.bed", "copy.number"),
                 do.intersect("bands.bed", f, "band"))
    res <- subset(res, end - start > 1000)
    dups <- aggregate(copy.number~band, res, function (x) any(x > 2))
    colnames(dups)[2] <- "dup"
    dels <- aggregate(copy.number~band, res, function (x) any(x < 2))
    colnames(dels)[2] <- "del"
    res <- merge(dups, dels)
    res$sample <- strsplit(basename(f), ".", fixed=T)[[1]][1]
    if (res$sample == "R3")
        res$patient <- "cell-line"
    else
        res$patient <- paste0("01-", strsplit(basename(f), "-")[[1]][3])
    res
}))
write.table(bands, "tmp.dat", row.names=F, col.names=F)

bands <- read.table("tmp.dat")
colnames(bands) <- c("band", "dup", "del", "sample", "patient")

dups <- aggregate(dup~band+patient, bands, any)
dels <- aggregate(del~band+patient, bands, any)
cnas <- merge(dups, dels)

dups <- aggregate(dup~band, cnas, sum)
dels <- aggregate(del~band, cnas, sum)
cnas <- merge(dups, dels)

cnas
