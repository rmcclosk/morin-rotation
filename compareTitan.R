#!/usr/bin/env Rscript

library(parallel)

source(file="/home/rmccloskey/morin-rotation/settings.conf")
source(file="/home/rmccloskey/morin-rotation/merge.R")
titan.dir <- file.path(WORK_DIR, "14_titan")

seg.classes <- c("NULL", "factor", rep("numeric", 2), rep("NULL", 5), "integer",
                 rep("NULL", 4))
seg.names <- c("chr", "start", "end", "copy.number")
read.segs <- function (f) {
    segs <- read.table(f, header=T, colClasses=seg.classes)
    colnames(segs) <- seg.names
    segs
}

dist.segs <- function (segs1, segs2) {
    colnames(segs1)[4] <- "copy.num.1"
    colnames(segs2)[4] <- "copy.num.2"
    merged <- do.call(rbind, lapply(c(1:22, "X", "Y"), function (by.chr) {
        s1 <- subset(segs1, chr==by.chr)[,2:4]
        s2 <- subset(segs2, chr==by.chr)[,2:4]
        if (nrow(s1) == 0 | nrow(s2) == 0)
            return (NULL)
        s <- merge.segs(s1, s2)
        s$chr <- by.chr
        s
    }))
    merged <- merged[!is.na(merged$copy.num.1) & !is.na(merged$copy.num.2),]
    merged$start <- as.numeric(merged$start)
    merged$end <- as.numeric(merged$end)
    merged$copy.num.1 <- as.integer(merged$copy.num.1)
    merged$copy.num.2 <- as.integer(merged$copy.num.2)
    with(merged, sqrt(sum((end-start+1)*abs(copy.num.1-copy.num.2)^2)))
}

d <- read.table(METADATA, header=T)
dirs <- file.path(titan.dir, d$tumor.sample)
d <- d[file.exists(dirs),]
files <- Sys.glob(file.path(dirs, "titan", "*cluster_*_segs.txt"))
clusters <- data.frame(file=files)
dirs <- gsub(paste0(titan.dir, .Platform$file.sep), "", files)
clusters$sample <- sapply(strsplit(dirs, "/"), "[[", 1)
patients <- sapply(strsplit(clusters$sample, "-"), "[[", 3)
clusters$patient <- paste0("01-", patients)

all.segs <- lapply(files, read.segs)
names(all.segs) <- clusters$file

pairs <- expand.grid(clusters$file, clusters$file)
pairs <- pairs[1:2,]
pairs <- merge(pairs, clusters, by.x=c("Var1"), by.y=c("file"))
colnames(pairs)[colnames(pairs) == "sample"] <- "sample.1"
colnames(pairs)[colnames(pairs) == "patient"] <- "patient.1"
colnames(pairs)[colnames(pairs) == "Var1"] <- "file.1"
pairs <- merge(pairs, clusters, by.x=c("Var2"), by.y=c("file"))
colnames(pairs)[colnames(pairs) == "sample"] <- "sample.2"
colnames(pairs)[colnames(pairs) == "patient"] <- "patient.2"
colnames(pairs)[colnames(pairs) == "Var2"] <- "file.2"

pairs$dist <- apply(pairs, 1, function (row) {
    dist.segs(all.segs[[row["file.1"]]], all.segs[[row["file.2"]]])
})

pairs <- pairs[,c("sample.1", "patient.1", "sample.2", "patient.2",
                  "dist")]
write.table(pairs, file="titan_pairs.csv", col.names=T, row.names=F, quote=F)
