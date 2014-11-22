#!/usr/bin/env Rscript

library(ggplot2)

source(file="/home/rmccloskey/morin-rotation/merge.R")

round0.05 <- function (x) {
    res <- gsub("0$", "", as.character(floor(x/0.05)*0.05))
    ifelse(res == "", 0, res)
}

seg.error <- function (merged.segs) {
    diff <- abs(merged.segs$copy.number.titan - merged.segs$copy.number.exomecnv)
    sqrt(sum(diff*(merged.segs$end-merged.segs$start))^2)
}

read.exomecnv.segs <- function (admix, sample) {
    exomecnv.dir <- "/extscratch/morinlab/shared/rmccloskey/colorectal_exome/03_cnv"
    seg.file <- paste0(admix, ".segment.copynumber.txt.bz2")
    seg.file <- file.path(exomecnv.dir, sample, seg.file)
    exomecnv.segs <- read.table(bzfile(seg.file), col.names=c("chr", "start", "end", "copy.number.exomecnv"))
    exomecnv.segs$chr <- factor(sub("chr", "", exomecnv.segs$chr), levels=c(1:22, "X", "Y"))
    exomecnv.segs
}

# read purity estimates from pathologist
d <- read.table("/home/rmccloskey/morin-rotation/metadata_colorectal.csv", header=T)
d$admix.path <- 1-d$purity/100

# read titan segments
titan.dir <- "/extscratch/morinlab/shared/rmccloskey/colorectal_exome/14_titan"
d <- d[file.exists(file.path(titan.dir, d$tumor.sample)),]
d$best.cluster <- sapply(d$tumor.sample, function (s) {
    files <- paste0("*cluster_", 1:5, "_params.txt.bz2")
    files <- Sys.glob(file.path(titan.dir, s, "titan", files))
    which.min(sapply(files, function (f) {
        params <- read.table(bzfile(f), sep="\t", stringsAsFactors=F)
        as.numeric(subset(params, V1 == "S_Dbw dens.bw:", select=V2)[1])
    }))
})

segs <- do.call(rbind, lapply(1:nrow(d), function (i) {
    sample <- d[i, "tumor.sample"]
    best.cluster <- d[i, "best.cluster"]

    #param.file <- paste0("*cluster_", best.cluster, "_params.txt.bz2")
    #param.file <- Sys.glob(file.path(titan.dir, sample, "titan", param.file))
    #params <- read.table(bzfile(param.file), sep="\t", stringsAsFactors=F)
    #titan.admix <- as.numeric(subset(params, V1=="Normal contamination estimate:", select=V2)[1])

    seg.file <- paste0("*cluster_", best.cluster, "_segs.txt.bz2")
    seg.file <- Sys.glob(file.path(titan.dir, sample, "titan", seg.file))[1]
    titan.segs <- read.table(bzfile(seg.file), header=T,
                             colClasses=c("NULL", "factor", rep("integer", 2), rep("NULL", 5), "integer", rep("NULL", 4)),
                             col.names=c("a", "chr", "start", "end", "b", "c", "d", "e", "f", "copy.number.titan", "g", "h", "i", "j"))
    colnames(titan.segs)[colnames(titan.segs) == "copy.number.titan"] <- "copy.number"
    titan.segs <- cbind(sample=sample, titan.segs)
    return(titan.segs)
    
    best.idx <- which.min(sapply(seq(0, 0.85, 0.05), function (admix) {
        exomecnv.segs <- read.exomecnv.segs(admix, sample)
        
        merged.segs <- do.call(rbind, lapply(c(1:22, "X"), function (by.chr) {
            merged.segs <- merge.segs(subset(exomecnv.segs, chr==by.chr, select=c("start", "end", "copy.number.exomecnv")),
                                      subset(titan.segs, chr==by.chr, select=c("start", "end", "copy.number.titan")))
            merged.segs <- subset(merged.segs, !is.na(merged.segs$copy.number.exomecnv) & !is.na(merged.segs$copy.number.titan))
            merged.segs$start <- as.numeric(merged.segs$start)
            merged.segs$end <- as.numeric(merged.segs$end)
            merged.segs$copy.number.exomecnv <- as.integer(merged.segs$copy.number.exomecnv)
            merged.segs$copy.number.titan <- as.integer(merged.segs$copy.number.titan)
            merged.segs$chr <- by.chr
            merged.segs
        }))
        seg.error(merged.segs)
    }))
    best.admix <- 0.05*(best.idx-1)

    exomecnv.segs <- read.exomecnv.segs(best.admix, sample)
    colnames(exomecnv.segs)[colnames(exomecnv.segs) == "copy.number.exomecnv"] <- "copy.number"
    exomecnv.segs$method <- "ExomeCNV"
    titan.segs$method <- "TITAN"
    all.segs <- rbind(exomecnv.segs, titan.segs)
    
    chrs <- c(1:22, "X")
    
    pdf(file.path("/home/rmccloskey/morin-rotation/exome/compare-plots", paste0(sample, ".pdf")))
    sapply(seq(1, length(chrs), 4), function (i) {
        plot.segs <- droplevels(subset(all.segs, chr %in% chrs[i:(i+3)]))
        print(ggplot(plot.segs, aes(x=start, y=copy.number, color=method)) +
            geom_segment(aes(xend=end, yend=copy.number)) +
            facet_wrap(~chr))
    })
    dev.off()
}))

write.table(segs, "segments.tsv", sep="\t", row.names=F, quote=F)
