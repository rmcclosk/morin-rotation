#!/usr/bin/env Rscript

source(file="settings.conf")

WORK_DIR <- sub("dlbcl", "colorectal", WORK_DIR)

###############################################################################
# Functions to read and process data files
###############################################################################

# Read a particular exome segments file.
read.exome.segments <- function (file.path) {
    d <- read.table(file.path,
                    col.names=c("chr", "start", "end", "copy.number"),
                    colClasses=c("character", rep("numeric", 3)))
    d$chr <- factor(sub("chr", "", d$chr))
    admixture.rate <- as.numeric(sub(".segment.copynumber.txt", "", 
                                     basename(file.path)))
    d$purity <- 1-admixture.rate
    d
}

# Read all exome segments files for a particular patient/sample.
read.all.exome <- function (sample, patient.id, exome.dir) {
    cat("Reading exome segments for sample", sample, "... ")
    seg.dir <- file.path(exome.dir, patient.id, sample)
    seg.files <- Sys.glob(file.path(seg.dir, "*segment.copynumber.txt"))
    d <- do.call(rbind, lapply(seg.files, read.exome.segments))
    d$patient.id <- patient.id
    d$sample <- sample
    cat("done\n")
    d
}

# Read a particular genome segments file.
read.genome.segments <- function (file.path) {
    d <- read.table(file.path, header=T)
    d$copy.number <- d$state - 1
    d <- d[,c("chr", "start", "end", "copy.number")]
}

# Read all genome segments files for a particular patient/sample.
read.all.genome <- function (sample, patient.id, genome.dir) {
    cat("Reading genome segments for sample", sample, "... ")
    seg.dir <- file.path(genome.dir, patient.id)
    seg.file <- file.path(seg.dir, paste0(sample, "_segments.dat"))
    if (!file.exists(seg.file)) {
        cat("none found\n")
        return (NULL)
    }
    d <- read.genome.segments(seg.file)
    d$patient.id <- patient.id
    d$sample <- sample
    cat("done\n")
    d
}

###############################################################################
# Main
###############################################################################

# Read in sample annotations
sample.data <- read.csv(EXOME_PURITY, header=T, fill=T, stringsAsFactors=F)
colnames(sample.data)[1] <- "sample"
sample.data <- sample.data[!sample.data$sample %in% c("HT 29", "R3"),]

sample.data$patient.id <- paste0("01-", sapply(strsplit(sample.data$sample, "-"), "[[", 3))
sample.data$purity <- ifelse(is.na(sample.data$X..viable.neoplastic.cells),
                             sample.data$X..Tumor,
                             sample.data$X..viable.neoplastic.cells)
sample.data <- sample.data[,c("patient.id", "sample", "purity")]

exome.dir <- file.path(WORK_DIR, "03_cnv")
genome.dir <- file.path(WORK_DIR, "05_hmmcopy")

exome.seg <- do.call(rbind, mapply(read.all.exome, 
                                   sample.data$sample, 
                                   sample.data$patient.id, 
                                   exome.dir, 
                                   SIMPLIFY=FALSE))
genome.seg <- do.call(rbind, mapply(read.all.genome, 
                                    sample.data$sample, 
                                    sample.data$patient.id, 
                                    genome.dir, 
                                    SIMPLIFY=FALSE))
exome.seg <- exome.seg[exome.seg$sample %in% genome.seg$sample,]

first <- T
# by chromosome and purity
. <- by(exome.seg, list(exome.seg$sample, exome.seg$chr, exome.seg$purity), function (exome.chr.seg) {
    by.chr <- exome.chr.seg[1, "chr"]
    by.purity <- exome.chr.seg[1, "purity"]
    by.sample <- exome.chr.seg[1, "sample"]
    by.patient.id <- exome.chr.seg[1, "patient.id"]
    
    cat("Processing sample", by.sample, ", chromosome", by.chr, ", purity", by.purity, "...")

    # loop through segments (rows)
    d <- do.call(rbind, lapply(1:nrow(exome.chr.seg), function (i) {
        exome.start <- exome.chr.seg[i,"start"]
        exome.end <- exome.chr.seg[i,"end"]
        exome.copy <- exome.chr.seg[i, "copy.number"]
        
        # Find genome intervals overlapping this exome interval.
        overlap <- subset(genome.seg, sample==by.sample & chr==by.chr & 
                          start <= exome.end & end >= exome.start)
    
        # loop through overlapping genome segments
        do.call(rbind, lapply(1:nrow(overlap), function (j) {
            genome.start <- overlap[j,"start"]
            genome.end <- overlap[j,"end"]
            genome.copy <- overlap[j, "copy.number"]
            
            start <- ifelse(genome.start <= exome.start, exome.start, genome.start)
            end <- ifelse(genome.end >= exome.end, exome.end, genome.end)
            c(start=start, end=end, genome.copy=genome.copy, exome.copy=exome.copy)
        }))
    }))
    d <- as.data.frame(d)
    d$purity <- by.purity
    d$chr <- by.chr
    d$patient.id <- by.patient.id
    d$sample <- by.sample
    cat(" done\n")
    if (first) {
        write.table(d, file="annotations.dat", col.names=T, row.names=F, quote=F)
        first <<- F
    } else {
        write.table(d, file="annotations.dat", append=T, row.names=F, col.names=F, quote=F)
    }
})
