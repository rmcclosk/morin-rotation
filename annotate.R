#!/usr/bin/env Rscript

source(file="settings.conf")

purity <- read.csv(EXOME_PURITY, header=T, fill=T, stringsAsFactors=F)

colnames(purity)[1] <- "sample"
purity$patient.id <- paste0("01-", strsplit(purity$sample, "-")[[1]][3])
purity
    
exome.dir <- file.path(WORK_DIR, "03_cnv")
genome.dir <- file.path(WORK_DIR, "05_hmmcopy")

sapply(

quit()

exome.seg <- do.call(rbind, lapply(seq(0.1, 0.85, 0.05), function (admixture.rate) {
    exome.file <- file.path(exome.dir, 
                            paste0(admixture.rate, ".segment.copynumber.txt"))
    d <- read.table(exome.file, 
                    col.names=c("chr", "start", "end", "copy.number"),
                    colClasses=c("character", rep("numeric", 3)))
    d$chr <- factor(sub("chr", "", d$chr))
    d$purity <- 1-admixture.rate
    d
}))

genome.file <- file.path(genome.dir, "1409141605_segments.dat")
genome.seg <- read.table(genome.file, header=T)
genome.seg$copy.number <- genome.seg$state - 1

# by chromosome and purity
all.seg <- do.call(rbind, by(exome.seg, list(exome.seg$chr, exome.seg$purity), function (exome.chr.seg) {
    by.chr <- exome.chr.seg[1, "chr"]
    by.purity <- exome.chr.seg[1, "purity"]

    # loop through segments (rows)
    d <- do.call(rbind, lapply(1:nrow(exome.chr.seg), function (i) {
        exome.start <- exome.chr.seg[i,"start"]
        exome.end <- exome.chr.seg[i,"end"]
        exome.copy <- exome.chr.seg[i, "copy.number"]
        
        # Find genome intervals overlapping this exome interval.
        overlap <- subset(genome.seg, 
                          chr==by.chr & start <= exome.end & end >= exome.start)
    
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
    d
}, simplify=F))

write.table(all.seg, file="annotations.dat", row.names=F, col.names=T, quote=F)








