#!/usr/bin/env Rscript

chr.data <- read.table("cytoband.txt", 
                       col.names=c("chr", "start", "end", "locus", "g"))
chr.data$chr <- factor(sub("chr", "", chr.data$chr))
chr.data$band <- sub("[.].*", "", chr.data$locus)

# aggregate sub-bands into bands
chr.data <- do.call(rbind, by(chr.data, list(chr.data$chr, chr.data$band), function (x) {
    data.frame(chr=x[1,"chr"], 
               start=min(x$start), 
               end=max(x$end), 
               locus=x[1,"band"])
}, simplify=F))

seg.data <- read.table("annotations.dat", header=T)
seg.data$length <- seg.data$end - seg.data$start + 1

# use only the purity closest to the pathologist's estimate
purity.to.use <- do.call(rbind, by(seg.data, seg.data$sample, function (x) {
    x[which.min(abs(x$purity-x$report.purity)), c("sample", "purity")]
}, simplify=F))

# do each sample separately
first <- T
. <- by(seg.data, seg.data$sample, function (sample.seg.data) {
    
    cat("Processing sample", sample.seg.data[1, "sample"], "...")
    # loop through loci (should loop through chr to be faster)
    copy.by.chr <- do.call(rbind, lapply(1:nrow(chr.data), function (i) {
        by.chr <- chr.data[i, "chr"]
        chr.start <- chr.data[i, "start"]
        chr.end <- chr.data[i, "end"]
        locus <- chr.data[i, "locus"]
    
        # Find genome intervals overlapping this exome interval.
        overlap <- subset(sample.seg.data, 
                          chr==by.chr & 
                          start <= chr.end & 
                          end >= chr.start)
    
        if (nrow(overlap) == 0) return (NULL)
    
        # loop through overlapping genome segments
        overlap <- do.call(rbind, lapply(1:nrow(overlap), function (j) {
            seg.start <- overlap[j,"start"]
            seg.end <- overlap[j,"end"]
            
            start <- ifelse(seg.start <= chr.start, chr.start, seg.start)
            end <- ifelse(seg.end >= chr.end, chr.end, seg.end)
            c(length=end-start+1, genome.copy=overlap[j, "genome.copy"], 
              exome.copy=overlap[j,"exome.copy"])
        }))

        # take state composing the longest part of the segment
        exome.agg <- aggregate(length~exome.copy, overlap, sum)
        genome.agg <- aggregate(length~genome.copy, overlap, sum)
        exome.copy <- exome.agg[which.max(exome.agg$length), "exome.copy"]
        genome.copy <- genome.agg[which.max(genome.agg$length), "genome.copy"]
        data.frame(chr=by.chr, locus=locus, exome.copy=exome.copy, genome.copy=genome.copy)
    }))
    cat(" done\n")

    copy.by.chr$sample <- seg.data[1, "sample"]

    if (first) {
        write.table(copy.by.chr, file="regions.dat", col.names=T, quote=F, row.names=F)
        first <<- F
    } else {
        write.table(copy.by.chr, file="regions.dat", col.names=F, quote=F, row.names=F, append=T)
    }
})
