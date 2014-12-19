#!/usr/bin/env Rscript

# make a list of locations of chromosome bands

band.url <- paste0("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/",
                   "cytoBand.txt.gz")
con <- textConnection(readLines(gzcon(url(band.url))))
cyto.band <- read.table(con, col.names=c("chr", "start", "end", "locus", "g"))

cyto.band$chr <- factor(sub("chr", "", cyto.band$chr), levels=c(1:22, "X", "Y"))
cyto.band$band <- sub("[.].*", "", cyto.band$locus)

# aggregate sub-bands into bands
band.starts <- aggregate(start~chr+band, cyto.band, min)
band.ends <- aggregate(end~chr+band, cyto.band, max)
bands <- merge(band.starts, band.ends)
bands <- bands[order(bands$chr, bands$start), c("chr", "start", "end", "band")]
write.table(bands, "bands.bed", sep="\t", col.names=F, row.names=F, quote=F)
