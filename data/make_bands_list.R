#!/usr/bin/env Rscript

# make a list of locations of chromosome bands

band.url <- paste0("http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/",
                   "cytoBand.txt.gz")
con <- textConnection(readLines(gzcon(url(band.url))))
cyto.band <- read.table(con, col.names=c("chr", "start", "end", "locus", "g"))

cyto.band$chr <- sub("chr", "", cyto.band$chr)
cyto.band$band <- sub("[.].*", "", cyto.band$locus)

# aggregate sub-bands into bands
band.starts <- aggregate(start~chr+band, cyto.band, min)
band.ends <- aggregate(end~chr+band, cyto.band, max)
bands <- merge(band.starts, band.ends)
write.table(bands, "bands.tsv", sep="\t", row.names=F, quote=F)
