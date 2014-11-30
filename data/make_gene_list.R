#!/usr/bin/env Rscript

# Fetch a list of genes (HUGO symbols) and coordinates.

url.options <- paste("col=gd_app_sym",
                     "status=Approved",
                     "status_opt=2", 
                     paste0(paste0("chr=", c(1:22, "X", "Y")), collapse="&"),
                     "order_by=gd_hgnc_id",
                     "format=text",
                     "submit=submit",
                     sep="&")
base.url <- "http://www.genenames.org/cgi-bin/download?"
names.url <- paste0(base.url, url.options)
gene.names <- read.table(names.url, header=T, sep="\t", col.names=c("name"))

loci.url <- paste0("http://hgdownload.cse.ucsc.edu/goldenPath/",
                   "hg19/database/refGene.txt.gz")
classes <- c(rep("NULL", 2), "character", "NULL", rep("integer", 2), 
             rep("NULL", 6), "character", rep("NULL", 3))
names <- c(rep("NULL", 2), "chr", "NULL", "start", "end",
           rep("NULL", 6), "name", rep("NULL", 3))

# stackoverflow.com/questions/9548630/read-gzipped-csv-directly-from-a-url-in-r
con <- textConnection(readLines(gzcon(url(loci.url))))
gene.loci <- read.table(con, colClasses=classes, col.names=names)
gene.loci <- merge(gene.loci, gene.names)
gene.loci$chr <- factor(sub("chr", "", gene.loci$chr), levels=c(1:22, "X", "Y"))
gene.loci <- merge(aggregate(start~chr+name, gene.loci, min),
                   aggregate(end~chr+name, gene.loci, max))
gene.loci <- gene.loci[,c("chr", "start", "end", "name")]
gene.loci <- gene.loci[order(gene.loci$chr, gene.loci$start),]

write.table(gene.loci, file="genes.bed", col.names=F, row.names=F, quote=F, sep="\t")
