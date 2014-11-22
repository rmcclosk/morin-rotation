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
classes <- c(rep("NULL", 2), "character", rep("NULL", 3), rep("integer", 2), 
             rep("NULL", 4), "character", rep("NULL", 3))
names <- c(rep("NULL", 2), "chr", rep("NULL", 3), "start", "end",
           rep("NULL", 4), "name", rep("NULL", 3))

# stackoverflow.com/questions/9548630/read-gzipped-csv-directly-from-a-url-in-r
con <- textConnection(readLines(gzcon(url(loci.url))))
gene.loci <- read.table(con, colClasses=classes, col.names=names)
gene.loci <- merge(gene.loci, gene.names)
gene.loci <- merge(gene.loci, aggregate(start~name, gene.loci, min))
gene.loci <- merge(gene.loci, aggregate(end~name, gene.loci, max))
gene.loci <- unique(gene.loci)
gene.loci$chr <- sub("chr", "", gene.loci$chr)

write.table(gene.loci, file="genes.tsv", row.names=F, quote=F, sep="\t")
