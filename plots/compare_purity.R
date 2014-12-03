#!/usr/bin/env Rscript

library(ggplot2)
library(GGally)

absolute <- read.delim("../ABSOLUTE/summary.PP-calls_tab.txt")
absolute <- absolute[,c("sample", "purity")]
colnames(absolute) <- c("sample", "ABSOLUTE")

pathologist <- read.delim("../metadata_colorectal.csv")
pathologist$purity <- pathologist$purity/100
pathologist$sample <- pathologist$tumor.sample
pathologist <- pathologist[,c("sample", "purity")]
colnames(pathologist) <- c("sample", "pathologist")

maf.files <- list.files("../MAF/by-sample", full.names=T)
mafs <- lapply(maf.files, read.delim)
names(mafs) <- basename(maf.files)
vafs <- do.call(rbind, lapply(names(mafs), function (m) {
    maf <- mafs[[m]]
    if (nrow(maf) == 0) return (NULL)
    data.frame(sample=strsplit(m, ".", fixed=T)[[1]][1],
               vaf=maf$t_alt_count*2/(maf$t_ref_count+maf$t_alt_count))
}))
peak.vaf <- aggregate(vaf~sample, vafs, function (x) {
    dens <- density(x)
    dens$x[which.max(dens$y)]
})
colnames(peak.vaf) <- c("sample", "peak.VAF")

titan <- read.table("../TITAN/purity.tsv")
colnames(titan) <- c("sample", "TITAN")

purity <- Reduce(merge, list(absolute, pathologist, titan, peak.vaf))
pdf("compare-purity.pdf")
ggpairs(purity[,2:5]) 
dev.off()
