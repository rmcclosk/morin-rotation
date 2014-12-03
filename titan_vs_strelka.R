#!/usr/bin/env Rscript

strelka.dir <- "/extscratch/morinlab/shared/rmccloskey/colorectal_exome/11_snv"
titan.dir <- "/home/rmccloskey/morin-rotation/TITAN/titan-exome"

titan.files <- Sys.glob(file.path(titan.dir, "*_titan.txt.bz2"))
titan.file <- titan.files[2]
sample <- strsplit(basename(titan.file), "_")[[1]][1]
strelka.file <- file.path(strelka.dir, paste0(sample, ".maf"))

titan <- read.table(bzfile(titan.file), header=T, sep="\t")
strelka <- read.table(strelka.file, header=T, sep="\t")

chrs <- c(1:22, "X")
titan <- subset(titan, Chr %in% chrs, select=c("Chr", "Position"))
strelka <- subset(strelka, Chromosome %in% chrs, select=c("Chromosome", "Start_Position"))
colnames(strelka) <- c("chr", "pos")
colnames(titan) <- c("chr", "pos")
strelka$chr <- factor(strelka$chr, levels=chrs)
titan$chr <- factor(titan$chr, levels=chrs)

merge(titan, strelka)

titan$pos <- titan$pos + 1
merge(titan, strelka)

titan$pos <- titan$pos - 2
merge(titan, strelka)
