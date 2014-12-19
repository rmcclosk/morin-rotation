#!/usr/bin/env Rscript

d <- read.delim("variants.tsv", stringsAsFactors=F)
d <- subset(d, grepl("049", sample))
d <- d[order(d$sample, d$chr, d$start),]

d$seq <- ifelse(d$alt.count > 2 & d$depth > 20, d$alt, d$ref)
seqs <- aggregate(seq~sample, d, paste, collapse="")
fasta <- paste0(">", seqs$sample, "\n", seqs$seq, "\n")

tf <- tempfile()
cat(fasta, sep="", file=tf)
cat("Aligning... ")
align <- system(paste("mafft --auto", tf), intern=T)
cat("done\n")
align
cat(align, sep="\n", file=tf)
system("rm RAxML*")
cat("Tree building... ")
system(paste("raxmlHPC -V -m GTRCAT -n 049 -p 42 -s", tf))
cat("done\n")
