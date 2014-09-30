d <- read.csv("freqs.dat", header=T)
d <- d[d$chrom %in% c(1:22, "X", "Y"),]
d$chrom <- factor(d$chrom, levels=c(1:22, "X", "Y"))
