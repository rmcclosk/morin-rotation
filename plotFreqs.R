#!/usr/bin/env Rscript

d <- read.csv("freqs.dat", header=T)
by(d, d$patient, function (d.patient) {
})

#pdf("freqs.pdf")

#dev.off()
