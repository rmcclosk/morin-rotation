#!/usr/bin/env Rscript

sample.data <- read.table("sample_info.dat", header=T)
sample.data$purity <- floor(sample.data$purity/0.05)*0.05
normal.data <- sample.data[is.na(sample.data$purity),]
tumor.data <- sample.data[!is.na(sample.data$purity),]
sample.data <- merge(normal.data, tumor.data, by=c("patient.id"), 
                     suffixes=c(".normal", ".tumor"))
sample.data$purity.tumor <- 0.5

seg.data <- read.table("annotations.dat", header=T)
seg.data <- merge(seg.data, sample.data, by.x=c("sample", "purity"), by.y=c("sample.tumor", "purity.tumor"))
seg.data$arm <- substr(seg.data$locus, 1, 1)
seg.data <- subset(seg.data, (! chr %in% c(13, 14, 15, 21, 22)) | arm != "p")
seg.data <- subset(seg.data, chr != "Y")
seg.data$length <- seg.data$end-seg.data$start+1
seg.data$cnv.exome <- (seg.data$exome.copy != 2)*seg.data$length
seg.data$cnv.genome <- (seg.data$genome.copy != 2)*seg.data$length

total.length <- aggregate(length~sample, seg.data, sum)
exome.cnv.length <- aggregate(cnv.exome~sample, seg.data, sum)
genome.cnv.length <- aggregate(cnv.genome~sample, seg.data, sum)
cnv.data <- merge(total.length, exome.cnv.length, by=c("sample"))
cnv.data <- merge(cnv.data, genome.cnv.length, by=c("sample"))
cnv.data$cnv.genome <- cnv.data$cnv.genome/cnv.data$length
cnv.data$cnv.exome <- cnv.data$cnv.exome/cnv.data$length
cnv.data <- cnv.data[,c("sample", "cnv.genome", "cnv.exome")]

sample.data <- merge(sample.data, cnv.data, by.x=c("sample.tumor"),
                     by.y=c("sample"))
sample.data$coverage.ratio.exome <- sample.data$coverage.mean.exome.tumor/
                                    sample.data$coverage.mean.exome.normal
sample.data$coverage.ratio.genome <- sample.data$coverage.mean.genome.tumor/
                                     sample.data$coverage.mean.genome.normal
sample.data <- sample.data[order(sample.data$cnv.genome),]
sample.data
quit()

fit <- glm(cnv.exome~prop.unmapped.exome.tumor +
                     prop.dup.exome.tumor +
                     prop.fail.exome.tumor +
                     coverage.mean.exome.tumor +
                     coverage.sd.exome.tumor +
                     coverage.ratio.exome +
                     purity.tumor,
           data=sample.data)
summary(fit)

fit <- glm(cnv.genome~prop.unmapped.genome.tumor +
                      prop.dup.genome.tumor +
                      prop.fail.genome.tumor +
                      coverage.mean.genome.tumor +
                      coverage.sd.genome.tumor +
                      coverage.ratio.genome +
                      purity.tumor,
           data=sample.data)
summary(fit)

