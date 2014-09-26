#!/usr/bin/env Rscript

source(file="/home/rmccloskey/morin-rotation/settings.conf")

d <- read.csv(GENOME_METADATA, header=T)[,c("sample_id", "patient_id")]
d <- d[!grepl("(WB|BF)", d$sample_id),]

hmmcopy.dir <- file.path(WORK_DIR, "05_hmmcopy", "10000", d$patient_id) 
d$seg.file <- file.path(hmmcopy.dir, paste0(d$sample_id, "_default_segments.dat"))
d <- d[file.exists(d$seg.file),]

# average copy number across genome
d$avg.copy <- sapply(d$seg.file, function (f) {
    seg <- read.table(f, header=T)
    seg$state <- as.integer(sub("\"", "", seg$state))
    seg$length <- seg$end-seg$start+1
    seg$ampl <- seg$state > 3
    seg$del <- seg$state < 3
    sum(seg$length*(seg$state-1))/sum(seg$length)
    #sum(seg$del*seg$length)/sum(seg$length)
    #sum(seg$del+seg$ampl)/nrow(seg)
})

sample.data <- read.csv("BC Genome Profiling-HQC.csv", header=T)
colnames(sample.data)
purity <- sample.data$X..tumor/100 *
          (sample.data$X..viable.neoplastic.cells + sample.data$X..necrosis)
sample.data$purity <- ifelse(is.na(purity), sample.data$X..Tumor, purity)
sample.data <- sample.data[,c("Sample.ID", "purity")]

d <- merge(d, sample.data, by.x=c("sample_id"), by.y=c("Sample.ID"))

cor.test(d$avg.copy, d$purity)

pdf("purityVsCNV.pdf")
plot(d$purity, d$avg.copy, main="Average copy number vs. tumor purity",
     xlab="estimated purity", ylab="average copy number")
dev.off()
