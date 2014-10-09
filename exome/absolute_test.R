#!/usr/bin/env Rscript

library(ABSOLUTE)
library(foreach)
library(doMC)

source(file="/home/rmccloskey/morin-rotation/settings.conf")
options("scipen"=100)

coverage.dir <- file.path(WORK_DIR, "02_coverage")
cnv.dir <- file.path(WORK_DIR, "03_cnv")
maf.dir <- file.path(WORK_DIR, "11_snv")
vcf.dir <- file.path(WORK_DIR, "12_vcf")
out.dir <- file.path(WORK_DIR, "13_absolute")
dir.create(out.dir, showWarnings=F)

platform <- "Illumina_WES"
primary.disease <- "colorectal_cancer"
sigma.p <- 0.01
max.sigma.h <- 0.03
min.ploidy <- 1
max.ploidy <- 6
max.as.seg.count <- 1500
max.non.clonal <- 0.1
max.neg.genome <- 0.01
min.mut.af <- 0.1
copy_num_type <- "total"

DoAbsolute <- function(sample) {
    registerDoSEQ()
  
    # read segments 
    seg.file <- file.path(cnv.dir, sample, "0.segment.lrr.txt")
    seg <- read.table(seg.file, header=T,
                      colClasses=c("character", rep("numeric", 2), "numeric"),
                      col.names=c("Chromosome", "Start", "End", "Segment_Mean"))
    seg$Chromosome <- sub("chr", "", seg$Chromosome)
    seg$Num_Probes <- floor((seg$End-seg$Start)/1000)
    
    seg$Chromosome <- factor(seg$Chromosome, levels=c(1:22, "X", "Y"))
    seg <- seg[order(seg$Chromosome, seg$Start),]
    seg.file <- file.path(out.dir, paste0(sample, ".seg"))
    write.table(seg, file=seg.file, col.names=T, row.names=F, quote=F, sep="\t")
    read.maf.vcf(sample)

    maf.file <- file.path(out.dir, paste0(sample, ".maf"))
    log.file <- file.path(out.dir, paste0(sample, ".log"))
  
    cat(paste0("Running absolute on sample ", sample, "..."))
    sink(file=log.file)
    RunAbsolute(seg.file, sigma.p, max.sigma.h, min.ploidy, max.ploidy,
                primary.disease, platform, sample, out.dir, max.as.seg.count,
                max.non.clonal, max.neg.genome, copy_num_type, verbose=TRUE,
                maf.fn=maf.file, min.mut.af=min.mut.af)
    sink()
    cat(" done\n")
}

# read mutations from MAF, ref and alt counts from VCF
read.maf.vcf <- function (sample) {
    maf <- read.table(file.path(maf.dir, paste0(sample, ".maf")), header=T, sep="\t", fill=T)
    vcf <- read.table(file.path(vcf.dir, paste0(sample, ".vcf")), header=F, sep="\t", fill=T,
                      col.names=c("CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                                  "FILTER", "INFO", "FORMAT", "NORMAL", "TUMOR"))
    chroms <- c(1:22, "X", "Y")
    nucs <- c("A", "C", "G", "T")
    maf$Chromosome <- factor(maf$Chromosome, levels=chroms)
    vcf$CHROM <- factor(vcf$CHROM, levels=chroms)
    maf$Reference_Allele <- factor(maf$Reference_Allele, levels=nucs)
    vcf$REF <- factor(vcf$REF, levels=nucs)
    
    new.maf <- do.call(rbind, lapply(1:nrow(maf), function (i) {
        row <- maf[i,]
        counts <- subset(vcf, CHROM==row$Chromosome & 
                              POS==row$Start_Position &
                              REF==row$Reference_Allele &
                              ALT==row$Tumor_Seq_Allele1)
        if (nrow(counts) == 0) return (NULL)
        counts <- counts$TUMOR
        counts <- strsplit(as.character(counts), ":", fixed=T)[[1]][5:8]
        counts <- as.integer(sapply(counts, function (x) strsplit(x, ",")[[1]][1]))
        ref.count <- counts[which(nucs==row$Reference_Allele)]
        alt.count <- counts[which(nucs==row$Tumor_Seq_Allele1)]
        cbind(row, t_ref_count=ref.count, t_alt_count=alt.count)
    }))
    colnames(new.maf)[colnames(new.maf) == "Start_Position"] <- "Start_position"
    
    write.table(new.maf, file.path(out.dir, paste0(sample, ".maf")), col.names=T, 
                row.names=F, sep="\t", quote=F)
}

d <- read.table(METADATA, header=T)
DoAbsolute(d$tumor.sample[1])
quit()
#registerDoMC(1)
foreach (sample=d$tumor.sample, .combine=c) %dopar% {
  DoAbsolute(sample)
}

obj.name <- "summary"
absolute.files <- file.path(out.dir, paste0(d$tumor.sample, ".ABSOLUTE.RData"))
CreateReviewObject(obj.name, absolute.files, out.dir, copy_num_type, verbose=TRUE)
## At this point you'd perform your manual review and mark up the file 
## output/abs_summary/DRAWS_summary.PP-calls_tab.txt by prepending a column with
## your desired solution calls. After that (or w/o doing that if you choose to accept
## the defaults, which is what running this code will do) run the following command:
calls.path = file.path(out.dir, "summary.PP-calls_tab.txt")
modes.path = file.path(out.dir, "summary.PP-modes.data.RData")
ExtractReviewedResults(calls.path, "absolute", modes.path, out.dir, obj.name, copy_num_type)
