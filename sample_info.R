#!/usr/bin/env Rscript

source(file="/home/rmccloskey/morin-rotation/settings.conf")

WORK_DIR <- sub("dlbcl", "colorectal", WORK_DIR)
bam.dir <- file.path(WORK_DIR, "00_bams")
flagstat.dir <- file.path(WORK_DIR, "09_flagstat")
coverage.dir <- file.path(WORK_DIR, "02_coverage")

################################################################################
# Functions to read in various files
################################################################################

read.coverage <- function (file.name) {
    cat("Reading coverage for ", file.name, "... ")
    coverage.file <- file.path(coverage.dir, paste0(file.name, ".sample_interval_summary"))
    coverage <- read.table(coverage.file, header=T, fill=T, na.strings="N/A")
    int.range <- strsplit(sub(".*:", "", coverage$Target), "-")
    coverage$start <- as.numeric(sapply(int.range, "[[", 1))
    coverage$end <- as.numeric(sapply(int.range, "[[", 2))
    coverage$length <- coverage$end-coverage$start+1
    coverage$total_coverage <- as.numeric(coverage$total_coverage)
    coverage.mean <- sum(as.numeric(coverage$total_coverage))/sum(as.numeric(coverage$length))
    coverage.sd <- sqrt(sum(coverage$total_coverage^2/coverage$length)/sum(coverage$length) - coverage.mean^2)
    cat("done\n")
    data.frame(coverage.mean=coverage.mean, coverage.sd=coverage.sd,
               file.name=file.name)
}

read.flagstat <- function (file.name) {
    cat("Reading flagstat for ", file.name, "... ")
    flagstat.file <- file.path(flagstat.dir, paste0(file.name, ".txt"))
    ncols <- max(count.fields(flagstat.file))
    classes <- c("numeric", "NULL", "numeric", rep("NULL", ncols-3))
    d <- read.table(flagstat.file, fill=T, colClasses=classes)
    rownames(d) <- c("total", "duplicates", "mapped", "seq.paired", "read1", 
                     "read2", "prop.paired", "self.mate.mapped", "singleton",
                     "mate.mapped.diff.chr", "mate.mapped.diff.chr.mapq.5")
    colnames(d) <- c("qc.pass", "qc.fail")
    prop.fail <- d["total", "qc.fail"]/(d["total", "qc.pass"] + d["total", "qc.fail"])
    prop.dup <- d["duplicates", "qc.pass"]/d["total", "qc.pass"]
    prop.unmapped <- (d["total", "qc.pass"]-d["mapped", "qc.pass"])/d["total", "qc.pass"]
    cat("done\n")
    data.frame(file.name=file.name, prop.fail=prop.fail, prop.dup=prop.dup, 
               prop.unmapped=prop.unmapped)
}

################################################################################
# Main
################################################################################

# read in tumor content estimates
sample.data <- read.csv(EXOME_PURITY, header=T, fill=T, stringsAsFactors=F)
colnames(sample.data)[1] <- "sample"
sample.data <- sample.data[!sample.data$sample %in% c("HT 29", "R3"),]

purity <- (sample.data$X..tumor/100) *
          (sample.data$X..viable.neoplastic.cells/100 + 
           sample.data$X..necrosis/100)
sample.data$purity <- ifelse(is.na(purity), sample.data$X..Tumor/100, purity)
sample.data <- sample.data[,c("sample", "purity")]

# get file names for samples
genome.file.data <- read.csv(GENOME_METADATA, header=T)
exome.file.data <- read.csv(EXOME_METADATA, header=T)
file.data <- merge(genome.file.data, exome.file.data, 
                   by=c("patient_id", "sample_id"),
                   suffixes=c(".genome", ".exome"))

file.data$file.genome <- paste(file.data$patient_id, 
                               file.data$gsc_genome_library_id,
                               "genome.GRCh37-lite.aln", sep="_")
file.data$file.exome <- paste(file.data$patient_id, 
                              file.data$gsc_exome_library_id,
                              "exome.GRCh37-lite.aln", sep="_")
file.data <- file.data[,c("patient_id", "sample_id", "file.exome", "file.genome")]
colnames(file.data)[1:2] <- c("patient.id", "sample")
sample.data <- merge(sample.data, file.data, all.y=T)

# read coverage info
ex.coverage <- do.call(rbind, lapply(sample.data$file.exome, read.coverage))
sample.data <- merge(sample.data, ex.coverage, 
                     by.x=c("file.exome"), by.y=c("file.name"))
gen.coverage <- do.call(rbind, lapply(sample.data$file.genome, read.coverage))
sample.data <- merge(sample.data, gen.coverage, 
                     by.x=c("file.genome"), by.y=c("file.name"),
                     suffixes=c(".exome", ".genome"))

# read flagstat info
ex.flagstat <- do.call(rbind, lapply(sample.data$file.exome, read.flagstat))
sample.data <- merge(sample.data, ex.flagstat, 
                     by.x=c("file.exome"), by.y=c("file.name"))
gen.flagstat <- do.call(rbind, lapply(sample.data$file.genome, read.flagstat))
sample.data <- merge(sample.data, gen.flagstat, 
                     by.x=c("file.genome"), by.y=c("file.name"),
                     suffixes=c(".exome", ".genome"))

write.table(sample.data, file="sample_info.dat", col.names=T, row.names=F)
