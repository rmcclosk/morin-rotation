#!/usr/bin/env Rscript

source(file="/home/rmccloskey/morin-rotation/settings.conf")

out.dir <- file.path(WORK_DIR, "12_segcoverage")
seg.dir <- file.path(WORK_DIR, "05_hmmcopy", "10000")
bam.dir <- file.path(WORK_DIR, "01_fixbams")
dir.create(out.dir, showWarnings=F)

# get list of samples
d <- read.csv(GENOME_METADATA, header=T)
d <- d[,c("patient_id", "sample_id", "gsc_genome_library_id")]
d <- d[!grepl("(BF|WB)", d$sample_id),]
d$bam.file <- file.path(bam.dir, paste(d$patient_id, d$gsc_genome_library_id, 
                        "genome.GRCh37-lite.aln.bam", sep="_"))
d$list.file <- paste0(file.path(out.dir, d$sample_id), ".list")
d$seg.file <- file.path(seg.dir, d$patient_id, 
                       paste0(d$sample_id, "_kmeans_segments.dat"))
d$outfile.stem <- file.path(out.dir, d$sample_id)
d$coverage.file <- paste0(d$outfile.stem, ".sample_interval_summary")
d <- d[file.exists(d$seg.file) & !file.exists(d$coverage.file),]

# read lengths of each contig
genome.dict <- read.table(sub("fa", "dict", HUMAN_REF), fill=T)
genome.dict <- genome.dict[genome.dict$V1 == "@SQ",]
genome.dict$chr <- sub("SN:", "", genome.dict$V2)
genome.dict$length <- as.integer(sub("LN:", "", genome.dict$V3))
genome.dict <- genome.dict[,c("chr", "length")]

# make segments lists
file.remove("jobs.txt")
options("scipen"=100, "digits"=4)
apply(d, 1, function (row) {
    seg <- read.table(row["seg.file"], header=T)
    seg <- merge(seg, genome.dict)
    seg$end <- ifelse(seg$end > seg$length, seg$length, seg$end)
    seg$chr <- factor(seg$chr, c(1:22, "X", "Y"))
    seg <- seg[order(seg$chr, seg$start),]
    
    seg.list <- data.frame(seg=paste0(seg$chr, ":", seg$start, "-", seg$end))
    write.table(seg.list, file=row["list.file"],
                col.names=F, row.names=F, quote=F)

    cat(JAVA_BIN, "-Xmx2048m -jar", GATK_JAR, "-T DepthOfCoverage",
        "-omitBaseOutput -omitLocusTable -R", HUMAN_REF, "-I", row["bam.file"],
        "-L", row["list.file"], "-o", row["outfile.stem"], "-K", GATK_KEY, "\n",
        file="jobs.txt", append=T)
})

if (file.exists("jobs.txt")) {
    system('mqsub --file jobs.txt --name segcoverage --chdir qsub-logs --qsub "-l h_vmem=3G -l mem_free=3G -l mem_token=3G"')
}
