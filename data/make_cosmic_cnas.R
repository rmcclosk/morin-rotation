#!/usr/bin/env Rscript

temp.bedfile <- function (df, dir=tempdir()) {
    if (!file.exists(dir))
        dir.create(dir)
    f <- tempfile(tmpdir=dir)
    write.table(df, file=f, col.names=F, row.names=F, quote=F, sep="\t")
    f
}

do.intersect <- function (bed.file.1, bed.file.2) {
    tmpdir <- tempdir()
    dir.create(tmpdir)
    if (class(bed.file.1) == "data.frame")
        bed.file.1 <- temp.bedfile(bed.file.1, tmpdir)
    if (class(bed.file.2) == "data.frame")
        bed.file.2 <- temp.bedfile(bed.file.2, tmpdir)

    cmd <- paste("bedtools intersect -a", bed.file.1, "-b", bed.file.2)
    cat(cmd, "\n")
    res <- read.table(textConnection(system(cmd, intern=T)))
    colnames(res)[1:3] <- c("chr", "start", "end")
    unlink(tmpdir, recursive=T)
    res
}

data.file <- "CosmicCompleteCNA.tsv"
chrs <- c(1:22, "X")

if (!file.exists(data.file)) {
    url <- "http://cancer.sanger.ac.uk/cancergenome/projects/cosmic/download"
    cat("Please download", data.file, "from", url, "\n")
    quit(status=1)
}

# select only colorectal cancer, and relevant columns
cmd <- paste("grep large_intestine", data.file, "| cut -f 11,13")
con <- textConnection(system(cmd, intern=T))
cna.data <- read.table(con, sep="\t", col.names=c("mutation", "locus"))

# extract chromosome and coordinates
ptn <- "([0-9]+):([0-9]+)[.][.]([0-9]+)"
cna.data$chr <- factor(gsub("23", "X", gsub(ptn, "\\1", cna.data$locus)))
cna.data$start <- as.numeric(gsub(ptn, "\\2", cna.data$locus))
cna.data$end <- as.numeric(gsub(ptn, "\\3", cna.data$locus))
cna.data <- cna.data[,c("chr", "start", "end", "mutation")]

# aggregate losses and gains
losses <- do.intersect("genes.bed", subset(cna.data, mutation == "LOSS"))
colnames(losses) <- c("chr", "start", "end", "gene")
losses <- setNames(aggregate(chr~gene, losses, length), c("gene", "losses"))

gains <- do.intersect("genes.bed", subset(cna.data, mutation == "GAIN"))
colnames(gains) <- c("chr", "start", "end", "gene")
gains <- setNames(aggregate(chr~gene, gains, length), c("gene", "gains"))

# combine losses and gains, and write to file
cnas <- merge(gains, losses, all=T)
cnas[is.na(cnas$gains), "gains"] <- 0
cnas[is.na(cnas$losses), "losses"] <- 0
write.table(cnas, "cosmic-cnas.tsv", row.names=F, col.names=T, sep="\t", quote=F)
