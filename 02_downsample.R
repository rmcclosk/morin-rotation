#!/home/rmccloskey/bin/Rscript

# Downsample sequences from tumor or normal to make the same average coverage
# per sample.

source(file="settings.R")
matching <- read.table(matching.file, sep="\t", fill=T, na.strings=c("", "Insufficient DNA"))
colnames(matching) <- c("normal", "diagnosis", "relapse")

# Keep only rows where at least one comparison can be made.
matching <- matching[rowSums(is.na(matching)) <= 1,]

# Extract file names.
matching <- data.frame(apply(matching, 2, function (col) {
	ifelse(is.na(col), NA, sapply(strsplit(col[!is.na(col)], " "), "[[", 2))
}))

# Read mean coverage for each sample.
coverage.files <- Sys.glob(file.path(coverage.dir, "*.sample_summary"))
bam.files <- sub(".sample_summary$", ".bam", basename(coverage.files))
mean.coverage <- sapply(coverage.files, function (f) {
	bam.file <- sub(".sample_summary$", ".bam", basename(f))
	coverage <- read.table(f, header=T, fill=T, sep="\t", na.strings=c("", "N/A"))
	coverage[coverage$sample_id == "Total", "mean"]
})
coverage <- data.frame(coverage=mean.coverage, row.names=bam.files)

# Find factor to down-sample by, to equate means.
scale.factors <- do.call(rbind, apply(matching, 1, function (row) {
	row <- row[!is.na(row)]
	data.frame(scale.factor=min(coverage[row,], na.rm=T)/coverage[row,],
		   bam.file=row)
}))

# Do the down-sampling.
dir.create(sample.dir, showWarnings=F)
main <- apply(scale.factors, 1, function (row) {
	scale.factor <- as.numeric(row["scale.factor"])
	in.bam <- file.path(bam.dir, row["bam.file"])
	out.bam <- file.path(sample.dir, row["bam.file"])
	if (!file.exists (out.bam)) {
		if (scale.factor < 1) {
			print('system(sprintf("samtools view -s %f -b %s > %s", scale.factor, in.bam, out.bam))')
		} else {
			coverage.file <- sub(".bam$", ".sample_summary", row["bam.file"])
			in.coverage <- file.path(coverage.dir, coverage.file)
			out.coverage <- file.path(sample.coverage.dir, coverage.file)
			file.symlink(in.bam, out.bam)
			file.symlink(in.coverage, out.coverage)
		}
	} else {
		cat(sprintf("Output file %s exists\n", out.bam))
	}
})
