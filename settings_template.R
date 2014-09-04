# File listing matched tumor-normal pairs.  The file should contain six
# columns: normal ID, normal file, pre ID, pre file, post ID, post file.
matching.file <- "sample_matching.txt"

# Location of BAM files and indices.
bam.dir <- "."

# Where to put coverage files output by GATK.
coverage.dir <- "coverage"

# Where to put downsampled data.
sample.dir <- "downsampled"

# Where to put coverage files for downsampled data.
sample.coverage.dir <- "coverage_downsampled"

# GATK jar file.
gatk.jar <- "GenomeAnalysisTK.jar"

# Where is java?
java.bin <- "/usr/bin/java"

# Human genome reference.
human.fasta <- "human_g1k_v37.fasta"

# List of exome regions.
exome.region.list <- "SureSelect_regions.list"
