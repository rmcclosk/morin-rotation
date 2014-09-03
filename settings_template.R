# File listing matched tumor-normal pairs.  The file should contain four
# columns: Normal, Pre, ?, Post. Normal, Pre, and Post should have a file name.
# I don't know what the third column is.
matching.file <- "sample_matching.txt"

# Location of BAM files and indices.
bam.dir <- "."

# GATK jar file.
gatk.jar <- "/home/rmccloskey/bin/Genome..."

# Human genome reference.
human.fasta <- "human_g1k_v37.fasta"

# List of exome regions.
exome.region.list <- "SureSelect_regions.list"
