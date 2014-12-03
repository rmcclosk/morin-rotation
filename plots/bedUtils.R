# A couple of useful functions for working with BED files.

# Write a data frame to a temporary file, in BED format.
# The data frame must already have the correct columns (ie. the first three
# columns must be chromosome, start, end).
write.bed <- function (bed) {
    f <- tempfile()
    write.table(bed, f, col.names=F, row.names=F, quote=F, sep="\t")
    f
}

# Call a bedtools command on a pair of BED data frames.
bedtools <- function (cmd, b1, b2, args=NULL) {
    f1 <- write.bed(b1)
    f2 <- write.bed(b2)
    cmd <- paste("bedtools", cmd, args, "-a", f1, "-b", f2)
    output <- system(cmd, intern=T)
    if (length(output) == 0) return (NULL)
    res <- setNames(read.table(textConnection(output)), colnames(b1))
    res$chr <- factor(res$chr, levels=c(1:22, "X"))
    res
}
