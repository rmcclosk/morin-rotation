#!/usr/bin/env Rscript

source(file="../settings.conf")

d <- read.table(METADATA, header=T, stringsAsFactors=F)

sr.dir <- file.path(WORK_DIR, "07_yaha")
pe.dir <- file.path(WORK_DIR, "08_discordant")
hist.dir <- file.path(WORK_DIR, "10_pairedend")
out.dir <- file.path(WORK_DIR, "11_lumpy")

read.length <- 100
min.non.overlap <- 100 
discordant.z <- 4
back.distance <- 20
weight <- 1
min.mapping.threshold <- 20
exclude.file <- "/home/rmccloskey/morin-rotation/genome/btu356_LCR-hs37d5.bed"
min.weight <- 4
trim.threshold <- 0.0

lumpy.template <- "lumpy -mw %d -tt %f %s > %s"
pe.template <- paste0("-pe ",
                      "bam_file:%s,",
                      "histo_file:%s,",
                      "mean:%d,",
                      "stdev:%d,",
                      "read_length:%d,",
                      "min_non_overlap:%d,",
                      "discordant_z:%d,",
                      "back_distance:%d,",
                      "weight:%d,",
                      "id:%s,",
                      "min_mapping_threshold:%d")
sr.template <- paste0("-sr ",
                      "bam_file:%s,",
                      "back_distance:%d,",
                      "weight:%d,",
                      "id:%s,",
                      "min_mapping_threshold:%d")

cmds <- by(d, d$patient, function (pd) {
    samples <- unique(c(pd$normal.sample, pd$tumor.sample))
    samples <- samples[!is.na(samples)]
    args <- paste(sapply(1:length(samples), function (i) {
        s <- samples[i]
        sr.file <- file.path(sr.dir, paste0(s, ".bam"))
        pe.file <- file.path(pe.dir, paste0(s, ".bam"))
        hist.file <- file.path(hist.dir, paste0(s, ".histo"))
        stat.file <- file.path(hist.dir, paste0(s, ".dat"))

        if (!file.exists(stat.file)) return ("")

        stat <- read.table(stat.file, stringsAsFactors=F)
        mean <- round(as.numeric(strsplit(stat[1,1], ":")[[1]][2]))
        stdev <- round(as.numeric(strsplit(stat[1,2], ":")[[1]][2]))

        pe.id <- 2*i - 1
        sr.id <- 2*i
        cat(s, pe.id, sr.id, "\n")

        pe.cmd <- sprintf(pe.template, pe.file, hist.file, mean, stdev,
                          read.length, min.non.overlap, discordant.z,
                          back.distance, weight, pe.id, min.mapping.threshold)
        sr.cmd <- sprintf(sr.template, sr.file, back.distance, weight, sr.id,
                          min.mapping.threshold)
        paste(pe.cmd, sr.cmd)
    }), collapse=" ")
    out.file <- file.path(out.dir, paste0(pd[1,"patient"], ".bedpe"))
    if (!file.exists(out.file))
        sprintf(lumpy.template, min.weight, trim.threshold, args, out.file)
    else ""
})

#cmds <- paste("source /home/rmccloskey/.bash_profile;", cmds)
cmds <- paste("", cmds)
write.table(cmds, file="jobs.txt", quote=F, row.names=F, col.names=F)
