#!/usr/bin/env Rscript

#library(ABSOLUTE)
library(numDeriv)

sink("/dev/null")
sapply(Sys.glob("ABSOLUTE/R/*"), function (f) source(file=f))
sink()
print("starting")
RunAbsolute("mix250K_seg_out.txt",
            min.ploidy=0.95, 
            max.ploidy=10, 
            max.sigma.h=0.02, 
            sigma.p=0, 
            platform="Illumina_WES", 
            copy_num_type="total",
            results.dir="test", 
            primary.disease="cancer", 
            sample.name="foo", 
            max.as.seg.count=1500)
print("done")
load("test/foo.ABSOLUTE.RData")
seg.dat[["mode.res"]][["mode.flag"]]

#CreateReviewObject("test", file.path("test", "foo.ABSOLUTE.RData"), "test-out", "total", verbose=T)
