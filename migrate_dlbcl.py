#!/usr/bin/env python

# Migrate exome data from apollo to genesis.

import csv
import glob
import sys
import imp

settings = imp.load_source("settings", "settings.conf")

from_dir = "/projects/rmorin/analysis/aligned/DLBCL_EXOME/all_bams"
to_dir = "/genesis/extscratch/morinlab/shared/rmccloskey/dlbcl/00_bams"

for from_file in glob.glob(from_dir + "/*.ba[im]"):
    to_file = "%s/%s" % (to_dir, from_file)
    print(to_file)
    break
    to_file = "/genesis/extscratch/morinlab/shared/qcroc_colorectal/exome/%s" % to_file
    try:
        from_file = glob.glob("%s/*.bam" % row["gsc_alignment_path"])[0]
    except:
        print(row["gsc_alignment_path"])
        sys.exit()
    command = "cp %s %s" % (from_file, to_file)
    print(command)
