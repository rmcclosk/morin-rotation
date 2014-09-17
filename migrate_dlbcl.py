#!/usr/bin/env python

# Migrate exome data from apollo to genesis.

import os
import glob

from_dir = "/projects/rmorin/analysis/aligned/DLBCL_EXOME/all_bams"
to_dir = "/genesis/extscratch/morinlab/shared/rmccloskey/dlbcl/00_bams"

for from_file in glob.glob(from_dir + "/*.ba[im]"):
    to_file = "%s/%s" % (to_dir, os.path.basename(from_file))
    command = "cp %s %s" % (from_file, to_file)
    print(command)
