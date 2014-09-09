#!/usr/bin/env python

# Migrate exome data from apollo to genesis.

import csv
import glob
import sys
import imp

settings = imp.load_source("settings", "settings.conf")

libs = ["A45012", "A45020", "A45035", "A45005", "A45028", "A44997"]

reader = csv.DictReader(open("exome_metadata.csv"))
for row in reader:
    if not row["gsc_exome_library_id"] in libs:
        continue
    to_file = "%s_%s_exome.GRCh37-lite.aln.bam.orig" % (row["patient_id"], row["gsc_exome_library_id"])
    to_file = "/genesis/extscratch/morinlab/shared/qcroc_colorectal/exome/%s" % to_file
    try:
        from_file = glob.glob("%s/*.bam" % row["gsc_alignment_path"])[0]
    except:
        print(row["gsc_alignment_path"])
        sys.exit()
    command = "cp %s %s" % (from_file, to_file)
    print(command)
