#!/usr/bin/env python3

import imp
import csv
import os
import re
import glob
import subprocess
import shutil

settings = imp.load_source("settings", "/home/rmccloskey/morin-rotation/settings.conf")
titan_dir = os.path.join(settings.WORK_DIR, "14_titan")
out_dir = "tmp"

for row in csv.DictReader(open(settings.METADATA), delimiter="\t"):
    results_dir = os.path.join(titan_dir, row["tumor.sample"], "titan")
    param_files = glob.glob(os.path.join(results_dir, "*params.txt"))
    best_sdbw = 10E10
    best_cluster = 0
    for i, f in enumerate(param_files, start=1):
        proc = subprocess.Popen(["tail", "-n", "1", f], stdout=subprocess.PIPE)
        line = proc.communicate()[0]
        sdbw = float(line.strip().split()[-1])
        if (sdbw < best_sdbw):
            best_sdbw = sdbw
            best_cluster = i

    if (best_cluster == 0): continue
    param_ptn = "*cluster_{}_params.txt".format(best_cluster)
    param_file = glob.glob(os.path.join(results_dir, param_ptn))[0]
    admix = float(next(open(param_file)).strip().split()[-1])
    purity = 1-admix
    print(row["tumor.sample"], purity)
    continue

    seg_ptn = "*cluster_{}*.seg".format(best_cluster)
    seg_file = glob.glob(os.path.join(results_dir, seg_ptn))[0]

    out_file = "{}.seg".format(os.path.join(out_dir, row["tumor.sample"]))
    sf = open(seg_file)
    next(sf)
    with open(out_file, "w") as f:
        f.write("Sample\tChromosome\tStart\tEnd\tNum_Probes\tSegment_Mean\n")
        for line in sf:
            f.write(line)
    sf.close()

    continue
    sample_dir = os.path.join(out_dir, row["tumor.sample"])
    if not os.path.exists(sample_dir):
        os.mkdir(sample_dir)

    txt_ptn = "*cluster_{}*.txt".format(best_cluster)
    image_ptn = "*cluster_{}/*.png".format(best_cluster)
    files = glob.glob(os.path.join(results_dir, txt_ptn))
    files += glob.glob(os.path.join(results_dir, image_ptn))
    for f in files:
        new_name = f.split("_")[-1]
        shutil.copyfile(f, os.path.join(sample_dir, new_name))
