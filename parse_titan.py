#!/usr/bin/env python3

import imp
import csv
import os
import glob

settings = imp.load_source("settings", "/home/rmccloskey/morin-rotation/settings.conf")
titan_dir = os.path.join(settings.WORK_DIR, "14_titan")

for line in csv.DictReader(open(settings.METADATA), delimiter="\t"):
    results_dir = os.path.join(titan_dir, line["tumor.sample"], "titan")
    param_files = glob.glob(os.path.join(results_dir, "*params.txt"))
    if len(param_files) != 5:
        continue
    normal_est = 0
    for f in param_files:
        normal_est += float(next(open(f)).strip().split()[-1])
    normal_est /= 5
    tumor_est = round((1-normal_est)*100)
    print(line["tumor.sample"], tumor_est, line["purity"])
