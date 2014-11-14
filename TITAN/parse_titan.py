#!/usr/bin/env python3

# Get the TITAN results for the best cluster, copy to an output directory, and
# compress.

import os.path
import glob
import argparse
import re
import subprocess
import shutil

parser = argparse.ArgumentParser("Parse TITAN output")
folder_help = "Folder containing TITAN results"
parser.add_argument("folder", metavar="F", help=folder_help)
stem_help = "File name stem for results"
parser.add_argument("output_stem", metavar="O", help=stem_help)
args = parser.parse_args()

ptn = os.path.join(args.folder, "titan", "*params.txt")
param_files = glob.glob(ptn)

cluster_ptn = re.compile(".*cluster_(?P<cluster>\d+)_.*")

best_cluster = 0
best_sdbw = 10E10
for param_file in param_files:
    cluster = int(cluster_ptn.search(param_file).group("cluster"))
    with open(param_file) as f:
        for line in f:
            if line.startswith("S_Dbw validity index"):
                sdbw = float(line.strip().split()[-1])
                break
    print(sdbw)
    if sdbw < best_sdbw:
        best_sdbw = sdbw
        best_cluster = cluster

print(best_cluster)
file_ptn = re.compile(".*cluster_{}(?P<suffix>[_.].+)$".format(best_cluster))
for f in os.listdir(os.path.join(args.folder, "titan")):
    print(f)
    match = file_ptn.match(f)
    if match:
        new_name = "{}{}".format(args.output_stem, match.group("suffix"))
        shutil.copyfile(os.path.join(args.folder, "titan", f), new_name)
        subprocess.call(["bzip2", new_name])
