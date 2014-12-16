#!/usr/bin/env python

import subprocess
import argparse
import tempfile
import shutil
import multiprocessing
import os
import sys

def _do_view(args):
    bam, regions, tmpdir = args
    fd, fn = tempfile.mkstemp(dir=tmpdir, suffix=".bam")
    view_command = ["samtools", "view", "-u", "-b", "-F", "0x0600", bam]
    view_command += regions
    subprocess.Popen(view_command, stdout=fd).communicate()
    os.close(fd)
    return fn

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract exome reads from a BAM")
    parser.add_argument("bam", help="input BAM file")
    parser.add_argument("intervals", help="list of intervals")
    parser.add_argument("outfile", help="output BAM file")
    parser.add_argument("nthreads", type=int, help="number of threads")
    args = parser.parse_args()
    
    nper = 1000
    read_length = 100

    regions = []
    buf = []
    prev_end = 0
    prev_chr = None
    dump = False
    for i, line in enumerate(open(args.intervals)):
        chr, bounds = line.strip().split(":")
        start, end = [int(b) for b in bounds.split("-")]

        if i % nper == 0 and i > 0:
            dump = True

        if dump:
            if chr != prev_chr or start - prev_end > read_length:
                regions.append(buf)
                buf = []
                dump = False

        buf.append(line.strip())
        prev_end = end
        prev_chr = chr

    if len(buf) > 0:
        regions.append(buf)
    
    try:
        tmpdir = tempfile.mkdtemp()
        pool = multiprocessing.Pool(args.nthreads)
        pool_args = []
        for r in regions:
            pool_args.append((args.bam, r, tmpdir))

        files = [x for x in pool.map(_do_view, pool_args)]

        print("Concatenating")
        fd, fn = tempfile.mkstemp(suffix=".bam")
        cat_command = ["samtools", "cat", "-o", fn] + files
        subprocess.Popen(cat_command).communicate()
        os.close(fd)

        print("Sorting")
        sort_command = ["samtools", "sort", "-T", fn, "-O", "bam", fn]
        fh = open(args.outfile, "w")
        subprocess.Popen(sort_command, stdout=fh).communicate()
        fh.close()

    finally:
        shutil.rmtree(tmpdir)
