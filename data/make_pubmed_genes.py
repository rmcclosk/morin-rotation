#!/usr/bin/env python3

import csv
import time
from pprint import pprint
from Bio import Entrez
Entrez.email = "rmccloskey@alumni.ubc.ca"

reader = csv.reader(open("genes.bed"), delimiter="\t")
out_handle = open("pubmed.tsv", "a")
writer = csv.writer(out_handle, delimiter="\t")

start = False
for i, row in enumerate(reader):
    gene = row[-1]
    if gene == "FAM201A":
        start = True

    if not start: continue

    print(gene)

    term = "colorectal cancer {}".format(gene)
    handle = Entrez.esearch(db="pubmed", term=term)
    record = Entrez.read(handle)
    handle.close()

    try:
        record["ErrorList"]
        writer.writerow([gene, 0])
    except KeyError:
        writer.writerow([gene, record["Count"]])
    
    if i % 10 == 0:
        out_handle.flush()
    time.sleep(1)
