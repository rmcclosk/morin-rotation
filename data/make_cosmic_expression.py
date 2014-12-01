#!/usr/bin/env python3

import csv
import os
import sys
import itertools

cna_file = "CosmicCompleteCNA.tsv"
expr_file = "CosmicCompleteGeneExpression.tsv"

if not (os.path.exists(cna_file) and os.path.exists(expr_file)):
    url = "cancer.sanger.ac.uk/cancergenome/projects/cosmic/download"
    msg = "Please download {} and {} from {}"
    sys.exit(msg.format(cna_file, expr_file))

reader = csv.DictReader(open(cna_file), delimiter="\t")
filter_fun = lambda r: r["Primary site"] == "large_intestine"
group_fun = lambda r: r["ID_SAMPLE"]

filtered_reader = filter(filter_fun, reader)
grouped_reader = itertools.groupby(reader, key=group_fun)
sample_ids = set([k for k, _ in grouped_reader])

reader = csv.DictReader(open(expr_file), delimiter="\t")
reader = filter(lambda r: r["ID_SAMPLE"] in sample_ids, reader)
gene_counts = {}
for row in reader:
    regulation = row["REGULATION"] 
    if regulation != "normal" and row["ID_SAMPLE"] in sample_ids:
        gene = row["GENE_NAME"]
        idx = 0 if regulation == "under" else 1
        try:
            gene_counts[gene][idx] += 1
        except KeyError:
            gene_counts[gene] = [0, 0]
            gene_counts[gene][idx] += 1

writer = csv.writer(open("cosmic-expression.tsv", "w"), delimiter="\t")
writer.writerow(["gene", "under", "over"])
for gene, counts in sorted(gene_counts).items():
    writer.writerow([gene, counts[0], counts[1]])
