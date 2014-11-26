#!/usr/bin/env python

# Import MAF files into the database.

import argparse
import cancerGenome
import sys
import csv

parser = argparse.ArgumentParser()
parser.add_argument("maf", metavar="M", help="MAF file to add to database")
parser.add_argument("library_id", metavar="L", type=int, help="library ID")
parser.add_argument("db_host", metavar="P", help="database host")
parser.add_argument("db_name", metavar="D", help="database name")
parser.add_argument("db_user", metavar="U", help="database user")
parser.add_argument("db_pass", metavar="P", help="database password")
args = parser.parse_args()

db = cancerGenome.cancerGenomeDB(
    database_name=args.db_name,
    database_host=args.db_host,
    database_user=args.db_user,
    database_password=args.db_pass
)

cursor = db.db.cursor()
reader = csv.DictReader(open(args.maf), delimiter="\t")

for row in reader:
    params = {"library_id": library_id, "type": "snv", "region": "CDS"}
    params["chromosome"] = row["Chromosome"]
    params["annotation"] = row["Protein_Change"]
    transcript = row["Annotation_Transcript"]
    ref_base = row["Reference_Allele"]
    nref_base = row["Tumor_Seq_Allele1"]
    effect = row["Variant_Classification"].lower()

    # SNV
    if len(ref_base) == len(nref_base):
        params["position"] = row["Start_Position"]
        params["ensembl_gene_id"] = row["Entrez_Gene_Id"]
        params["base_change"] = "{}>{}".format(ref_base, nref_base)
        params["protein_altering"] = row["Protein_Change"] != ""
        params["cdna_change"] = row["cDNA_Change"]
        params["identifiers"] = row["dbSNP_RS"].replace("&", ",")
        params["splice_site"] = row["Variant_Classification"] == "Splice_Site"
        params["tumour_ref"] = row["t_ref_count"]
        params["tumour_nref"] = row["t_alt_count"]
        params["normal_ref"] = row["n_ref_count"]
        params["normal_nref"] = row["n_alt_count"}
        params["mutation_seq_probability"] = 0
        params["sift_score"] = 0
        params["polyphen_score"] = 0
        params["mutation_ass_score"] = 0
        db.addMutation(**params)

    # indel
    else:
        params["start"] = row["Start_position"]
        params["end"] = row["End_Position"]
        params["ensembl_id"] = row["Entrez_Gene_Id"]
        addIndel(library_id, chromosome, start, end, ref_base, alt_base, effect,
                 annotation, ensembl_id=ensembl_gene_id)
