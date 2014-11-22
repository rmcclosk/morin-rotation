#!/usr/bin/env python

# Import MAF files into the database.
#NOTE: This is the MAF format
#Hugo_Symbol Entrez_Gene_Id Center NCBI_Build Chromosome Start_position End_Position Strand Variant_Classification Variant_Type Reference_Allele Tumor_Seq_Allele1 Tumor_Seq_Allele2 dbSNP_RS dbSNP_Val_Status Tumor_Sample_Barcode Matched_Norm_Sample_Barcode Match_Norm_Seq_Allele1 Match_Norm_Seq_Allele2 Tumor_Validation_Allele1 Tumor_Validation_Allele2 Match_Norm_Validation_Allele1 Match_Norm_Validation_Allele2 Verification_Status Validation_Status Mutation_Status Sequencing_Phase Sequence_Source Validation_Method Score BAM_File Sequencer Tumor_Sample_UUID Matched_Norm_Sample_UUID Annotation_Transcript Transcript_Position cDNA_Change Protein_Change effect categ t_ref_count t_alt_count

import argparse
import cancerGenome
import sys
import csv
from db_utils import *

db = cancerGenome.cancerGenomeDB(
    database_name="colorectal",
    database_host="jango.bcgsc.ca",
    database_user="rmccloskey",
    database_password="rmccloskey"
)

parser = argparse.ArgumentParser()
parser.add_argument("maf", metavar="M", help="MAF file to add to database")
parser.add_argument("sample", metavar="S", help="sample ID")
parser.add_argument("type", metavar="L", help="library type", choices=LIBRARY_TYPES)
args = parser.parse_args()

cursor = db.db.cursor()
library_id = get_library(cursor, args.sample, args.type)
reader = csv.DictReader(open(args.maf), delimiter="\t")

query = """INSERT INTO mutation (library_id, type, region, chromosome,
                                 position, gene, protein_altering, base_change,
                                 annotation, external_gene_id, transcripts,
                                 reference_base_count, nonreference_base_count,
                                 cdna_change, known_identifiers, ref_base, 
                                 nref_base)
           VALUES (%(library_id)s, %(type)s, %(region)s, %(chromosome)s,
                   %(position)s, %(gene)s, %(protein_altering)s, 
                   %(base_change)s, %(annotation)s, %(external_gene_id)s,
                   %(transcripts)s, %(reference_base_count)s, 
                   %(nonreference_base_count)s, %(cdna_change)s,
                   %(known_identifiers)s, %(ref_base)s, %(nref_base)s)"""

for row in reader:
    params = {"library_id": library_id, "type": "snv", "region": "CDS"}
    params["chromosome"] = row["Chromosome"]
    params["position"] = row["Start_position"]
    params["gene"] = row["Entrez_Gene_Id"]
    params["protein_altering"] = "no" if row["Protein_Change"] == "" else "yes"
    params["base_change"] = "{}>{}".format(row["Reference_Allele"], row["Tumor_Seq_Allele1"])
    params["annotation"] = row["Protein_Change"]
    params["external_gene_id"] = row["Hugo_Symbol"]
    params["transcripts"] = row["Annotation_Transcript"]
    params["reference_base_count"] = row["t_ref_count"]
    params["nonreference_base_count"] = row["t_alt_count"]
    params["cdna_change"] = row["cDNA_Change"]
    params["known_identifiers"] = row["dbSNP_RS"].replace("&", ",")
    params["ref_base"] = row["Reference_Allele"]
    params["nref_base"] = row["Tumor_Seq_Allele1"]
    cursor.execute(query, params)
