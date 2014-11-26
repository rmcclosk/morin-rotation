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
        params["position"] = row["Start_position"]
        params["ensembl_gene_id"] = row["Entrez_Gene_Id"]
        params["base_change"] = "{}>{}".format(ref_base, nref_base)
        params["protein_altering"] = row["Protein_Change"] != ""
        params["cdna_change"] = row["cDNA_Change"]
        params["identifiers"] = row["dbSNP_RS"].replace("&", ",")
        params["splice_site"] = row["Variant_Classification"] == "Splice_Site"
        params["tumour_ref"] = row["t_ref_count"]
        params["tumour_nref"] = row["t_alt_count"]
        params["normal_ref"] = 0
        params["normal_nref"] = 0
    #def addMutation(self, library_id, chromosome, position, ensembl_gene_id, base_change, status=None, protein_altering=None, annotation=None, validation_outcome=None, cdna_change = None, to_validate=None,identifiers=None, splice_site=None, mutation_seq_probability=None, triplet=None, tumour_ref = None, tumour_nref = None, normal_ref = None, normal_nref = None,transcript=None,sift_score=None,polyphen_score=None,mutation_ass_score=None,ref_base=None,nref_base=None):
        db.addMutation(library_id, chromosome, position, ensembl_gene_id,
                       base_change, protein_altering=protein_altering,
                       annotation=annotation, cdna_change=cdna_change,
                       identifiers=identifiers, splice_site=splice_site,
                       mutation_seq_probability=0, tumour_ref=tumour_ref,
                       tumour_nref=tumour_nref, transcript=transcript,
                       ref_base=ref_base, nref_base=nref_base)
    # indel
    else:
        params["start"] = row["Start_position"]
        params["end"] = row["End_Position"]
        params["ensembl_id"] = row["Entrez_Gene_Id"]
        addIndel(library_id, chromosome, start, end, ref_base, alt_base, effect,
                 annotation, ensembl_id=ensembl_gene_id)
