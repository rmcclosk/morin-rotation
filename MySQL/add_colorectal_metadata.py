#!/usr/bin/env python

# Add colorectal metadata to the database.
# This script is full of BADNESS (queries with string replacement) because I am
# lazy.


import cancerGenome
import csv
import sys

db = cancerGenome.cancerGenomeDB(
    database_name='colorectal',
    database_host='jango.bcgsc.ca',
    database_user='rmccloskey',
    database_password='rmccloskey'
)
cursor = db.db.cursor()

f = "/extscratch/morinlab/shared/qcroc_colorectal/exome/sample_metadata.txt"
reader = csv.DictReader(open(f), delimiter="\t")
for row in reader:
    if row["patient_id"] == "R3":
        row["patient_id"] = "HT29"

    db.addPatient(row["patient_id"])
    query = "UPDATE patient SET sex='{}' WHERE res_id='{}'".format(row["sex"][0], row["patient_id"])
    cursor.execute(query)

    sample_id = db.addSample(format(row["sample_id"]), "'{}'".format(row["patient_id"]), "primary")
    db.addLibrary(format(row["gsc_exome_library_id"]), "exome", sample_id)

f = "/extscratch/morinlab/shared/qcroc_colorectal/genome/sample_metadata.txt"
reader = csv.DictReader(open(f), delimiter="\t")
for row in reader:
    if row["patient_id"] == "R3":
        row["patient_id"] = "HT29"

    sample_id = db.addSample(format(row["sample_id"]), "'{}'".format(row["patient_id"]), "primary")
    db.addLibrary(format(row["gsc_genome_library_id"]), "genome", sample_id)
