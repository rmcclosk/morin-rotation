#!/usr/bin/env python

# Import a SEG file containing CNVs into the database.
# The SEG file must have the format specified here: 
# http://www.broadinstitute.org/igv/SEG
# In particular, the fifth column should be the copy number.

import argparse
import cancerGenome
import csv
from db_utils import *

def addLOH(cursor, library_id, chr, start, end, copy_number, loh_state):
    """Add a LOH event to the database"""
    event_id = addEvent(cursor, library_id, "loh")
    query = """INSERT INTO loh (event_id, chromosome, segment_start, 
                                segment_end, size, copy_number, loh_state)
               VALUES (%s, %s, %s, %s, %s, %s, %s)"""
    size = int(end)-int(start)
    params = (event_id, chr, start, end, size, copy_number, loh_state)
    try:
        cursor.execute(query, params)
        return cursor.lastrowid
    except:
        print("Couldn't add LOH")
        query = "DELETE FROM event WHERE id = %s"
        cursor.execute(query, event_id)
        return None

def addCNV(cursor, library_id, chr, start, end, copy_number, cnv_type):
    """Add a copy number variation to the database"""
    event_id = addEvent(cursor, library_id, "CNV")
    query = """INSERT INTO cnv (event_id, chromosome, segment_start, 
                                segment_end, segment_state, size, type)
               VALUES (%s, %s, %s, %s, %s, %s, %s)"""
    size = int(end)-int(start)
    params = (event_id, chr, start, end, copy_number, size, cnv_type)
    try:
        cursor.execute(query, params)
        return cursor.lastrowid
    except:
        print("Couldn't add CNV")
        query = "DELETE FROM event WHERE id = %s"
        cursor.execute(query, event_id)
        return None

if __name__ == "__main__":
    desc = "Insert CNVs from seg file into database"
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("segfile", help="titan *_segs.txt file")
    parser.add_argument("library_id", help="library ID")
    parser.add_argument("database", help="database name")
    parser.add_argument("host", help="database host")
    parser.add_argument("user", help="database user")
    parser.add_argument("password", help="database password")
    args = parser.parse_args()
    
    db = cancerGenome.cancerGenomeDB(
        database_name=args.database,
        database_host=args.host,
        database_user=args.user,
        database_password=args.password
    )
    cnv_type = "somatic"
    cursor = db.db.cursor()
    reader = csv.DictReader(open(args.segfile), delimiter="\t")
    next(reader) # skip header
    for row in reader:
        sample = row["Sample"]
        chr = row["Chromosome"]
        start = row["Start_Position(bp)"]
        end = row["End_Position(bp)"]
        copy_number = row["Copy_Number"]
        loh_state = row["TITAN_call"]
        addCNV(cursor, args.library_id, chr, start, end, copy_number, cnv_type)
        addLOH(cursor, args.library_id, chr, start, end, copy_number, loh_state)
