#!/usr/bin/env python

# Import a titan *_segs.txt file containing CNVs and LOH segments into the
# database.

import argparse
import cancer_db as cancerGenome
import csv
import itertools

if __name__ == "__main__":
    desc = "Insert CNV and LOH from TITAN *_segs.txt into database"
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

    cursor = db.db.cursor()
    query = "SELECT library_name FROM library WHERE id = %s"
    cursor.execute(query, args.library_id)
    library_name = cursor.fetchone()[0]
    print(args.library_id)
    db.loadCNVs(args.segfile, library_name=library_name, file_format="titan")

    reader = csv.DictReader(open(args.segfile), delimiter="\t")
    next(reader) # skip header

    params = {"library_id": args.library_id}

    filter_fun = lambda r: r["TITAN_call"] != "HET"
    for row in itertools.ifilter(filter_fun, reader):
        params.update({"chromosome": row["Chromosome"],
                       "segment_start": row["Start_Position(bp)"],
                       "segment_end": row["End_Position(bp)"],
                       "size": row["Length(bp)"],
                       "bin_count": 0, # don't know what this is
                       "copy_number": row["Copy_Number"],
                       "loh_state": row["TITAN_call"],
                       "Major_allele_count": row["MajorCN"],
                       "minor_allele_count": row["MinorCN"]})
        db.addLohSegment(**params)
