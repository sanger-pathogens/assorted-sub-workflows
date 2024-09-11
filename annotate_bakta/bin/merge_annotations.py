#!/usr/bin/python3

import gffutils
import argparse
import os
from Bio import SeqIO
from io import StringIO

#inspired from general practices from here: https://biopython.org/DIST/docs/tutorial/Tutorial.html

def convert_genbank_to_gff3(genbank_file):
    gff3_output = StringIO()
    for record in SeqIO.parse(genbank_file, "genbank"):
        for feature in record.features:
            if feature.type in ['gene', 'CDS', 'mRNA', 'exon', 'intron']:
                gff3_output.write(f"{record.id}\tGenBank\t{feature.type}\t{int(feature.location.start)+1}\t{int(feature.location.end)}\t.\t{feature.strand}\t.\tID={feature.id};")
                if 'gene' in feature.qualifiers:
                    gff3_output.write(f"Name={feature.qualifiers['gene'][0]};")
                if 'product' in feature.qualifiers:
                    gff3_output.write(f"product={feature.qualifiers['product'][0]};")
                gff3_output.write("\n")
    gff3_output.seek(0)
    return gff3_output

def merge_gff_gtf(main_file, additional_files, output, db_dir=':memory:'):

    # Run in memory mode
    if db_dir == ':memory:':
        db = gffutils.create_db(main_file, dbfn=db_filename, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    
    #or not
    else:
        # Make a directory to store the DB
        if not os.path.exists(db_dir):
            os.makedirs(db_dir)
        # Name the DB after the input file
        db_filename = os.path.join(db_dir, os.path.splitext(os.path.basename(main_file))[0] + "_db")

        # add something so it can load if its there (someone else using this)
        if not os.path.exists(db_filename):
            db = gffutils.create_db(main_file, dbfn=db_filename, force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
        else:
            db = gffutils.FeatureDB(db_filename)



    # Merge the additional GFF/GTF files into the existing database

    #could do better than extention maybe we look for LOCUS or something in the genbank file
    for file in additional_files:
        if file.endswith('.gbk') or file.endswith('.gb'):
            # Convert GenBank to GFF3 format
            gff3_file = convert_genbank_to_gff3(file)
            db.update(gffutils.DataIterator(gff3_file), merge_strategy='merge')
        else:
            db.update(gffutils.DataIterator(file), merge_strategy='merge')

    # Write the merged database to a GFF3 file
    with open(output, 'w') as gff_out:
        for feature in db.all_features(order_by=None):
            gff_out.write(str(feature) + '\n')


if __name__ == "__main__":
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description="Merge multiple GFF/GTF files into a single GFF3 file using a main GFF as the base.")
    parser.add_argument('main_file', type=str, help='Main GFF file to use as the base for the database')
    parser.add_argument('additional_files', metavar='F', type=str, nargs='+', help='List of additional GFF/GTF files to merge')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output merged GFF3 file')
    parser.add_argument('-d', '--dbdir', type=str, help='Directory to save the database file. leave blank for an in-memory database.')

    args = parser.parse_args()

    merge_gff_gtf(args.main_file, args.additional_files, args.output, db_dir=args.dbdir)
