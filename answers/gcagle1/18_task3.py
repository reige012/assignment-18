#! /usr/bin/env python
# encoding: utf-8

'''
Grace Cagle
Assignment 18 Task 3
Using task 2 to perform a BLAST search against the nucleotide database
and writing output text files to a directory
'''

import argparse
# import glob
import os
from Bio import Entrez
#from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import time
#import csv


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-n', "--name",
        required=True,
        nargs='+',
        help="The name of the taxa to acquire sequence data for"
    )
    parser.add_argument(
        '-o', "--outfile",
        required=True,
        help="The output directory name"
    )
    return parser.parse_args()


def get_tax_id(args):
    '''searches NCBI taxonomy for the input organism'''
    esearch_query = Entrez.esearch(db='taxonomy', term=args.name, retmode='xml')
    esearch_result = Entrez.read(esearch_query)
    org_id = esearch_result['IdList']
    return org_id


def get_genbank_entries(org_id):
    '''returns a list of GIs'''
    org_id2 = "txid{}".format(org_id[0])
    genbank_entries = Entrez.esearch(db='nucleotide', term=org_id2, retmode='xml')
    esearch_result = Entrez.read(genbank_entries)
    seq_id = esearch_result['IdList']
    return seq_id


def blast_entries(genbank_entries):
    '''blasts GIs and returns a list of records'''
    records = []
    for record in genbank_entries:
        time.sleep(1)
        try:
            result_handle = NCBIWWW.qblast("blastn", "nt", record)
            blast_record = NCBIXML.read(result_handle)
            records.append(blast_record)
        except:
            pass
    return records


def write_file(args, records):
    for rec in records:
        filename = (os.path.join(args.outfile, str(rec.query_id)))
        with open(filename, 'w') as f:
            f.write(rec)


def main():
    Entrez.email = 'gcagle1@lsu.edu'
    args = get_args()
    org_id = get_tax_id(args)
    genbank_entries = get_genbank_entries(org_id)
    records = blast_entries(genbank_entries)
    write_file(args, records)


if __name__ == '__main__':
    main()
