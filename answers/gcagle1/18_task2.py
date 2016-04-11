#! /usr/bin/env python
# encoding: utf-8

'''
Grace Cagle
Assignment 18 Task 2
A proram taking the name of an organism, querying genbank for sequence data
for the organism, and writing the sequence data to a fasta file of a name
specified by the user.
'''

import argparse
# import glob
# import os
from Bio import Entrez
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord
import time
import csv


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
        help="The output file name"
    )
    return parser.parse_args()


def get_tax_id(args):
    esearch_query = Entrez.esearch(db='taxonomy', term=args.name, retmode='xml')
    esearch_result = Entrez.read(esearch_query)
    org_id = esearch_result['IdList']
    return org_id


def get_genbank_entries(org_id):
    entries = []
    org_id2 = "txid{}".format(org_id[0])
    genbank_entries = Entrez.esearch(db='nucleotide', term=org_id2, retmode='xml')
    esearch_result = Entrez.read(genbank_entries)
    seq_id = esearch_result['IdList']
    for record in seq_id:
        time.sleep(1)
        genbank_record = Entrez.efetch(db='nucleotide', id=record, rettype='gb', retmode='text')
        record = SeqIO.read(genbank_record, 'genbank')
        entries.append(record.seq)
    return(entries)


def write_file(args, genbank_entries):
    with open(args.outfile, 'w') as f:
        for entry in genbank_entries:
            seq_record_object = SeqRecord.SeqRecord(entry)
            f.write(seq_record_object.format('fasta'))


def main():
    Entrez.email = 'gcagle1@lsu.edu'
    args = get_args()
    org_id = get_tax_id(args)
    genbank_entries = get_genbank_entries(org_id)
    write_file(args, genbank_entries)

if __name__ == '__main__':
    main()
