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
# import os
from Bio import Entrez
#from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import GenBank
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
        help="The output directory name"
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
    return(seq_id)


def blast_entries(genbank_entries):
    rec = []
    for record in genbank_entries:
        try:
            result_handle = NCBIWWW.qblast("blastn", "nt", record)
            blast_records = NCBIXML.read(result_handle)
            b_parser = NCBIWWW.BlastParser()
            b_record = b_parser.parse(blast_records)
            rec.append(b_record)
            print(b_record)
            import pdb; pdb.set_trace()
        except:
            pass

    return rec


def write_file(blast_results):
    pass

def main():
    Entrez.email = 'gcagle1@lsu.edu'
    args = get_args()
    org_id = get_tax_id(args)
    genbank_entries = get_genbank_entries(org_id)
    blast_results = blast_entries(genbank_entries)

if __name__ == '__main__':
    main()
