#!/usr/bin/env python
# utf-8

"""
BIOL 7800 Assignment 18, Task 3
Austen T. Webber
2016_4_12
"""

from Bio import Entrez
import argparse
import os
import sys
from Bio import SeqIO
import time
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


def askingforfiles():
    parser = argparse.ArgumentParser(
        description="get name and output directory")
    parser.add_argument(
        "--directory",
        required=True,
        help="Directory path for output .txtfile",
        type=str
    )
    parser.add_argument(
        "--name",
        required=True,
        help="What is the scientific name of the organism?",
        type=str
    )
    return parser.parse_args()


def makedirectory(directory_file):
    newdir = os.path.join(directory_file, 'BLAST')
    print(newdir)
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    return newdir


def find_it(species):
    test = []
    Entrez.email = "awebbe4@lsu.edu"
    esearch_query = Entrez.esearch(db="nucleotide",
                                   term=species, retmode="xml")
    esearch_record = Entrez.read(esearch_query)
    print(esearch_record)
    bb_seq_id = esearch_record['IdList']
    print(bb_seq_id)
    for record in bb_seq_id:
        genbank_record = Entrez.efetch(db="nucleotide",
                                       id=record, rettype="gb", retmode="text")
        record = SeqIO.read(genbank_record, 'genbank')
        result_handle = NCBIWWW.qblast("blastn", "nt", record.seq)
        blast_records = NCBIXML.parse(result_handle)
        save_file = open("blast_file.txt", "a")
        save_file.write(result_handle.read())
        save_file.close()
        result_handle.close()
    return test


def main():
    time.sleep(1)
    path = askingforfiles()
    X = makedirectory(path.directory)
    os.chdir(X)
    taxid = find_it(path.name)

if __name__ == '__main__':
    main()
