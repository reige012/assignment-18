#!/usr/bin/env python
# utf-8

"""
BIOL 7800 Assignment 18, Task 2
Austen T. Webber
2016_4_12
"""

from Bio import Entrez
import argparse
import os
import sys
from Bio import SeqIO
import time


def askingforfiles():
    parser = argparse.ArgumentParser(
        description="get species name and directory for fasta output")
    parser.add_argument(
        "--directory",
        required=True,
        help="Directory path for output .fasta file",
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
    newdir = os.path.join(directory_file, 'Sequence')
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
        time.sleep(1)
        genbank_record = Entrez.efetch(db="nucleotide",
                                       id=record, rettype="gb", retmode="text")
        record = SeqIO.read(genbank_record, 'genbank')
        Seq_output = open("query_seq.fasta", "a")
        SeqIO.write(record, Seq_output, "fasta")
    return test


def main():
    path = askingforfiles()
    X = makedirectory(path.directory)
    os.chdir(X)
    taxid = find_it(path.name)


if __name__ == '__main__':
    main()
