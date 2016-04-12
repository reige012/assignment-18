#!/usr/bin/env python
# utf-8

"""
My 2nd program for Assignment 18

Created by Michael Henson on 11 April 2016.
Copyright 2016 Michael W Henson. All rights reserved.

With help from Austen and AJ

"""

from Bio import Entrez
import argparse
import os
from Bio import SeqIO
import time


def askingforfiles():
    parser = argparse.ArgumentParser(
        description="Looking for the directory and name of the species you want?")
    parser.add_argument(
        "--D",
        required=True,
        help="Where do you want to put the .csv file?",
        type=str
    )
    parser.add_argument(
        "--I",
        required=True,
        help="Provide the scientific name of which you want to search with + inbetween rather than space" ,
        type=str
    )
    return parser.parse_args()


def makedirectory(directory):
    new = os.path.join(directory, 'Requested Sequence')
    if not os.path.exists(new):
        os.makedirs(new)
    return new


def get_it(species):
    test = []
    query = Entrez.esearch(db = "nucleotide", term = species, retmode = "xml")
    record = Entrez.read(query)
    seq_id = record['IdList']
    for item in seq_id:
        genbank = Entrez.efetch(db = "nucleotide", id = item, rettype = "gb", retmode = "text")
        record = SeqIO.read(genbank, 'genbank')
        test.append(record)
    return test


def main():
    time.sleep(1)
    Entrez.email = "mhenso4@lsu.edu"
    path = askingforfiles()
    X = makedirectory(path.D)
    os.chdir(X)
    y = get_it(path.I)
    Seq_output = open("seq.fasta", "w")
    SeqIO.write(y, Seq_output, "fasta")


if __name__ == '__main__':
    main()
