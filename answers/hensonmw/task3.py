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
import time
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


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
        help="Provide the scientific name of which you want to search",
        type=str
    )
    return parser.parse_args()


def makedirectory(directory):
    new = os.path.join(directory, 'Requested Sequence')
    if not os.path.exists(new):
        os.makedirs(new)
    return new


def get_it(species):
    query = Entrez.esearch(db = "nucleotide", term = species, retmode = "xml")
    record = Entrez.read(query)
    seq_id = record['IdList']
    return seq_id


def blastin(GI):
    results = []
    for one_GI in GI:
        try:
            result = NCBIWWW.qblast("blastn", "nr", one_GI)
        except:
            pass
        results.append(result)
    return results


def main():
    time.sleep(1)
    Entrez.email = "mhenso4@lsu.edu"
    path = askingforfiles()
    X = makedirectory(path.D)
    os.chdir(X)
    input_GI = get_it(path.I)
    blasting = blastin(input_GI)
    parser = NCBIWWW.BlastParser()
    record = parser.parse(blasting)
    with open("blast_results.txt", 'w') as fun:
        fun.write(record)


if __name__ == '__main__':
    main()
