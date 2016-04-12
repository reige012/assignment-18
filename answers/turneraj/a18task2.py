#!/usr/bin/env python
# encoding: utf-8

"""
My second task for Assignment 18.

Created by A.J. Turner on April 11, 2016
Copyright 2016 A.J. Turner. All rights reserved. Collaboration/helpful hints
provided by Mikey Henson and Austen Webber.
"""
from Bio import Entrez
import argparse
import os
from Bio import SeqIO
import time


def file_info():
    """input species name and create name of output file directory """
    parser = argparse.ArgumentParser()
    parser.add_argument("--org_in", help="name of organism for search", type=str)
    parser.add_argument("--path_out", help="Director for output file", type=str)
    return parser.parse_args()


def create_directory(directory):
    new = os.path.join(directory, 'GenBank Sequences')
    if not os.path.exists(new):
        os.makedirs(new)
    return new


def genbank_info(species):
    all_seq = []
    info = Entrez.esearch(db="nucleotide", term=species, retmode="xml")
    record = Entrez.read(info)
    get_id = record['IdList']
    for the_id in get_id:
        genbank = Entrez.efetch(db="nucleotide", id=the_id, rettype="gb", retmode="text")
        record = SeqIO.read(genbank, 'genbank')
        all_seq.append(record)
    return all_seq


def main():
    time.sleep(1)
    Entrez.email = "aturn59@lsu.edu"
    path = file_info()
    new_directory = create_directory(path.path_out)
    os.chdir(new_directory)
    organism = genbank_info(path.org_in)
    Seq_output = open("seq.brown-recluse", "w")
    SeqIO.write(organism, Seq_output, "fasta")


if __name__ == '__main__':
    main()
