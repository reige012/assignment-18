#!/usr/bin/env python
# encoding: utf-8

"""
Assignment 18
Task 2 Program: Writing sequence data from GenBank to a FASTA file
Jessie Salter
9 April 2016
"""

import argparse
from Bio import Entrez
from Bio import SeqIO
import time


def parser():
    '''Takes user input'''
    parser = argparse.ArgumentParser(
        description='''Gets the path to the input file''')
    parser.add_argument(
        "--orgn",
        required=True,
        help="Enter the name of the species you want to search for",
        type=str
        )
    parser.add_argument(
        "--output",
        required=True,
        help="Enter the name of the output fasta file you want to create",
        type=str
        )
    return parser.parse_args()


def esearch(orgn):
    '''Collects the Identifiers from NCBI'''
    query = Entrez.esearch(db="nucleotide", term=orgn, retmode="xml")
    result = Entrez.read(query)
    iden = result['IdList']
    return iden


def efetch(id_list):
    '''Collects all sequence data for the given organism from GenBank'''
    seqs = []
    for iden in id_list:
        query = Entrez.efetch(
            db='nucleotide',
            id=iden,
            rettype='gb',
            retmode='text'
            )
        data = SeqIO.read(query, 'genbank')
        seqs.append(data)
        time.sleep(1)
    return seqs


def results_writer(seqs, output):
    '''Writes the sequences to a FASTA file'''
    with open(output, 'w') as outfile:
        SeqIO.write(seqs, outfile, 'fasta')


def main():
    Entrez.email = "jsalte5@lsu.edu"
    args = parser()
    id_list = esearch(args.orgn)
    seqs = efetch(id_list)
    results_writer(seqs, args.output)


if __name__ == '__main__':
    main()
