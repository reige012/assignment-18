#!/usr/bin/env python
# utf-8

"""
Assignment 18, Task 2
Write a program that queries GenBank for all the sequence data for Loxosceles reclusa and then writes the sequence data to a file in FASTA format. Use argparse to take the name of the organism as input and the name of some FASTA file as output. Be sure that your program is formatted correctly (PEP8). Also be sure that you pass your email to NCBI (Entrez.email) and that you limit your requests to 1 per second."""

import argparse
from Bio import Entrez
from Bio import SeqIO
import time

Entrez.email = 'jonnations@gmail.com'


def in_n_out():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=True, type=str,
                        help="give input taxon name for NCBI to search.")
    parser.add_argument('-o', '--outfile', required=True, help="give output file (.fasta format) where GenBank data is going.")
    return parser.parse_args()


def search(args):
    time.sleep(1)
    taxon = args.input
    search_query = Entrez.esearch(db="taxonomy", term=taxon, retmode="xml")
    search_result = Entrez.read(search_query)
    # print(search_result)
    speciesID = search_result['IdList']
    # speciesID = 'txid{}'.format(speciesID[0])
    nuc_search = Entrez.esearch(db='nucleotide', term=speciesID, retmode='xml')
    nuc_result = Entrez.read(nuc_search)
    # print(nuc_result)
    nucs = nuc_result['IdList']
    nuc_fetch = Entrez.efetch(db='nucleotide', id=nucs, rettype='gb',
                              retmode='text')
    fetch_results = SeqIO.read(nuc_fetch, 'genbank')
    # print(fetch_results)
    seqs = fetch_results
    return seqs


def conversion(args, seqs):
    with open(args.outfile, 'w') as outfile:
        outfile.write(seqs.format('fasta'))


def main():
    args = in_n_out()
    seqs = search(args)
    conversion(args, seqs)

if __name__ == '__main__':
    main()
