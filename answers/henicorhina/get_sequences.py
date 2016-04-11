# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BIOL7800 assignment 18
Oscar Johnson 11 April 2016

Copyright Oscar Johnson 2016

provide input txt file of species names
returns NCBI taxonomy as csv file
"""

import os
import time
import argparse
from Bio import Entrez
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser(
            description="""get sequence data for a species""")
    parser.add_argument('--species',
                        type=str,
                        required=True,
                        help="enter a genus + species name",
                        )
    parser.add_argument('--out_file',
                        type=str,
                        required=True,
                        help="enter an output file name for fasta file",
                        )
    parser.add_argument('--out_dir',
                        type=str,
                        required=True,
                        help="enter an output directory for fasta file",
                        )
    return parser.parse_args()


def ncbi_fasta(args):
    """
    takes input species name and outfile
    queries NCBI for all sequence data
    and writes the results to a fasta file
    """
    esearch_query = Entrez.esearch(db="taxonomy",
                                   term=args.species,
                                   retmode="xml")
    esearch_result = Entrez.read(esearch_query)
    time.sleep(1)
    # print(esearch_result)
    ident = "txid{}".format(esearch_result['IdList'][0]) 
    genbank_entries = Entrez.esearch(db="nucleotide",
                                     term=ident,
                                     retmode="xml")
    esearch_seqs = Entrez.read(genbank_entries)
    time.sleep(1)
    for record in esearch_seqs['IdList']:    
        seq_entry = Entrez.efetch(db="nucleotide",
                                  id=record,
                                  rettype="gb",
                                  retmode="text")
        time.sleep(1)
        record = SeqIO.read(seq_entry, 'genbank')


def main():
    args = get_args()
    # rename outfile to .csv if needed
    if args.out_file[-6:] != '.fasta':
        args.out_file += '.fasta'
    else:
        pass
    Entrez.email = "ojohns7@lsu.edu"
    os.chdir(os.path.dirname(args.out_dir))
    args.species = args.species.lower()
    ncbi_fasta(args)

if __name__ == '__main__':
    main()
