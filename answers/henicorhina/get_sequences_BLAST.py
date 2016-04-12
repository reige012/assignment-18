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
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


def get_args():
    parser = argparse.ArgumentParser(
            description="""get sequence data for a species
                        and blast against GenBank nucleotide database""")
    parser.add_argument('--species',
                        type=str,
                        required=True,
                        help="enter a genus + species name enclosed in parentheses",
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
    then BLASTs them against NCBI
    and writes the results to an xml file
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
        records = SeqIO.read(seq_entry, 'genbank')
        print('i am working on record #: {}'.format(records.name))
        if records.seq[0] == 'N':
            pass
        else:
            result = NCBIWWW.qblast('blastn',
                                    'nt',
                                    records.format("fasta"))
            # blast_records = NCBIXML.read(result)
            name = "{}.xml".format(records.name)
            with open(name, 'w') as outfile:
                outfile.write(result.read())


def main():
    args = get_args()
    Entrez.email = "ojohns7@lsu.edu"
    os.chdir(os.path.abspath(args.out_dir))
    args.species = args.species.lower()
    ncbi_fasta(args)

if __name__ == '__main__':
    main()
