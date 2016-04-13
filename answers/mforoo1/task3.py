#!/usr/bin/env python
# encoding: utf-8
"""
" created by me for task3 to write program to get GI record of each
sequence from NCBI and then send those records to NCBI BLAST to perform a
 blastn search of each sequence
against the nt (nucleotide) database "
"""
import time
import argparse
from Bio import Entrez
from Bio.Blast import NCBIWWW
import os

Entrez.email = "mforoo1@lsu.edu"


def get_parser():
    """
   using argparse to takes the list  as input
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputname", type=str, required=True)
    parser.add_argument("--outputdir", required=True)
    args = parser.parse_args()
    return args


def genebank_sequence(name):
    esearch_query = Entrez.esearch(db="nucleotide",
                                   term=name,
                                   retmode="xml")
    esearch_result = Entrez.read(esearch_query)
    sequenc_entry = esearch_result['IdList']
    print(sequenc_entry)
    for i in sequenc_entry:
        try:
            result_blast = NCBIWWW.qblast("blastn", "nt", i,
                                          format_type='Text')
            output = result_blast.read()
            time.sleep(1)
            with open("outputfile.txt", "a") as outfile:
                outfile.write(output)
        except ValueError:
            output = ''


def main():
    args = get_parser()
    organism = args.inputname
    os.makedirs(args.outputdir)
    os.chdir(args.outputdir)
    genebank_sequence(organism)


if __name__ == '__main__':
    main()
