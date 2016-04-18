#!/usr/bin/env python
# encoding: utf-8
"""
created by me for task2 to Write a program that queries GenBank for all
the sequence data for this spider and then writes the sequence data to
 a file in FASTA format
"""
import time
import argparse
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "mforoo1@lsu.edu"


def get_parser():
    """
   using argparse to takes the list  as input
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputname", type=str, required=True)
    parser.add_argument("--outputfile", required=True)
    args = parser.parse_args()
    return args


def genebank_sequence(name, outfile):
    esearch_query = Entrez.esearch(db="nucleotide",
                                   term=name,
                                   retmode="xml")
    esearch_result = Entrez.read(esearch_query)
    sequenc_entry = esearch_result['IdList']
    print(sequenc_entry)
    sequence = []
    for i in sequenc_entry:
        result_sequence_entry = Entrez.efetch(db="nucleotide", id=i,
                                              rettype="gb", retmode="text")
        sequences = SeqIO.read(result_sequence_entry, 'genbank')
        print(sequences)
        sequence.append(sequences)
        time.sleep(1)
    # outfile.write("%s" % (sequence))
    SeqIO.write(sequence, outfile, "fasta")


def main():
    args = get_parser()
    outfile = open(args.outputfile, "w")
    organism = args.inputname
    genebank_sequence(organism, outfile)


if __name__ == '__main__':
    main()
