#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on Apr 12, 2016 Building on Task 2, modify your program above to get GI
record of each sequence from NCBI (e.g. do not work from the FASTA file you
created) and the send those records to NCBI BLAST to perform a blastn search of
each sequence against the nt (nucleotide) database. Write each BLAST result (the
result for each sequence, which may consist of several hits) out to text files
in BLAST format (use argparse to get the output directory name from the user).
@author: York
'''


from Bio import Entrez
import argparse
import time
from Bio import SeqIO


Entrez.email = "catchdonovan@gmail.com"


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputdir", required=True)
    parser.add_argument("--outputdir", required=True)
    args = parser.parse_args()
    return args


def esearch(spider):
    genbank_data = Entrez.esearch(db="nucleotide", term=spider, retmode="xml")
    data = Entrez.read(genbank_data)
    print(data)
    return data


def spi_fetch(record):
    genbank_record = Entrez.efetch(
        db="nucleotide", id=record, rettype="gb", retmode="text")
    data2 = SeqIO.read(genbank_record, 'genbank')
    return data2


def main():
    args = get_parser()
    outfile = open(args.outputdir, "w+")
    inputspi = args.inputdir
    search = esearch(inputspi)
    bb_id = search['IdList']
    print(bb_id)
    for item in bb_id:
        time.sleep(1)
        efetch = spi_fetch(item)
        print(efetch)
        outfile.write(str(efetch))


if __name__ == '__main__':
    main()
