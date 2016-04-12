#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
A18 task2
Created on Apr 12, 2016 Write a program that queries GenBank for all the
sequence data for this spider and then writes the sequence data to a file in
FASTA format. Use argparse to take the name of the organism as input and the
name of some FASTA file as output. 
@author: York
'''


from Bio import Entrez
import argparse
import time
from Bio import SeqIO


Entrez.email = "catchdonovan@gmail.com"


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputname", required=True)
    parser.add_argument("--outputfile", required=True)
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
    outfile = open(args.outputfile, "w+")
    inputspi = args.inputname
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
