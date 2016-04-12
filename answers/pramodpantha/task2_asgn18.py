#!/usr/bin/env python
# utf-8


"""
task 2 of assignment 18
Created working on group by Pramod Pantha and Mukesh Maharjan on April 10, 2016
Copyright 2016. All right reserved.
reference:
http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc65
http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc23
class slide
"""


from Bio import SeqIO
import time
from Bio import Entrez
import argparse
Entrez.email = "ppanth1@lsu.edu"


def get_fastafile(organism, outputfile):
    esearch_query = Entrez.esearch(db="nucleotide", term=organism, retmode="xml")
    esearch_result = Entrez.read(esearch_query)
    bb_seq_id = esearch_result['IdList']
    all_records = []
    for record in bb_seq_id:
        genbank_record = Entrez.efetch(db="nucleotide", id=record, rettype="gb", retmode="text")
        records = SeqIO.read(genbank_record, 'genbank')
        all_records.append(records)
        time.sleep(1)
    SeqIO.write(all_records, outputfile, "fasta")


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("genus", help="give name of your genus")
    parse.add_argument("species", help="give name of your species")
    parse.add_argument("outputfile", help="give your output fasta file name")
    args = parse.parse_args()
    organismname = args.genus + ' ' + args.species
    fastafile = args.outputfile
    get_fastafile(organismname, fastafile)


if __name__ == '__main__':
    main()
