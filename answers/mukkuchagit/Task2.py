#!/usr/bin/env python
# encoding: utf-8
"""
assignment 18.

Copyright 2016 Mukesh Maharjan. All rights reserved.
"""
from Bio import SeqIO
import argparse
import time
from Bio import Entrez
Entrez.email = "mmahar4@lsu.edu"


def get_seqeuence_data(organism, output_filename):
    # organism = 'Loxosceles reclusa'
    # output_filename= "mk4.fasta"
    searchquery = Entrez.esearch(db="nucleotide", term=organism, retmode="xml")
    esearch_result = Entrez.read(searchquery)
    bb_seq_id = esearch_result['IdList']
    all_records = []
    for record in bb_seq_id:
        genbank_record = Entrez.efetch(db="nucleotide", id=record,
                                       rettype="gb", retmode="text")
        records = SeqIO.read(genbank_record, 'genbank')
        all_records.append(records)
        time.sleep(1)
    SeqIO.write(all_records, output_filename, "fasta")


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("organism_name", help="Give the name of the organism" +
                       " strictly inside \" \" ")
    parse.add_argument("outfile_name",
                       help="Give the name of the output file.")

    file = parse.parse_args()
    print(file.organism_name)
    org_name = file.organism_name
    out_file1 = file.outfile_name
    out_file = out_file1 + '.fasta'

    get_seqeuence_data(org_name, out_file)


if __name__ == '__main__':
    main()
