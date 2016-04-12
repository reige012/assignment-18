#!/usr/bin/env python
# encoding: utf-8
"""
assignment 18.

Copyright 2016 Mukesh Maharjan. All rights reserved.
"""
import os
import argparse
from Bio.Blast import NCBIWWW
import time
from Bio import Entrez
Entrez.email = "mmahar4@lsu.edu"


def get_blast_result(organism, direcory_name):
    esearch_query = Entrez.esearch(db="nucleotide",
                                   term=organism, retmode="xml")
    esearch_result = Entrez.read(esearch_query)
    bb_seq_id = esearch_result['IdList']

    for gi in bb_seq_id:
        filename = direcory_name+"/"

        try:
            result_handle = NCBIWWW.qblast("blastn", "nt", gi,
                                           format_type='Text')
            output = result_handle.read()
        except ValueError:
            output = ''
        filename = filename+'gi_'+gi+'.txt'
        os.makedirs(os.path.dirname(filename), exist_ok=True)

        save_file = open(filename, "w")
        save_file.write(output)
        time.sleep(1)

        save_file.close()
        result_handle.close()


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("organism_name", help="Give the name of the organism" +
                       " strictly inside \" \" ")
    parse.add_argument("output_directory",
                       help="Give the name of the output file.")

    file = parse.parse_args()
    org_name = file.organism_name
    out_file = file.output_directory

    get_blast_result(org_name, out_file)


if __name__ == '__main__':
        main()
