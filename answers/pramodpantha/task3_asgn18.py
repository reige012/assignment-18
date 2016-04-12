#!/usr/bin/env python
# utf-8


"""
task 3 of assignment 18
Created working on group by Pramod Pantha and Mukesh Maharjan on April 10, 2016
Copyright 2016. All right reserved.
reference:
http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc65
https://docs.python.org/3.5/tutorial/errors.html
class slide
"""

from Bio.Blast import NCBIWWW
import argparse
import time
from Bio import Entrez
Entrez.email = "ppanth1@lsu.edu"
organism = 'Loxosceles reclusa'


def get_fastafile(organisname, output):
    esearch_query = Entrez.esearch(db="nucleotide", term=organism, retmode="xml")
    esearch_result = Entrez.read(esearch_query)
    bb_seq_id = esearch_result['IdList']
    time.sleep(1)
    # limits the speed of request
    filename = 'gi_'
    # it will extract the GI list for all IDs of the given species
    for gi in bb_seq_id:
        try:
            result_handle = NCBIWWW.qblast("blastn", "nt", gi, format_type='Text')
            print('Result Obtained')
            output = result_handle.read()
        except ValueError:
            output = ''
        # try except will allow to run program even though it prints ValueError
        # for some of the IDs, in our case second ID have NNN as a sequence.so
        # it will be useful
        save_file = open(filename+gi+'.txt', "w")
        save_file.write(output)
        save_file.close()
        result_handle.close()


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("outputfiledirectory", help="give output directory")
    args = parse.parse_args()
    fastafile = args.outputfiledirectory
    get_fastafile(fastafile)

if __name__ == '__main__':
    main()
