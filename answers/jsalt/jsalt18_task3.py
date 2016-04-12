#!/usr/bin/env python
# encoding: utf-8

"""
Assignment 18
Task 3 Program: GenBank search + blastn
Jessie Salter
9 April 2016
"""

import argparse
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import time
import os.path


def parser():
    '''Takes user input'''
    parser = argparse.ArgumentParser(
        description='''Gets the path to the input file''')
    parser.add_argument(
        "--orgn",
        required=True,
        help="Enter the name of the species you want to search for",
        type=str
        )
    parser.add_argument(
        "--output",
        required=True,
        help="Enter the path to the output directory",
        type=str
        )
    return parser.parse_args()


def esearch(orgn):
    '''Collects the Identifiers from NCBI'''
    query = Entrez.esearch(db="nucleotide", term=orgn, retmode="xml")
    result = Entrez.read(query)
    iden = result['IdList']
    return iden


def efetch(id_list):
    '''Collects all sequence data for the given organism from GenBank'''
    seqs = []
    for iden in id_list:
        query = Entrez.efetch(
            db='nucleotide',
            id=iden,
            rettype='gb',
            retmode='text'
            )
        data = SeqIO.read(query, 'genbank')
        seqs.append(data)
        time.sleep(1)
    return seqs


def gi_records(seqs):
    '''Creates a list of GI record numbers from the GenBank sequences'''
    gis = []
    for record in seqs:
        gis.append(record.id)
    gis.pop(1)
    return gis


def blast(gis):
    '''Performs a BLAST search using the GI records, and returns an iterator of
    the blastn results for each species in a list'''
    results = []
    for record in gis:
        # Something in here is not working, not sure what the exact problem is,
        # I keep getting ID#32 error message: "query contains no sequence data"
        blastn = NCBIWWW.qblast('blastn', 'nt', str(record))
        blastn_results = NCBIXML.parse(blastn)
        results.append(blastn_results)
        time.sleep(1)
    return results


def results_writer(results, output, gis):
    '''Parses the results of the blastn searches and writes an output text file
    for each blastn result'''
    for x in range(len(results)):
        # This generates a new filename using the GI number:
        filename = str(gis[x])+'.txt'
        newpath = os.path.join(output, filename)
        with open(newpath, 'w') as outfile:
            for record in results[x]:
                outfile.write(record)


def main():
    Entrez.email = "jsalte5@lsu.edu"
    args = parser()
    id_list = esearch(args.orgn)
    seqs = efetch(id_list)
    gis = gi_records(seqs)
    results = blast(gis)
    results_writer(results, args.output, gis)


if __name__ == '__main__':
    main()
