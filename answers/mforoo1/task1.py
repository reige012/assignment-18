#!/usr/bin/env python
# encoding: utf-8
"""
created by me for task1 to Write a program to gets the taxonomic information
 from NCBI for each sample, and writes that to a new file as output
"""
import csv
import time
import argparse
from Bio import Entrez

Entrez.email = "mforoo1@lsu.edu"


def get_parser():
    """
   using argparse to takes the list  as input
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputfile", required=True)
    parser.add_argument("--outputfile", required=True)
    args = parser.parse_args()
    return args


def main():
    args = get_parser()
    with open(args.inputfile, "r") as infile:
        for term in infile:
            esearch_query = Entrez.esearch(db="taxonomy", term=term,
                                           retmode="xml")
            esearch_result = Entrez.read(esearch_query)
            # (esearch_result)
            # print(esearch_result["IdList"])
            for iden in esearch_result['IdList']:
                taxonomy_entry = Entrez.efetch(db="taxonomy", id=iden,
                                               retmode="xml")
                taxonomy_entery_result = Entrez.read(taxonomy_entry)
                result = taxonomy_entery_result[0]['LineageEx']
                data = []
                names = ['superclass', 'class', 'subclass', 'infraclass',
                         'superorder', 'order', 'superfamily', 'family',
                         'genus']
                for i in range(len(names)):
                    found = 0
                    for d in result:
                        if d['Rank'] == names[i]:
                            data.append(d['ScientificName'])
                            found = 1
                    if found == 0:
                        data.append(' ')
            else:
                for j in range(9):
                    data.append(' ')
            # print(data)
            with open(args.outputfile, 'w') as myfile:
                wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
                for value in data:
                    wr.writerow(data)
    time.sleep(1)


if __name__ == '__main__':
    main()
