# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
BIOL7800 assignment 18
Oscar Johnson 11 April 2016

Copyright Oscar Johnson 2016

provide input txt file of species names
returns NCBI taxonomy as csv file
"""

import os
import csv
import time
import argparse
from Bio import Entrez


def get_args():
    parser = argparse.ArgumentParser(
            description="""provide input txt file""")
    parser.add_argument('--in_file',
                        type=str,
                        required=True,
                        help="enter a .txt file",
                        )
    parser.add_argument('--out_file',
                        type=str,
                        required=True,
                        help="enter an output file name for csv file",
                        )
    return parser.parse_args()


def ncbi(args):
    """
    takes data from input txt file
    queries NCBI taxonomy for lineage information
    and writes the results to a .csv file
    """
    full_tax = ['superclass', 'class', 'subclass',
                'infraclass', 'superorder', 'order',
                'superfamily', 'family', 'genus']
    with open(args.out_file, 'w') as outfile:
        writer = csv.writer(outfile, delimiter=',',
                            quotechar='|',
                            quoting=csv.QUOTE_MINIMAL)
        writer.writerow(full_tax)  # add header row
        writer = csv.DictWriter(outfile, fieldnames=full_tax)
        with open(args.in_file, 'r') as infile:
            for sp in infile:
                species = sp.strip('\n')  # remove newlines
                print('i am working on species: {}'.format(species))
                esearch_query = Entrez.esearch(db="taxonomy",
                                               term=species,
                                               retmode="xml")
                esearch_result = Entrez.read(esearch_query)
                time.sleep(1)
                # print(esearch_result)
                for iden in esearch_result['IdList']:
                    tax_entry = Entrez.efetch(db="taxonomy",
                                              id=iden,
                                              retmode="xml")
                    result = Entrez.read(tax_entry)
                    taxon_dict = {}
                    for rank in result[0]['LineageEx']:
                        if rank["Rank"] in full_tax:
                            taxon_dict[rank["Rank"]] = rank["ScientificName"]
                    writer.writerow(taxon_dict)


def main():
    args = get_args()
    # rename outfile to .csv if needed
    if args.out_file[-4:] != '.csv':
        args.out_file += '.csv'
    else:
        pass
    Entrez.email = "ojohns7@lsu.edu"
    os.chdir(os.path.dirname(args.in_file))
    ncbi(args)

if __name__ == '__main__':
    main()
