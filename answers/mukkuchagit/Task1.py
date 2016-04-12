#!/usr/bin/env python
# encoding: utf-8
"""
assignment 18.

Copyright 2016 Mukesh Maharjan. All rights reserved.
"""
import time
import argparse
import csv
from Bio import Entrez
Entrez.email = "mmahar4@lsu.edu"


def get_data(species_name):
    esearch_query = Entrez.esearch(db="taxonomy",
                                   term=species_name, retmode="xml")
    esearch_result = Entrez.read(esearch_query)
    bb_id = esearch_result['IdList']
    # print(species_name)

    data = []
    names = ['superclass', 'class', 'subclass', 'infraclass',
             'superorder', 'order', 'superfamily', 'family', 'genus']

    if bb_id:
        genbank_entries = Entrez.efetch(db="taxonomy", id=bb_id, retmode="xml")
        esearch_result = Entrez.read(genbank_entries)
        for i in range(len(names)):
            found = 0
            for d in esearch_result[0]['LineageEx']:
                if d['Rank']== names[i]:
                    data.append(d['ScientificName'])
                    found = 1
            if found == 0:
                data.append(' ')
    else:
        for j in range(9):
            data.append(' ')
    return data


def save_taxonomy(outfilename, infilename):
    writefile = open(outfilename, 'w')
    writer = csv.DictWriter(writefile, fieldnames=['superclass','class','subclass','infraclass','superorder','order','superfamily',
                                                   'family','genus'])
    writer.writeheader()
    f1 = []
    # datawriter = csv.writer(writefile, dialect='mydialect')
    with open(infilename, 'r') as readfile:
        for line in readfile:
            f1.append(line[:-1])
            row_name = (get_data(line[:-1]))
            writer.writerow({'superclass':row_name[0],'class':row_name[1],'subclass':row_name[2],'infraclass':row_name[3],'superorder':row_name[4],'order':row_name[5],
                             'superfamily':row_name[6],'family':row_name[7],'genus':row_name[8]})
            time.sleep(1)

    writefile.close()
    readfile.close()


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("input_filename", help="Give the" +
                       "filename of the list of organisms")
    parse.add_argument("output_filename",
                       help="Give the name of the output file.")

    file = parse.parse_args()
    in_file1 = file.input_filename
    infilename = in_file1 + '.txt'
    out_file1 = file.output_filename
    outfilename = out_file1 + '.csv'

    save_taxonomy(outfilename, infilename)

if __name__ == '__main__':
    main()
