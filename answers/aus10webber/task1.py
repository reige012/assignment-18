#!/usr/bin/env python
# utf-8

"""
BIOL 7800 Assignment 18, Task 1
Austen T. Webber
2016_4_12
"""

import argparse
import os
from Bio import Entrez
import time


def askingforfiles():
    parser = argparse.ArgumentParser(
        description="Get files and directories")
    parser.add_argument(
        "--directory",
        required=True,
        help="Directory path for output .csv file",
        type=str
    )
    parser.add_argument(
        "--spectext",
        required=True,
        help=".txt file containing list of species",
        type=str
    )
    return parser.parse_args()


def makedirectory(directory_file):
    newdir = os.path.join(directory_file, 'Taxonomy Names')
    print(newdir)
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    return newdir


def get_taxID(species):
    species_ID_list = []
    for names in species:
        names = names.replace(' ', "+")
        items = Entrez.esearch(term=names, db="taxonomy", retmode="xml")
        record = Entrez.read(items)
        x = record['IdList']
        species_ID_list.append(x)
    return species_ID_list


def get_taxdata(taxid):
    taxonomy_ids = []
    for ids in taxid:
        search = Entrez.efetch(id=ids, db="taxonomy", retmode="xml")
        tax = Entrez.read(search)
        lineage = {d['Rank']: d['ScientificName'] for d in tax[0]['LineageEx']}
        taxonomy_ids.append(lineage)
    return taxonomy_ids


def main():
    time.sleep(1)
    Entrez.email = "awebbe4@lsu.edu"
    path = askingforfiles()
    X = makedirectory(path.directory)
    os.chdir(X)
    names = open(path.spectext)
    taxID = get_taxID(names)
    taxonomy = get_taxdata(taxID)
    needs = ['superclass', 'class', 'subclass', 'infraclass', 'superorder',
             'order', 'superfamily', 'family', 'genus']
    with open("test.csv", "w") as tests:
        for species in taxonomy:
            for types in needs:
                try:
                    tests.write(species[types] + ",")
                except:
                    tests.write("NA" + ",")
            tests.write("\n")


if __name__ == '__main__':
    main()
