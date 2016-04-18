#!/usr/bin/env python
# encoding: utf-8

"""
My 1st program for Assignment 18

Created by Michael Henson on 11 April 2016.
Copyright 2016 Michael W Henson. All rights reserved.

"""
import argparse
import os
from Bio import Entrez
import time


def askingforfiles():
    parser = argparse.ArgumentParser(
        description="Input of Species of Interest")
    parser.add_argument(
        "--I",
        required=True,
        help="Provide the desired text file of organisms names ",
        type=str
    )
    parser.add_argument(
        "--O",
        required=True,
        help="Provide the desired output file name",
        type=str
    )
    return parser.parse_args()


def get_taxID(species):
    species_ID_list = []
    for names in species:
        names = names.replace(' ', "+")
        items = Entrez.esearch(term = names, db = "taxonomy", retmode = "xml")
        record = Entrez.read(items)
        x = record['IdList']
        species_ID_list.append(x)
    return species_ID_list


def get_taxdata(taxid):
    taxonomy_ids = []
    for ids in taxid:
        search = Entrez.efetch(id = ids, db = "taxonomy", retmode = "xml")
        tax = Entrez.read(search)
        lineage = {d['Rank']:d['ScientificName'] for d in tax[0]['LineageEx']}
        taxonomy_ids.append(lineage)
    return taxonomy_ids


def main():
    time.sleep(1)
    Entrez.email = "mhens04@lsu.edu"
    path = askingforfiles()
    names = open(path.I)
    taxID = get_taxID(names)
    taxonomy = get_taxdata(taxID)
    needs = ['superclass', 'class', 'subclass', 'infraclass', 'superorder', 'order', 'superfamily', 'family', 'genus']
    #tested = []
    with open(path.O+".csv", "w") as tests:
        for species in taxonomy:
            for types in needs:
                try:
                    tests.write(species[types] + ",")
                except:
                    tests.write("NA" + ",")
            tests.write("\n")


if __name__ == '__main__':
    main()
