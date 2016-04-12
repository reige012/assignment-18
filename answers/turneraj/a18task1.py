#!/usr/bin/env python
# encoding: utf-8

"""
My first task for Assignment 18.

Created by A.J. Turner on April 10, 2016
Copyright 2016 A.J. Turner. All rights reserved. Collaboration/helpful hints
provided by Mikey Henson and Austen Webber. Other information on task provided
by:
http://stackoverflow.com/questions/16504238/attempting-to-obtain-taxonomic-information-from-biopython

"""
import argparse
from Bio import Entrez
import time


def file_info():
    """input file with species names """
    parser = argparse.ArgumentParser()
    parser.add_argument("--path_in", help="path of input file", type=str)
    parser.add_argument("--path_out", help="Directory (folder) for output file", type=str)
    return parser.parse_args()


def my_tax_id(my_species):
    all_species = []  # creating empty list for taxonomic info to parse into
    for species in my_species:
        species = species.replace(' ', "+")
        info = Entrez.esearch(term=species, db="taxonomy", retmode="xml")
        record = Entrez.read(info)
        species_id = record['IdList']
        all_species.append(species_id)
    return all_species


def my_tax_info(tax_id):
    """getting records using the taxonomic IDs that were extracted"""
    all_taxonomy = []
    for the_id in tax_id:
        search = Entrez.efetch(id=the_id, db="taxonomy", retmode="xml")
        taxonomy = Entrez.read(search)
        lineage = {d['Rank']: d['ScientificName'] for d in taxonomy[0]['LineageEx']}
        all_taxonomy.append(lineage)
    return all_taxonomy


def main():
    Entrez.email = "aturn59@lsu.edu"
    time.sleep(1)
    my_path = file_info()
    my_species = open(my_path.path_in)
    tax_id = my_tax_id(my_species)
    # print(tax_id) # to see if ids were retrieved
    tax_data = my_tax_info(tax_id)
    # print(tax_data) # to see if you retrieved it
    specifiy_it = ['superclass', 'class', 'subclass', 'infraclass', 'superorder','order', 'superfamily', 'family', 'genus']
    with open("test.csv", "w") as tests:
        for species in specifiy_it:
            tests.write(species+",")
        tests.write("\n")
        for taxa in tax_data:
            for species in specifiy_it:
            # print(species) to make sure accessing above list items
                try:
                    tests.write(taxa[species]+",")
                except:
                    tests.write("NA"+",")
            tests.write("\n")


if __name__ == '__main__':
    main()
