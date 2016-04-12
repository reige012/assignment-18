#!/usr/bin/env python
# encoding: utf-8

"""
Assignment 18
Task 1 Program: Using the NCBI Taxonomy Database
Jessie Salter
9 April 2016
"""

import argparse
import re
from Bio import Entrez
import time


def parser():
    '''Takes user input'''
    parser = argparse.ArgumentParser(
        description='''Gets the path to the input file''')
    parser.add_argument(
        "--input",
        required=True,
        help="Enter the path to the input file of species names",
        type=str
        )
    parser.add_argument(
        "--output",
        required=True,
        help="Enter the name of the output file you want to create as a .csv",
        type=str
        )
    return parser.parse_args()


def reader(input_file):
    '''Reads in the list of species names from the input file'''
    with open(input_file, 'r') as infile:
        species = infile.read()
        spp_list = re.split('\n', species)
        spp_list = list(filter(None, spp_list))
    return spp_list


def esearch(spp_list):
    '''Collects the Taxonomy Identifiers from NCBI'''
    id_list = []
    for entry in spp_list:
        query = Entrez.esearch(db="taxonomy", term=entry, retmode="xml")
        result = Entrez.read(query)
        iden = result['IdList']
        iden = "txid{}".format(iden[0])
        id_list.append(iden)
        time.sleep(1)
    return id_list


def efetch(id_list):
    '''Collects the taxonomic info for each species in the list from the NCBI
    Taxonomy database'''
    results = []
    for iden in id_list:
        query = Entrez.efetch(db="taxonomy", id=iden, retmode="xml")
        data = Entrez.read(query)
        results.append(data)
        time.sleep(1)
    return results


def lineageex(results):
    '''Creates a new list with the LineageEx fields excised from the total
    search results'''
    linex = []
    for x in range(len(results)):
        # This includes a list of the full taxonomy and full scientific name:
        linex.append([results[x][0]['LineageEx'], results[x][0]['ScientificName']])
    return linex


def taxonomy_dict(linex, taxonomy):
    '''Creates a dictionary with each taxonomic level as a key and the
    scientific name as the corresponding value for each species'''
    tax_dict = dict()
    all_spp = []
    # Each entry is a new species:
    for entry in linex:
        # Each record is a different rank dictionary:
        for record in entry[0]:
            for key, value in record.items():
                # This iterates over every taxonomic level:
                for level in taxonomy:
                    if key == 'Rank' and value == level:
                        tax_dict[level] = record['ScientificName']
                    else:
                        pass
            # This adds the full scientific name:
            tax_dict['scientific name'] = entry[1]
        # dict.copy() will keep the results from overwriting:
        all_spp.append(dict.copy(tax_dict))
    return all_spp


def results_writer(output, taxonomy, all_spp):
    '''Writes the results of the query to an output CSV file '''
    with open(output, 'w') as outfile:
        # This writes the header:
        for level in taxonomy:
            outfile.write('{},'.format(level))
        # This puts a line break before the species data start:
        outfile.write('\n')
        # This writes out the corresponding value or blank for each column:
        for species in all_spp:
            for level in taxonomy:
                if level in species:
                    outfile.write('{},'.format(species[level]))
                else:
                    outfile.write(',')
            # This starts a new line before the next species:
            outfile.write('\n')


def main():
    Entrez.email = "jsalte5@lsu.edu"
    taxonomy = ['scientific name', 'superclass', 'class', 'subclass', 'infraclass', 'superorder', 'order', 'superfamily', 'family', 'genus']
    args = parser()
    spp_list = reader(args.input)
    id_list = esearch(spp_list)
    results = efetch(id_list)
    linex = lineageex(results)
    all_spp = taxonomy_dict(linex, taxonomy)
    results_writer(args.output, taxonomy, all_spp)


if __name__ == '__main__':
    main()
