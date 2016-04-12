#!/usr/bin/env python
# encoding: utf-8

"""
BIOL7800 Assignment 18 Task 1

Amie Settlecowski
12 Apr. 2016

This program searches the ncbi taxonomy database for a list of scientific names
input as a text file, and output their taxonomy for a defined rank_list as a
csv file.
default rank_list:  'superclass', 'class', 'subclass', 'infraclass',
'superorder', 'order', 'superfamily', 'family', 'genus'

python ../task1_settlecowski.py --i ../pathto/input.txt --o output.csv
"""
import argparse
import os
import time
from Bio import Entrez


def get_args(parserr):
    '''
    Requires --i flag for user to specify path to input file
    Requires --o flag for user to name output file
    '''
    parserr.add_argument("--i",
        required=True,
        help="Path to input file",
        type=str)

    parserr.add_argument("--o",
        required=True,
        help="Name output file",
        type=str)


def query_search(in_file):
    '''Generate list of ncbi IDs for each species'''
    esearch_result = []
    for line in in_file:
        taxon = line
        esearch_query = Entrez.esearch(db="taxonomy",
                                        term=taxon,
                                        retmode="xml")
        esearch_result.append(Entrez.read(esearch_query)['IdList'][0])
    return esearch_result


def get_taxonomy(taxa_ID, rank_list):
    '''Mine ncbi taxonomy by species ID'''
    tax_entry = Entrez.efetch(db="taxonomy", id=taxa_ID, retmode="xml")
    result = Entrez.read(tax_entry)
    name = result[0]['ScientificName']
    taxonomy = {}
    for index in result[0]['LineageEx']:
        if index['Rank'] in rank_list:
            taxonomy[index['Rank']] = index['ScientificName']
    return name, taxonomy


def write_taxonomy(name, tax, rank_list, out_file):
    '''Writes one line of csv file for each species'''
    out_file.write('\n{},'.format(name))
    for rank in rank_list:
        if rank in tax.keys():
            out_file.write('{},'.format(tax[rank]))
        else:
            out_file.write(',')


def main():
    # Create Parser with arguments for input and output files
    file_name_parser = argparse.ArgumentParser()
    get_args(file_name_parser)
    file_args = file_name_parser.parse_args()

    # Change to directory with input file
    working_dir = os.path.split(file_args.i)[0]
    os.chdir(working_dir)

    # set email for Entrez
    Entrez.email = 'asettl1@lsu.edu'

    with open(file_args.i, 'r') as in_file:
        with open(file_args.o, 'w') as out_file:
            # define list of ranks to be mined and written out to file
            rank_list = ['scientific name', 'superclass', 'class', 'subclass', 'infraclass', 'superorder', 'order', 'superfamily', 'family', 'genus']
            # write header row of csv file
            for rank in rank_list:
                out_file.write('{},'.format(rank))
            rank_list.remove('scientific name')
            # generate list of species IDs from ncbi
            taxa_IDs = query_search(in_file)
            # for each species, mine taxonomy from ncbi and write out to file
            for ID in taxa_IDs:
                sci_name, taxonomy = get_taxonomy(ID, rank_list)
                write_taxonomy(sci_name, taxonomy, rank_list, out_file)
                time.sleep(1)


if __name__ == '__main__':
    main()
