#! /usr/bin/env python
# encoding: utf-8

'''
Grace Cagle
Assignment 18 Task 1

A program reading a file of species names, retriving taxonomic information
from NCBI and writing it to a CSV file
'''

import argparse
# import glob
# import os
from Bio import Entrez
import time
import csv

names = ['superclass', 'class', 'subclass', 'infraclass', 'superorder',
'order', 'superfamily', 'family', 'genus']


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', "--infile",
        required=True,
        help="A file containing species names"
    )
    parser.add_argument(
        '-o', "--outfile",
        required=True,
        help="The output file name"
    )
    return parser.parse_args()


def read_file(args):
    with open(args.infile, 'r') as f:
        sp = f.read().strip('\n').split('\n')
        # print(sp)
        return sp


def get_ID(sp):
    Entrez.email = 'gcagle1@lsu.edu'
    tax_id_list = []
    for s in sp:
        esearch_query = Entrez.esearch(db='taxonomy', term=s, retmode='xml')
        esearch_result = Entrez.read(esearch_query)
        tax_id = esearch_result['IdList']
        tax_id_list.append(tax_id)
        time.sleep(1)
    return tax_id_list


def get_tax(tax_id_list):
    taxonomy = []
    for each_id in tax_id_list:
        search = Entrez.efetch(id=each_id, db="taxonomy", retmode='xml')
        data = Entrez.read(search)
        taxonomy.append(data)
    return taxonomy


def make_dict(tax_id_list, taxonomy):
    data = []
    for x in range(len(tax_id_list)):
        data.append({d['Rank']: d['ScientificName'] for d in taxonomy[x][0]['LineageEx'] if d['Rank'] in names})
    return data


def write_csv(args, data):
    with open(args.outfile, 'w') as csvfile:
        writer = csv.DictWriter(csvfile, restval='', fieldnames=names)
        writer.writeheader()
        for d in data:
            writer.writerow(d)


def main():
    args = get_args()
    sp = read_file(args)
    tax_id_list = get_ID(sp)
    taxonomy = get_tax(tax_id_list)
    data = make_dict(tax_id_list, taxonomy)
    write_csv(args, data)


if __name__ == '__main__':
    main()
