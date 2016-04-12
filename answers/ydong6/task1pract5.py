#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on Apr 11, 2016
A18T1
There are two records are out of date, so delete the:"rhynchospiza strigiceps" and"peucaea ruficauda" to make script to work.
@author: York
'''

from Bio import Entrez
import argparse
import time


Entrez.email = "catchdonovan@gmail.com"


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputfile", required=True)
    parser.add_argument("--outputfile", required=True)
    args = parser.parse_args()
    return args


def spe_id(species):
    species = species.replace(" ", "+").strip()
    search = Entrez.esearch(term=species, db="taxonomy", retmode="xml")
    record = Entrez.read(search)
    return record['IdList'][0]


def fetch(idid):
    search = Entrez.efetch(id=idid, db="taxonomy", retmode="xml")
    # print(search)
    return Entrez.read(search)


def main():
    args = get_parser()
    outfile = open(args.outputfile, "w+")
    with open(args.inputfile, "r") as infile:
        for item in infile:
            time.sleep(1)
            idid = spe_id(item)
            outfile.write(item)
            efetch = fetch(idid)
            csv_result = {d['Rank']: d['ScientificName']
                          for d in efetch[0]['LineageEx']
                          if d['Rank'] in ['superclass', 'class', 'subclass', 'infraclass', 'superorder', 'family', 'order', 'superfamily', 'family', 'genus']}
            outfile.write(str(csv_result))
            print(csv_result)
            if item == 0:
                break


if __name__ == '__main__':
    main()
