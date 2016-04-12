#! /usr/bin/env python
# encoding UTF-8

'''
Assignment18Task1 biol7800
ZacCarver 04/12/2016
'''
import argparse
import csv
import time
from Bio import Entrez


def args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tax_list", type=str,
                        help="provide taxa list as text file")
    parser.add_argument("--tax_csv", )
    return parser.parse_args()


def Entrez_query(tax_list):
    #taxonomy = []
    with open(tax_list, 'r') as tl:
        for spp in iter(tl):
            #print(spp)
            time.sleep(3)
            e_query = Entrez.esearch(db="taxonomy",
                                     term=spp,
                                     retmode="xml")
            e_result = Entrez.read(e_query)
            #print(e_result["IdList"][0])
            ids = e_result["IdList"][0]
            #no for loop?
            e_fetched = Entrez.efetch(db="taxonomy",
                                      id=ids,
                                      retmode="xml")
            data = Entrez.read(e_fetched)
            #print(data[0]["Lineage"])
            '''stackoverflow.com/questions/16504238/attempting-to-obtain
            -taxonomic-information-from-biopython'''
            #this was rather tricky I thought...
            for spp in data:
                lineage = []
                l = {d['Rank']: d['ScientificName']
                     for d in data[0]['LineageEx']
                     if d['Rank'] in ['superclass', 'class', 'subclass',
                                      'infraclass', 'superorder', 'order',
                                      'superfamily', 'family', 'genus']}
                lineage.append(l)
                return lineage


def tocsv(l, tax_csv):
    '''stackoverflow.com/questions/28277150/write-a-list-in-a-\
    python-csv-file-one-new-row-per-list'''
    #the ugliest csv file ever created...
    with open(tax_csv, 'a') as f:
        out = csv.writer(f)
        out.writerow(l)


def main():
    Entrez.email = "zcarve1@lsu.edu"
    arg = args()
    lineage = Entrez_query(arg.tax_list)
    tocsv(lineage, arg.tax_csv)
    #keys = lineage[0].keys()
    '''print(lineage)
    with open(arg.tax_csv, 'w') as f:
        keynames = ['superclass', 'class', 'subclass',
                    'infraclass', 'superorder', 'order',
                    'superfamily', 'family', 'genus']
        writer = csv.DictWriter(f, fieldnames=keynames)
        writer.writeheader()
        writer.writerows(lineage)'''

if __name__ == '__main__':
    main()
