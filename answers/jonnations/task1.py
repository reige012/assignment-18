#!/usr/bin/env python
# utf-8

"""
Assignment 18, Task 1
Jon Nations
Write a general program that takes this list as input (using argparse), gets the taxonomic information from NCBI for each sample, and writes that to a new file as output (again, get the output file name using argparse). Format the output file as a CSV. The output file should contain the values for each species for the following: superclass, class, subclass, infraclass, superorder, order, superfamily, family, and genus. Be sure that your program is formatted correctly (PEP8) and be sure that your file works appropriately with any generic list of species. Also be sure that you pass your email to NCBI (Entrez.email) and that you limit your requests to 1 per second.
"""
import argparse
from Bio import Entrez
import time
# import pdb


def file_in():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', required=True, type=str,
                        help="give input file name.")
    parser.add_argument('-o', '--output', required=True, help="give output file   (.csv format) where taxonomy list is going.")
    return parser.parse_args()


def search(args):
    with open(args.infile, 'r') as infile:
        taxlist = infile.read().splitlines()
    for taxon in taxlist:
        Entrez.email = 'jonnations@gmail.com'
        time.sleep(1)
        search_query = Entrez.esearch(db="taxonomy", term=taxon, retmode="xml")
        s_result = Entrez.read(search_query)
        # print(s_result)
        result = s_result['IdList']
        # print(result)
        # result = 'taxonid{}'.format(result[0])
        s_fetch = Entrez.efetch(db='taxonomy', id=result, retmode='xml')
        f_result = Entrez.read(s_fetch)
        lin_ex = f_result[0]["LineageEx"]
        # print(lin_ex) This prints correctly
        return lin_ex


def superclass(lin_ex):
    dict1 = {}
    for level in lin_ex:
        for k, v in level.items():
            if k == ('Rank') and v == ('superclass'):
                dict1['superclass'] = level['ScientificName']
            else:
                pass
    return dict1


def reg_class(lin_ex):
    dict2 = {}
    for level in lin_ex:
        for k, v in level.items():
            if k == ('Rank') and v == ('class'):
                dict2['class'] = level['ScientificName']
            else:
                pass
    return dict2


def subclass(lin_ex):
    dict3 = {}
    for level in lin_ex:
        for k, v in level.items():
            if k == ('Rank') and v == ('subclass'):
                dict3['subclass'] = level['ScientificName']
            else:
                pass
    return dict3


def infraclass(lin_ex):
    dict4 = {}
    for level in lin_ex:
        for k, v in level.items():
            if k == ('Rank') and v == ('infraclass'):
                dict4['infraclass'] = level['ScientificName']
            else:
                pass
    return dict4


def superorder(lin_ex):
    dict5 = {}
    for level in lin_ex:
        for k, v in level.items():
            if k == ('Rank') and v == ('superorder'):
                dict5['superorder'] = level['ScientificName']
            else:
                pass
    return dict5


def order(lin_ex):
    dict6 = {}
    for level in lin_ex:
        for k, v in level.items():
            if k == ('Rank') and v == ('order'):
                dict6['order'] = level['ScientificName']
            else:
                pass
    return dict6


def superfamily(lin_ex):
    dict7 = {}
    for level in lin_ex:
        for k, v in level.items():
            if k == ('Rank') and v == ('superfamily'):
                dict7['superfamily'] = level['ScientificName']
            else:
                pass
    return dict7


def family(lin_ex):
    dict8 = {}
    for level in lin_ex:
        for k, v in level.items():
            if k == ('Rank') and v == ('family'):
                dict8['family'] = level['ScientificName']
            else:
                pass
    return dict8


def genus(lin_ex):
    dict9 = {}
    for level in lin_ex:
        for k, v in level.items():
            if k == ('Rank') and v == ('genus'):
                dict9['genus'] = level['ScientificName']
            else:
                pass
    return dict9


def write_tax(args, dict1, dict2, dict3, dict4, dict5, dict6, dict7, dict8, dict9):
    with open(args.output, 'a') as output:
        output.write(str('dict1', 'dict2', 'dict3', 'dict4', 'dict5', 'dict6', 'dict7', 'dict8', 'dict9'))
        output.write('\n\n')

# taxlist2 = [(d['Rank'],d['ScientificName']) for d in taxlist]


def main():
    args = file_in()
    file_in()
    lin_ex = search(args)
    # lin_ex = f_result[0]['LineageEx']
    dict1 = superclass(lin_ex)
    dict2 = reg_class(lin_ex)
    dict3 = subclass(lin_ex)
    dict4 = infraclass(lin_ex)
    dict5 = superorder(lin_ex)
    dict6 = order(lin_ex)
    dict7 = superfamily(lin_ex)
    dict8 = family(lin_ex)
    dict9 = genus(lin_ex)
    write_tax(args, dict1, dict2, dict3, dict4, dict5, dict6, dict7, dict8,
              dict9)


if __name__ == '__main__':
    main()
