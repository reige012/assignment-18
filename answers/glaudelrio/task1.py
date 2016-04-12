# !/usr/bin/env python
# -*- coding: utf-8 -*-
# Made with Marco Rego

import argparse
import re
import time

from Bio import Entrez


def arguments():
    '''
    Getting arguments
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('in_file', type=str,
                        help='input txt file name')
    parser.add_argument('out_file', type=str,
                        help='output csv file name')
    args = parser.parse_args()
    inp = args.in_file
    outp = args.out_file
    return inp, outp


def open_file(txt):
    '''a function that opens the list with species names'''
    with open(txt) as file:
        return [line.strip("\n") for line in file]


def taxonomic_information(email, species_list):
    '''Gets the NCBI information on the main taxonomic categories for a 
       series of species'''
    Entrez.email = email
    id_number_list = []
    final = []
    for taxon in species_list:
        try:
            species_information = Entrez.esearch(db="taxonomy", term=taxon, 
            	                                 retmode="xml")
            id_number_list.append(Entrez.read(species_information)['IdList'][0])
            time.sleep(0.5)
        except:
            pass
    for number in id_number_list:
        species_information = Entrez.efetch(db="taxonomy", id=number, 
        	                                retmode="xml")
        file = Entrez.read(species_information)
        sci_name = file[-1]['ScientificName']
        categories = []
        for item in range(0, len(file[0]['LineageEx'])):
            if (file[0]['LineageEx'][item]['Rank']) != 'no rank':
                categories.append(str(file[0]['LineageEx'][item]['Rank']) +
                                   ": " + str(file[0]['LineageEx'][item]
                                   ['ScientificName']))
        # Code adapted from http://www.cademuir.eu/blog/2011/10/20/python-
        #searching-for-a-string-within-a-list-list-comprehension:
        regex1 = re.compile(".*(^family.*$).*")
        regex2 = re.compile(".*(^superfamily.*$).*")
        regex3 = re.compile(".*(^subfamily.*$).*")
        index_family = [categories.index(family.group(0)) for elem1 in categories
                        for family in [regex1.search(elem1)] if family][0]
        index_superfamily = [superfamily.group(0) for elem2 in categories for 
                             superfamily in [regex2.search(elem2)] if superfamily]
        index_subfamily = [subfamily.group(0) for elem3 in categories for 
                           subfamily in [regex3.search(elem3)] if subfamily]
        if index_subfamily == []:
            categories.insert(index_family + 1, "subfamily: none")
        if index_superfamily == []:
            categories.insert(index_family, "superfamily: none")
        join_categories = ";".join(categories)
        final.append("{}; {}".format(str(sci_name), join_categories))
        time.sleep(0.5)
    return final


def taxonomic_categories(species_list, csv):
    '''creates the output and stores it in an csv file'''
    list1 = [line.split(";") for line in species_list]
    with open(csv, 'w') as csv_file:
        for x in list1:
            csv_file.write(",".join(x))
            csv_file.write("\n")


def main():
    inp, outp = arguments()
    species = open_file(inp)
    information = taxonomic_information('glaucia.ornito@gmail.com', species)
    taxonomic_categories(information, outp)


if __name__ == '__main__':
    main()
