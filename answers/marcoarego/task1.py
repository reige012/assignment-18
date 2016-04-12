# !/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Marco
# @Date:   2016-04-07 13:51:09
# @Last Modified by:   Marco
# @Last Modified time: 2016-04-12 09:35:41

import argparse
import re
import time
from Bio import Entrez


def parser_function():
    '''
    Function to parse arguments
    '''
    parser = argparse.ArgumentParser(description='''This program gets a
         text file with species names as input to create a taxonomic query on
         NCBI and retrieve taxonomic levels that have a defined Rank on NCBI
        ''')
    # Adding an argument to 'parse'
    parser.add_argument('in_file', type=str,
                        help='type the name of the input txt file')
    parser.add_argument('out_file', type=str,
                        help='type the name of the output csv file')
    parser.add_argument('email', type=str,
                        help='type your email address for Entrez.email')
    args = parser.parse_args()
    inp = args.in_file
    outp = args.out_file
    email = args.email
    return inp, outp, email


def get_list(file):
    '''oppening the list file with the species name'''
    with open(file) as f:
        return [line.strip("\n") for line in f]


def idlist_query(email, list_of_taxa):
    '''this function creates a query to get the Id number for each species'''
    Entrez.email = email
    idlist = []
    for taxon in list_of_taxa:
        try:
            my_query = Entrez.esearch(db="taxonomy", term=taxon, retmode="xml")
            idlist.append(Entrez.read(my_query)['IdList'][0])
            time.sleep(0.5)
        except:
            print("ATTENTION!!\n{} did not have an Id number in the"
                  " IdList field -- Please check spelling!"
                  .format(taxon.capitalize()))
    return idlist


def tax_query(email, id_list):
    '''
    This function gets any taxonomic category that has a stabilished rank
    for each taxa
    '''
    Entrez.email = email
    lis = []
    for ident in id_list:
        known_ranks = []
        my_query = Entrez.efetch(db="taxonomy", id=ident, retmode="xml")
        file_read = Entrez.read(my_query)
        sci_name = file_read[-1]['ScientificName']
        for x in range(0, len(file_read[0]['LineageEx'])):
            if (file_read[0]['LineageEx'][x]['Rank']) != 'no rank':
                known_ranks.append(str(file_read[0]['LineageEx'][x]['Rank']) +
                                   ": " + str(file_read[0]['LineageEx'][x]
                                              ['ScientificName']))
        # Next piece of code was adapted from:
        # http://www.cademuir.eu/blog/2011/10/20/python-searching-for-a-string-within-a-list-list-comprehension/
        regex = re.compile(".*(^family.*$).*")
        regex2 = re.compile(".*(^superfamily.*$).*")
        regex3 = re.compile(".*(^subfamily.*$).*")
        index_family = [known_ranks.index(a.group(0)) for elem1 in known_ranks
                        for a in [regex.search(elem1)] if a][0]
        index_superf = [b.group(0) for elem2 in known_ranks for b in
                        [regex2.search(elem2)] if b]
        index_subf = [c.group(0) for elem3 in known_ranks for c in
                      [regex3.search(elem3)] if c]
        if index_subf == []:
            known_ranks.insert(index_family+1, "subfamily: none")
        if index_superf == []:
            known_ranks.insert(index_family, "superfamily: none")
        # End of adapted code #
        rank_join = ";".join(known_ranks)
        lis.append("{}; {}".format(str(sci_name), rank_join))
        time.sleep(0.5)
    return lis


def tax_query_to_csv(my_query, csv_query_file):
    '''this function will get the idlist_query output, which is a list
    of strings containing the taxonomy of desired species, separeted by
    semi-colons'''
    splitter = [line.split(";") for line in my_query]
    with open(csv_query_file, 'w') as csv_f:
        for item in splitter:
            csv_f.write(",".join(item))
            csv_f.write("\n")


def main():
    inp, outp, email = parser_function()
    lis = get_list(inp)
    id_list = idlist_query(email, lis)
    final_lis = tax_query(email, id_list)
    tax_query_to_csv(final_lis, outp)


if __name__ == '__main__':
    main()
