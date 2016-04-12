# !/usr/bin/env python
# encoding: utf-8

"""
Assignment 18 task 1 program: Get taxonomic info from NCBI
Created by Andre Moncrieff on 11 April 2016.
Copyright 2016 Andre E. Moncrieff. All rights reserved.

This assignment asks for output values for superclass, class, subclass,
infraclass, superorder, order, superfamily, family, and genus. I check for the
availability of information for these ranks in NCBI and only return what is
available.

"""

import argparse
import csv
import time
from Bio import Entrez


Entrez.email = "amoncr2@lsu.edu"


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_filename", required=True,
                        help="Enter full filename containing species list.",
                        type=str)
    parser.add_argument("--out_filename", required=True,
                        help="Enter full filename for output file.",
                        type=str)
    args = parser.parse_args()
    return args


def read_in_sp_list(in_filename):
    species_list = []
    with open(in_filename, "r") as in_file:
        for line in in_file:
            species_list.append(line.replace("\n", ""))
    return species_list


def get_sp_identifiers(species_list):
    "Heavily based on Brant's Lecture 18 code"
    species_identifiers = []
    for species in species_list:
        esearch_query = Entrez.esearch(db="taxonomy",
                                       term=species,
                                       retmode="xml")
        esearch_result = Entrez.read(esearch_query)
        identifier = esearch_result.get("IdList", "No available data")[0]
        species_identifiers.append(identifier)
        time.sleep(1)
    return species_identifiers


def raw_taxonomies(species_identifiers):
    nested_list_of_raw_taxonomies = []
    for sp_id in species_identifiers:
        taxonomy_entry = Entrez.efetch(db="taxonomy",
                                       id=sp_id,
                                       retmode="xml")
        efetch_result = Entrez.read(taxonomy_entry)
        raw_taxonomies = list(efetch_result[0].get("LineageEx"))
        nested_list_of_raw_taxonomies.append(raw_taxonomies)
    return nested_list_of_raw_taxonomies


def mk_rank_name_dict(nested_list_of_raw_taxonomies):
    list_of_dicts = []
    for one_list in nested_list_of_raw_taxonomies:
        rank_name_dict = {}
        for dictionary in one_list:
            ScientificName = dictionary["ScientificName"]
            Rank = dictionary["Rank"]
            rank_name_dict[Rank] = ScientificName
        list_of_dicts.append(rank_name_dict)
    return list_of_dicts


def parse_dict(list_of_rank_name_dicts, species_list):
    wanted_ranks = ["superclass", "class", "subclass", "infraclass",
                    "superorder", "order", "superfamily", "family", "genus"]
    nested_list_of_species_taxonomies = []
    counter = 0
    for dictionary in list_of_rank_name_dicts:
        nested_list_of_species_taxonomies.append([species_list[counter]])
        for rank in wanted_ranks:
            if rank in dictionary.keys():
                name = dictionary[rank]
                a_string = "{}: {}".format(rank, name)
                nested_list_of_species_taxonomies[counter].extend([a_string])
        counter += 1
    return nested_list_of_species_taxonomies


def cvs_writer(nested_list_of_species_taxonomies, out_filename):
    writer = csv.writer(open(out_filename, "w"))
    for sp_tax_list in nested_list_of_species_taxonomies:
        writer.writerow(sp_tax_list)


def main():
    args = parser()
    species_list = read_in_sp_list(args.in_filename)
    # print(species_list)
    species_identifiers = get_sp_identifiers(species_list)
    # print(species_identifiers)
    nested_list_of_raw_taxonomies = raw_taxonomies(species_identifiers)
    # print(nested_list_of_raw_taxonomies)
    list_of_rank_name_dicts = mk_rank_name_dict(nested_list_of_raw_taxonomies)
    # print(list_of_rank_name_dicts)
    nested_list_of_species_taxonomies = parse_dict(list_of_rank_name_dicts,
                                                   species_list)
    cvs_writer(nested_list_of_species_taxonomies, args.out_filename)


if __name__ == '__main__':
    main()
