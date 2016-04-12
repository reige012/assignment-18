# !/usr/bin/env python
# encoding: utf-8

"""
Assignment 18 task 3 program:
Created by Andre Moncrieff on 11 April 2016.
Copyright 2016 Andre E. Moncrieff. All rights reserved.



"""

import argparse
import time
import os
from Bio import Entrez
from Bio.Blast import NCBIWWW


Entrez.email = "amoncr2@lsu.edu"


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", required=True,
                        help="Enter scientific name of the organism of" +
                        " in format Genus_species",
                        type=str)
    parser.add_argument("--out_dir", required=True,
                        help="Enter desired full path of a new directory for" +
                        " storing output files",
                        type=str)
    args = parser.parse_args()
    return args


def get_sp_identifier(name):
    "Heavily based on Brant's Lecture 18 code"
    sp_name_0 = name.split("_")
    sp_name = "{} {}".format(sp_name_0[0], sp_name_0[1])
    esearch_query = Entrez.esearch(db="taxonomy",
                                   term=sp_name,
                                   retmode="xml")
    esearch_result = Entrez.read(esearch_query)
    tax_identifier = esearch_result.get("IdList", "No available data")[0]

    return tax_identifier


def get_genbank_identifiers(tax_identifier):
    genbank_entries = Entrez.esearch(db="nucleotide",
                                     term="txid"+tax_identifier,
                                     retmode="xml")
    esearch_result = Entrez.read(genbank_entries)
    genbank_identifiers = esearch_result["IdList"]
    return genbank_identifiers


def get_genbank_data(genbank_identifiers):
    nested_data_list = []
    for record in genbank_identifiers:
        genbank_record = Entrez.efetch(db="nucleotide",
                                       id=record,
                                       retmode="xml")
        data = Entrez.read(genbank_record)
        nested_data_list.append(data)
        time.sleep(1)
    return nested_data_list


def get_GI_seqIDs(nested_data_list):
    GI_seqIDs = []
    for data_block in nested_data_list:
        GIs = data_block[0]["GBSeq_other-seqids"][1]
        GIs_num_only = GIs.replace("gi|", "")
        GI_seqIDs.append(GIs_num_only)
    return GI_seqIDs


def blast_with_GIs(GI_seqIDs):
    """
    Based in part on Biopython cookbook example.
    Try and except structure of function suggested by Subir
    """
    counter = 1
    for GI_ID in GI_seqIDs:
        try:
            result = NCBIWWW.qblast("blastn", "nt", GI_ID, format_type="Text")
            blast_results = result.read()
            with open("{}_{}.txt".format("blast_results", counter), "w") as outfile:
                outfile.write(blast_results)
            counter += 1
        except:
            print("No sequence available for gi|{}".format(GI_ID))
        time.sleep(1)


def main():
    args = parser()
    tax_identifier = get_sp_identifier(args.name)
    # print(tax_identifier)
    genbank_identifiers = get_genbank_identifiers(tax_identifier)
    # print(genbank_identifiers)
    nested_data_list = get_genbank_data(genbank_identifiers)
    # print(nested_data_list)
    GI_seqIDs = get_GI_seqIDs(nested_data_list)
    # print(GI_seqIDs)
    #GI_seqIDs = ['239819098', '815826725', '49458051', '49458049', '223037237', '219398693', '215260324', '57792506', '34484220', '225031058']
    os.makedirs(args.out_dir)
    os.chdir(args.out_dir)
    blast_with_GIs(GI_seqIDs)


if __name__ == '__main__':
    main()
