# !/usr/bin/env python
# encoding: utf-8

"""
Assignment 18 task 2 program:
Created by Andre Moncrieff on 11 April 2016.
Copyright 2016 Andre E. Moncrieff. All rights reserved.



"""

import argparse
import time
from Bio import Entrez
from Bio import SeqIO


Entrez.email = "amoncr2@lsu.edu"


def parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--name", required=True,
                        help="Enter scientific name of the organism of" +
                        " in format Genus_species",
                        type=str)
    parser.add_argument("--out_filename", required=True,
                        help="Enter filename (without extension) for fasta" +
                        " output file",
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


def get_genbank_data(genbank_identifiers, out_filename):
    counter = 1
    for record in genbank_identifiers:
        genbank_record = Entrez.efetch(db="nucleotide",
                                       id=record,
                                       rettype="gb",
                                       retmode="text")
        record = SeqIO.read(genbank_record, "genbank")
        with open("{}_{}.fasta".format(out_filename, counter), "w") as outfile:
            outfile.write(record.format("fasta"))
        counter += 1
        time.sleep(1)


def main():
    args = parser()
    tax_identifier = get_sp_identifier(args.name)
    #print(tax_identifier)
    genbank_identifiers = get_genbank_identifiers(tax_identifier)
    #print(genbank_identifiers)
    get_genbank_data(genbank_identifiers, args.out_filename)


if __name__ == '__main__':
    main()
