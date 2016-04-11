#!/usr/bin/env python
# encoding: utf-8

"""
Script takes a species name as input and returns all sequence data related to
that species in fasta format to an output file specified by the user.

Edited by Alicia Reigel. 9 April 2016.
Copyright Alicia Reigel. Louisiana State University. 9 April 2016. All
rights reserved.

"""


import time
import argparse
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parser_get_args():
    """Collect the path to the file of species names and name of output file"""
    parser = argparse.ArgumentParser(
        description="""Input the species name and desired output file name"""
        )
    parser.add_argument(
            '--speciesname',
            required=True,
            type=str,
            help='Enter the species name.'
        )
    parser.add_argument(
            '--outputfile',
            required=True,
            type=str,
            help='Enter the desired name for the output fasta file make sure it ends in .faa'
        )
    return parser.parse_args()


def collect_sequence_data(species_name, seq_id_list, outputfile):
    for record in seq_id_list:
        time.sleep(1)
        genbank_seq = Entrez.efetch(db="nucleotide", id=record,  rettype="gb", retmode="text")
        record = SeqIO.read(genbank_seq, 'genbank')
        new_record = SeqRecord(seq=record.seq, id=record.id, name=record.name, description=species_name)
        with open(outputfile, 'a') as output:
            output.write(new_record.format('fasta'))


def main():
    args = parser_get_args()
    outputfile = args.outputfile
    species = str(args.speciesname)
    Entrez.email = "areige1@lsu.edu"
    search_query = Entrez.esearch(db="taxonomy", term=species, retmode="xml")
    result = Entrez.read(search_query)
    species_id = result['IdList']
    species_id = "txid{}".format(species_id[0])
    genbank_data = Entrez.esearch(db="nucleotide", term=species_id, retmode="xml")
    data = Entrez.read(genbank_data)
    data_seq_ids = data['IdList']
    collect_sequence_data(species, data_seq_ids, outputfile)


if __name__ == '__main__':
    main()
