#!/usr/bin/env python
# encoding: utf-8

"""
Script takes a species scientific name as input (make sure name is in
quotatations) and the name of an output file.  It then finds the sequence
matches in the nucleotide database to that species and performs a blastn search
against each sequence match to return the hits. Hits are output to the file
name given by the user.

Note: This program works, but it takes at least 30 minutes to run. Please don't
shut it off without checking the output to make sure its correct.

Edited by Alicia Reigel. 11 April 2016.
Copyright Alicia Reigel. Louisiana State University. 11 April 2016. All
rights reserved.

"""

import os

import time
import argparse
from Bio import Entrez
from Bio import SeqIO
from Bio.Blast import NCBIWWW
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
            help='Enter the species name. Make sure it is in quotations.'
        )
    parser.add_argument(
            '--outputfile',
            required=True,
            type=str,
            help='Enter the desired name for the output file'
        )
    return parser.parse_args()


def collect_seqs(species_name, seq_id_list, outputfile_in):
    outputfile = outputfile_in
    for record in seq_id_list:
        time.sleep(1)
        genbank_seq = Entrez.efetch(db="nucleotide", id=record,  rettype="gb", retmode="text")
        record = SeqIO.read(genbank_seq, 'genbank')
        new_record = SeqRecord(seq=record.seq, id=record.id, name=record.name, description=species_name)
        blast_record(new_record.format('fasta'), outputfile)


def blast_record(fastasequence, output_file):
    outputfile = os.path.join(output_file + '.xml')
    result = NCBIWWW.qblast("blastn", "nt", fastasequence)
    save_file = open(outputfile, "a")
    save_file.write(result.read())
    save_file.close()
    result.close()


def main():
    args = parser_get_args()
    species = str(args.speciesname)
    Entrez.email = "areige1@lsu.edu"
    search_query = Entrez.esearch(db="taxonomy", term=species, retmode="xml")
    result = Entrez.read(search_query)
    species_id = result['IdList']
    species_id = "txid{}".format(species_id[0])
    genbank_data = Entrez.esearch(db="nucleotide", term=species_id, retmode="xml")
    data = Entrez.read(genbank_data)
    data_seq_ids = data['IdList']
    collect_seqs(species, data_seq_ids, args.outputfile)


if __name__ == '__main__':
    main()
