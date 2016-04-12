# !/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Marco
# @Date:   2016-04-10 08:16:09
# @Last Modified by:   Marco
# @Last Modified time: 2016-04-10 16:04:11


# import pdb
import time
import argparse

from Bio import Seq
from Bio import SeqIO
from Bio import Entrez
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


def parser_function():
    '''
    Function to parse arguments
    '''
    parser = argparse.ArgumentParser(description='''This program looks for
        sequence data in NCBI for a particular species and write them in a
        fasta file. The input shell arguments should be like this:
        python task2.py Loxosceles_reclusa my_file.fasta''')
    # Adding an argument to 'parse'
    parser.add_argument('in_taxon', type=str,
                        help='type the name of the species that you are'
                             ' querying. Please use underlines (_) between'
                             ' genus and species names')
    parser.add_argument('out_file', type=str,
                        help='type the name of the output .fasta file.')
    args = parser.parse_args()
    inp = args.in_taxon
    outp = args.out_file
    return inp, outp


def idlist_query(email, list_of_taxa):
    '''this function creates a query to get the Id number for each species'''
    Entrez.email = email
    idlist = []
    if type(list_of_taxa) == list:
        for taxon in list_of_taxa:
            try:
                my_query = Entrez.esearch(db="taxonomy", term=taxon,
                                          retmode="xml")
                idlist.append(Entrez.read(my_query)['IdList'][0])
                time.sleep(0.5)
            except:
                print("ATTENTION!!\n{} did not have an Id number in the"
                      " IdList field -- Please check spelling!"
                      .format(taxon.capitalize()))
        return idlist
    elif type(list_of_taxa) == str:
        try:
            my_query = Entrez.esearch(db="taxonomy", term=list_of_taxa,
                                      retmode="xml")
            idlist.append(Entrez.read(my_query)['IdList'][0])
        except:
            print("ATTENTION!!\n{} did not have an Id number in the"
                  " IdList field -- Please check spelling!"
                  .format(taxon.capitalize()))
        return idlist


def nucleotide_query(email, id_list):
    '''
    This function gets any taxonomic category that has a stabilished rank
    for each taxa
    '''
    Entrez.email = email
    seq_list = []
    for ident in id_list:
        entry = "txid"+str(ident)
        my_query = Entrez.esearch(db="nucleotide", term=entry, retmode="xml")
        file_read = Entrez.read(my_query)
        for record in file_read['IdList']:
            genbank_entry = Entrez.efetch(db="nucleotide", id=record,
                                          retmode="xml")
            data = Entrez.read(genbank_entry)
            seq_list.append(data)
            time.sleep(1)
    return seq_list


def write_to_fasta(list_of_seqs, output_name):
    '''
    This function creates Seq records and writes them to a fasta file
    '''
    with open(output_name, 'w') as outfile:
        for record in list_of_seqs:
            try:
                identifier = "".join(record[0]['GBSeq_other-seqids'])
                descr = (record[0]['GBSeq_source'] + '. Length:' +
                         record[0]['GBSeq_length'])
                string_sequence = (record[0]['GBSeq_sequence']).upper()
                seq_object = Seq.Seq(string_sequence, IUPAC.ambiguous_dna)
                new = SeqRecord(seq_object, id=identifier, name=identifier,
                                description=descr)
                SeqIO.write([new], outfile, 'fasta')
            except:
                pass


def main():
    inp, outp = parser_function()
    my_input = inp.replace("_", " ")
    idlist = idlist_query('mrego1@lsu.edu', my_input)
    nuc_query = nucleotide_query('mrego1@lsu.edu', idlist)
    # pdb.set_trace()
    write_to_fasta(nuc_query, outp)

if __name__ == '__main__':
    main()
