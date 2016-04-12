# !/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author: Marco
# @Date:   2016-04-10 16:09:29
# @Last Modified by:   Marco
# @Last Modified time: 2016-04-11 22:47:18


'''
Building on Task 2, modify your program above to get GI record of each sequence
from NCBI (e.g. do not work from the FASTA file you created) and the send
those records to NCBI BLAST to perform a blastn search of each sequence against
the nt (nucleotide) database. Write each BLAST result (the result for each
sequence, which may consist of several hits) out to text files in BLAST format
(use argparse to get the output directory name from the user). We didn't
cover BLAST searches in class, so you'll need to look at the BioPython Cookbook
and go from there. Be sure that your program is formatted correctly
(EP8). Also be sure that you pass your email to NCBI (Entrez.email) and that
you limit your requests to 1 per second.
'''

# import pdb
import os
import time
import argparse
from Bio import Entrez
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


def parser_function():
    '''
    Function to parse arguments
    '''
    parser = argparse.ArgumentParser(description=r'''This program looks for
        sequence data in NCBI for a particular species and write them in a
        fasta file. The input shell arguments should be like this:
        python task2.py Loxosceles_reclusa c:\users\xxxx\desktop''')
    # Adding an argument to 'parse'
    parser.add_argument('in_taxon', type=str,
                        help='type the name of the species that you are'
                             ' querying. Please use underlines (_) between'
                             ' genus and species names')
    parser.add_argument('out_file', type=str,
                        help='type the name of the output folder.')
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


def blaster_stuff(list_of_seqs, out_folder):
    for record in list_of_seqs:
        try:
            gi_number = (record[0]['GBSeq_other-seqids'][-1]).split("|")[-1]
            # pdb.set_trace()
            result = NCBIWWW.qblast("blastn", "nt", gi_number)
            blast_records = NCBIXML.parse(result)
            with open(os.path.join(out_folder, gi_number+".txt"), 'a') as txt:
                for rec in blast_records:
                    for alignment in rec.alignments:
                        for hsp in alignment.hsps:
                            txt.write('****Alignment****\n')
                            txt.write('sequence: ' + str(alignment.title))
                            txt.write('\nlength: ' + str(alignment.length))
                            txt.write('\ne value: ' + str(hsp.expect))
                            txt.write("\n"+hsp.query[0:75] + '...')
                            txt.write("\n"+hsp.match[0:75] + '...')
                            txt.write("\n"+hsp.sbjct[0:75] + '...\n')
        except:
            pass


def main():
    inp, output_folder = parser_function()
    my_input = inp.replace("_", " ")
    idlist = idlist_query('mrego1@lsu.edu', my_input)
    nuc_query = nucleotide_query('mrego1@lsu.edu', idlist)
    blaster_stuff(nuc_query, output_folder)

if __name__ == '__main__':
    main()
