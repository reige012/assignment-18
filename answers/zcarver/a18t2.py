#! /usr/bin/env python
# encoding UTF-8

'''
Assignment18Task2 biol7800
ZacCarver 04/12/2016
Thanks given to JNations for reminding me of my "p's and 'q's", by setting 'id'
for efetch instead of 'term'...
'''
import argparse
from Bio import Entrez
#from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord


def args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--spp", type=str,
                        help="provide, in this format, genus+species")
    parser.add_argument("--fasta",
                        help="provide name-for-file + .fasta")
    return parser.parse_args()


def Entrez_query(spp):
    #time limit here
    e_query = Entrez.esearch(db="taxonomy",
                             term=spp,
                             retmode="xml")
    e_result = Entrez.read(e_query)
    spp_id = e_result['IdList']
    #print(spp_id)
    spp_id = "txid{}".format(spp_id[0])
    gb_entries = Entrez.esearch(db="nucleotide",
                                term=spp_id,
                                retmode='xml')
    sp_data = Entrez.read(gb_entries)
    seq_data = sp_data['IdList']
    return seq_data


def Entrez_gb_search(seq_data):
    #biopython-cn.readthedocs.org/en/latest/en/chr09.html
    '''gb_entries = Entrez.esearch(db="nucleotide",
                                id=sp_id,
                                retmode='xml')
    e_result = Entrez.read(gb_entries)'''
    #print(e_result)
    #gb_id = e_result['IdList']
    '''for record in gb_id['IdList']:
        gb_ent = Entrez.esearch(db="nucleotide",
                                term=record,
                                retmode="xml")
    gb_data = Entrez.read(gb_ent)
    seq_id = gb_data['IdList'][0]'''
    for record in seq_data:
        gb = Entrez.efetch(db="nucleotide",
                           id=record,
                           rettype='gb',
                           retmode='text')
        rec = SeqIO.read(gb, 'genbank')
    #print(rec)
    return rec


def main():
    Entrez.email = "zcarve1@lsu.edu"
    arg = args()
    seq_data = Entrez_query(arg.spp)
    record = Entrez_gb_search(seq_data)
    #print(record)
    seq_record = SeqRecord.SeqRecord(record.seq,
                                     id=record.id,
                                     name=record.name,
                                     description=record.description
                                     )
    with open(arg.fasta, 'w') as o:
        SeqIO.write([seq_record], o, 'fasta')

if __name__ == '__main__':
    main()
