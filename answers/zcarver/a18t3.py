#! /usr/bin/env python
# encoding UTF-8

'''
Assignment18Task3 biol7800
ZacCarver 04/12/2016
You get a very long and ugly file instead of nice wee ones, joy.
'''
import argparse
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
#from Bio import Seq
#from Bio import SeqIO
#from Bio import SeqRecord


def args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--spp", type=str,
                        help="provide, in this format, genus+species")
    parser.add_argument("--blast_file",
                        help=" no destination directory, instead, filename.txt")
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
    return spp_id


def Entrez_gb_entries(spp_id):
    gb_entries = Entrez.esearch(db="nucleotide",
                                term=spp_id,
                                retmode='xml')
    sp_data = Entrez.read(gb_entries)
    seq_data = sp_data['IdList']
    return seq_data


def Entrez_gb_search(seq_data):
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
        #http://www.cbs.dtu.dk/courses/27624/IAH_2.pdf
        gb = Entrez.efetch(db="nucleotide",
                           id=record,
                           rettype='gb',
                           retmode='xml')
        rec = Entrez.read(gb)
    return rec[0]['GBSeq_other-seqids'][1]
    #return rec


def get_gi(rec):
    gi = rec.split('|')[1]
    #print(gi)
    b_result = NCBIWWW.qblast("blastn", "nt", gi)
    b_result = b_result.read()
    '''with open("qblast_result.xml", 'w') as b:
        b.write(b_result)
        with open("qblast_result.xml") as xml:
            b_record = NCBIXML.read(xml)
            return b_record'''


def main():
    Entrez.email = "zcarve1@lsu.edu"
    arg = args()
    spp_id = Entrez_query(arg.spp)
    seq_data = Entrez_gb_entries(spp_id)
    rec = Entrez_gb_search(seq_data)
    b_result = get_gi(rec)
    #print(b_result)
    with open(arg.blast_file, 'w') as b:
        b.write(str(b_result))
    '''blast_file = open(arg.blast_file, 'w')
    blast_file.write([b_result])
    blast_file.close()'''
    #print(record)
    """seq_record = SeqRecord.SeqRecord(record.seq,
                                     id=record.id,
                                     name=record.name,
                                     description=record.description
                                     )
    with open(arg.fasta, 'w') as o:
        SeqIO.write([seq_record], o, 'fasta')"""

if __name__ == '__main__':
    main()
