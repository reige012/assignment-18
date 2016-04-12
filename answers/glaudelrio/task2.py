# !/usr/bin/env python
# -*- coding: utf-8 -*-
# Made with Marco Rego

import argparse
import time

from Bio import Seq
from Bio import SeqIO
from Bio import Entrez
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


def arguments():
    '''
    Gets species name as input and file name for final fasta
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('species', type=str,
                        help='species name separated by underlines')
    parser.add_argument('output', type=str,
                        help='type the name of the output .fasta file.')
    args = parser.parse_args()
    inp = args.species
    outp = args.output
    return inp, outp


def getting_sequences(email, taxon):
    '''Gets sequences from NCBI and appends them to a list'''
    taxon = taxon.replace("_", " ")
    Entrez.email = email
    id_numbers = []
    seqs = []
    species_information = Entrez.esearch(db="taxonomy", term=taxon, retmode="xml")
    id_numbers.append(Entrez.read(species_information)['IdList'][0])
    for id_n in id_numbers:
        inte = "txid" + str(id_n)
        species_information = Entrez.esearch(db="nucleotide", term=inte, retmode="xml")
        records = Entrez.read(species_information)
        for record in records['IdList']:
            genbank = Entrez.efetch(db="nucleotide", id=record, retmode="xml")
            sequences = Entrez.read(genbank)
            seqs.append(sequences)
            time.sleep(1)
    return seqs


def fasta(sequences, output):
    '''
    Writes seq records to a fasta file
    '''
    with open(output, 'w') as out:
        for record in sequences:
            try:
                ide = "".join(record[0]['GBSeq_other-seqids'])
                description = (record[0]['GBSeq_source'] + '. Length:' +
                         record[0]['GBSeq_length'])
                sequence = (record[0]['GBSeq_sequence']).upper()
                seq = Seq.Seq(sequence, IUPAC.ambiguous_dna)
                seq_complete = SeqRecord(seq, id=ide, name=ide,
                                description=description)
                SeqIO.write([seq_complete], out, 'fasta')
            except:
                pass


def main():
    inp, outp = arguments()
    seqs = getting_sequences('glaucia.ornito@gmail.com', inp)
    fasta(seqs, outp)

if __name__ == '__main__':
    main()
