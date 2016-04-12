#!/usr/bin/env python
# utf-8


"""
task 1 of assignment 18
Created by Pramod Pantha on April 10, 2016.
Copyright 2016 Pramod Pantha. All right reserved.
reference:
http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc65
class slide
"""


import time
import csv
from Bio import Entrez
import argparse
Entrez.email = "ppanth1@lsu.edu"
species_name = 'species-list.txt'
outfilename = 'outfile.csv'

csv.register_dialect(
    'mydialect',
    delimiter=',',
    quotechar='"',
    doublequote=True,
    skipinitialspace=True,
    lineterminator='\n',
    quoting=csv.QUOTE_MINIMAL)
# reference: https://www.getdatajoy.com/examples/python-data-analysis/read-and-write-a-csv-file-with-the-csv-module
f1 = []
infilename = 'species-list.txt'
# species_list.txt is modified as per discussionin the slack channel(Genus of
# two species in line 23 and 27 has been changed). use file I have submitted

def get_data(species_name):
    esearch_query = Entrez.esearch(db="taxonomy", term=species_name, retmode="xml")
    esearch_result = Entrez.read(esearch_query)
    bb_id = esearch_result['IdList']
    data = []
    names = ['superclass','class','subclass','infraclass','superorder','order','superfamily','family','genus']

    if bb_id:
        genbank_entries = Entrez.efetch(db="taxonomy", id=bb_id, retmode="xml")
        esearch_result = Entrez.read(genbank_entries)
        for i in range(len(names)):
            found = 0
            for d in esearch_result[0]['LineageEx']:
                if d['Rank'] == names[i]:
                    data.append(d['ScientificName'])
                    found = 1
            if found == 0:
                data.append(' ')
    return data


def file_writer(outfilename):
    writefile = open(outfilename, 'w')
    writer = csv.DictWriter(writefile, fieldnames=['superclass','class','subclass','infraclass','superorder','order','superfamily',
                                                     'family','genus'])
    writer.writeheader()
    datawriter = csv.writer(writefile, dialect='mydialect')
    with open(infilename, 'r') as readfile:
        for line in readfile:
            f1.append(line[:-1])
            taxo = (get_data(line[:-1]))
            print(taxo)
            writer.writerow({'superclass':taxo[0],'class':taxo[1],'subclass':taxo[2],'infraclass':taxo[3],'superorder':taxo[4],'order':taxo[5]
                             ,'superfamily':taxo[6],'family':taxo[7],'genus':taxo[8]})
            time.sleep(1)
            # limits the speed of request
    writefile.close()
    readfile.close()


def main():
    parse = argparse.ArgumentParser()
    parse.add_argument("inputfile", help="give name of your txt filename")
    parse.add_argument("outputfile", help="give your output file name")
    args = parse.parse_args()
    specieslist = args.inputfile
    csvfile = args.outputfile + '.csv'
    get_data(specieslist, csvfile)


if __name__ == '__main__':
    main()
