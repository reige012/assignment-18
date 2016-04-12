#!/usr/bin/env python
# encoding: utf-8

"""
Script prints taxonomic information for each species in a given list. The
following information will be obtained from the NCBI Taxonomy database:
    Superclass, class, subclass, infraclass, superorder, order, superfamily,
    family, and genus
If the information is not available with certainty then an empty bracket will
appear where that information should be located in the hierarchy listed above.
In general there is a lot of missing information due to the lack of database
info being ranked appropriately.


Notes: Input file must list each species name on a separate line for this
script to work appropriately.Output file must end in .csv (Example: output.csv)

Edited by Alicia Reigel. 9 April 2016.
Copyright Alicia Reigel. Louisiana State University. 9 April 2016. All
rights reserved.

"""

import time
import argparse
from Bio import Entrez


def parser_get_args():
    """Collect the path to the file of species names and name of output file"""
    parser = argparse.ArgumentParser(
        description="""Input the full path to species name file and desired
            output file name"""
        )
    parser.add_argument(
            '--filepath',
            required=True,
            type=str,
            help='Enter the path to the species name file.'
        )
    parser.add_argument(
            '--outputfile',
            required=True,
            type=str,
            help='Enter the desired name for the output file. Must end in .csv'
        )
    return parser.parse_args()


def get_superclass(lineageex):
    superclass_dict = {}
    for entry in lineageex:
        for key, value in entry.items():
            # finds superclass and adds to dictionary
            # adds nothing if data not available
            if key == ("Rank") and value == ("superclass"):
                superclass_dict["Superclass"] = entry["ScientificName"]
            else:
                pass
    return superclass_dict


def get_class(lineageex):
    class_dict = {}
    for entry in lineageex:
            for key, value in entry.items():
                # finds class and adds to dictionary
                # adds nothing if data not available
                if key == ('Rank') and value == ('class'):
                    class_dict['Class'] = entry['ScientificName']
                else:
                    pass
    return class_dict


def get_subclass(lineageex):
    subclass_dict = {}
    for entry in lineageex:
        for key, value in entry.items():
            # finds subclass and adds to dictionary
            # adds nothing if not available
            if key == ("Rank") and value == ("subclass"):
                subclass_dict['Subclass'] = entry["ScientificName"]
            else:
                pass
    return subclass_dict


def get_infraclass(lineageex):
    infraclass_dict = {}
    for entry in lineageex:
        for key, value in entry.items():
            # finds infraclass and adds to dictionary
            # adds nothing if not available
            if key == ("Rank") and value == ("infraclass"):
                infraclass_dict['Infraclass'] = entry["ScientificName"]
            else:
                pass
    return infraclass_dict


def get_superorder(lineageex):
    superorder_dict = {}
    for entry in lineageex:
        for key, value in entry.items():
            # finds superorder and adds to dictionary
            # adds nothing if not available
            if key == ("Rank") and value == ("superorder"):
                superorder_dict['Superorder'] = entry["ScientificName"]
            else:
                pass
    return superorder_dict


def get_order(lineageex):
    order_dict = {}
    for entry in lineageex:
        for key, value in entry.items():
            # finds order and adds to dictionary
            # adds nothing if not available
            if key == ("Rank") and value == ("order"):
                order_dict['Order'] = entry["ScientificName"]
            else:
                pass
    return order_dict


def get_superfamily(lineageex):
    superfamily_dict = {}
    for entry in lineageex:
        for key, value in entry.items():
            # finds superfamily and adds to dictionary
            # adds nothing if not available
            if key == ("Rank") and value == ("superfamily"):
                superfamily_dict['Superfamily'] = entry["ScientificName"]
            else:
                pass
    return superfamily_dict


def get_family(lineageex):
    family_dict = {}
    for entry in lineageex:
        for key, value in entry.items():
            # finds family and adds to dictionary
            # adds nothing if not available
            if key == ("Rank") and value == ("family"):
                family_dict['Family'] = entry["ScientificName"]
            else:
                pass
    return family_dict


def get_genus(lineageex):
    genus_dict = {}
    for entry in lineageex:
        for key, value in entry.items():
            # finds genus and adds to dictionary
            # adds nothing if not available
            if key == ("Rank") and value == ("genus"):
                genus_dict['Genus'] = entry["ScientificName"]
            else:
                pass
    return genus_dict


def main():
    args = parser_get_args()
    with open(args.filepath, 'r') as species_names:
        species_list = species_names.read().splitlines()
        for species in species_list:
            time.sleep(1)
            Entrez.email = "areige1@lsu.edu"
            search_query = Entrez.esearch(db="taxonomy", term=species, retmode="xml")
            result = Entrez.read(search_query)
            species_id = result['IdList']
            species_id = "txid{}".format(species_id[0])
            taxonomy_data = Entrez.efetch(db="taxonomy", id=species_id, retmode="xml")
            data = Entrez.read(taxonomy_data)
            lineageex = data[0]["LineageEx"]
            species_dict = {}
            species_dict['Species'] = species
            superclass_dict = get_superclass(lineageex)
            class_dict = get_class(lineageex)
            subclass_dict = get_subclass(lineageex)
            infraclass_dict = get_infraclass(lineageex)
            superorder_dict = get_superorder(lineageex)
            order_dict = get_order(lineageex)
            superfamily_dict = get_superfamily(lineageex)
            family_dict = get_family(lineageex)
            genus_dict = get_genus(lineageex)
            with open(args.outputfile, 'a') as output_file:
                output_file.write(str((species_dict, superclass_dict, class_dict, subclass_dict, infraclass_dict, superorder_dict, order_dict, superfamily_dict, family_dict, genus_dict)))
                output_file.write('\n\n')

if __name__ == '__main__':
    main()
