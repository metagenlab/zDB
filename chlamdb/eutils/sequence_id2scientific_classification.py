#!/usr/bin/env python

from Bio import Entrez, SeqIO


Entrez.email = "trestan.pillonel@unil.ch"


def sequence_id2scientific_classification(ncbi_id, protein=False):
    '''
    :param ncbi_id
    :param key
    :return: complete taxonomic affiliation
    '''
    if protein:
        handle = Entrez.esearch(db="protein", term=ncbi_id)
    else:
        handle = Entrez.esearch(db="nucleotide", term=ncbi_id)
    record = Entrez.read(handle)
    try:
        uid = record["IdList"][0]
    except IndexError:
        return None
    if protein:
        handle = Entrez.elink(dbfrom="protein", db="taxonomy", id=uid)

    else:
        handle = Entrez.elink(dbfrom="nuccore", db="taxonomy", id=uid)
    record = Entrez.read(handle)
    try:
        link = record[0]["LinkSetDb"][0]["Link"][0]["Id"]
    except IndexError:
        print (record)

    handle = Entrez.efetch(db="taxonomy", id=link, rettype="xml")
    records = Entrez.parse(handle)
    for record in records:
        # get scientific names of all classification levels
        classification = {}

        for level in record['LineageEx']:

            rank = level["Rank"]
            scientific_name = level["ScientificName"]
            classification[rank] = scientific_name
        return classification

def taxon_id2scientific_classification(taxon_id_list, taxon2rank=False):
    '''
    :param ncbi taxon id
    :param key
    :return: dictionary with taxon id as primary key and a nested dictionnary as value. The nested dictionary contain
    the classification in the form:

    {'superkingdom': 'bacteria'; 'phylum': chlamydiae;...}

    NOTE: if taxon is not completely classified (i.e. unclassified bacteria), it will return the available data only

    {'superkingdom': 'bacteria'; 'no rank': 'unclassified bacteria'}

    '''

    if type(taxon_id_list) is not list:
        raise TypeError('Expect a list of taxon id(s)')

    # check if N/A are not presents in the taxon list (i.e can happen with blast results)
    new_taxon = [value for value in taxon_id_list if value != 'N/A']

    merged_taxons = ','.join(new_taxon)
    #print 'merged taxonss', merged_taxons
    handle = Entrez.efetch(db="taxonomy", id=merged_taxons, rettype="xml")
    records = Entrez.parse(handle)
    try:
        records = [i for i in records]
    except RuntimeError:
        return None
    all_classifications = {}

    if taxon2rank:

        taxon2rank = {}

        for taxon, record in zip(new_taxon, records):
            taxon2rank[taxon] = [record['Rank'], record['ScientificName']]
        return taxon2rank


    else:

        for taxon, record in zip(new_taxon, records):
            # get scientific names of all classification levels
            classification = {}
            #print record
            #print record['LineageEx']
            #print 'rank', record['Rank'], record['ScientificName']
            for level in record['LineageEx']:
                print('##', level)
                rank = level["Rank"]
                scientific_name = level["ScientificName"]
                classification[rank] = scientific_name
            classification[record['Rank']] = record['ScientificName']
            all_classifications[taxon] = classification
        return all_classifications


def taxon_id2scientific_classification_and_taxids(taxon_id_list, taxon2rank=False):
    '''
    :param ncbi taxon id
    :param key
    :return: dictionary with taxon id as primary key and a nested dictionnary as value. The nested dictionary contain
    the classification in the form:

    {'superkingdom': 'bacteria'; 'phylum': chlamydiae;...}

    NOTE: if taxon is not completely classified (i.e. unclassified bacteria), it will return the available data only

    {'superkingdom': 'bacteria'; 'no rank': 'unclassified bacteria'}

    '''

    if type(taxon_id_list) is not list:
        raise TypeError('Expect a list of taxon id(s)')

    # check if N/A are not presents in the taxon list (i.e can happen with blast results)
    new_taxon = [value for value in taxon_id_list if value != 'N/A']

    merged_taxons = ','.join(new_taxon)
    #print 'merged taxonss', merged_taxons
    handle = Entrez.efetch(db="taxonomy", id=merged_taxons, rettype="xml")
    records = Entrez.parse(handle)
    try:
        records = [i for i in records]
    except RuntimeError:
        return None
    all_classifications = {}

    if taxon2rank:

        taxon2rank = {}

        for taxon, record in zip(new_taxon, records):
            taxon2rank[taxon] = [record['Rank'], record['ScientificName']]
        return taxon2rank


    else:

        for taxon, record in zip(new_taxon, records):
            # get scientific names of all classification levels
            classification = {}
            #print record
            #print record['LineageEx']
            #print 'rank', record['Rank'], record['ScientificName']
            for level in record['LineageEx']:
                rank = level["Rank"]
                scientific_name = level["ScientificName"]
                level_taxid = level['TaxId']
                classification[rank] = [scientific_name, level_taxid]
            classification[record['Rank']] = [record['ScientificName'], record['TaxId']]
            all_classifications[taxon] = classification
        return all_classifications

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",'--seq_id', type=str, help="sequence ncbi id")

    parser.add_argument("-p", '--protein_seq', action="store_true", help="Protein sequence (default=False, search in nucleotide databases)")

    args = parser.parse_args()


    print (sequence_id2scientific_classification(args.seq_id, protein=args.protein_seq))


