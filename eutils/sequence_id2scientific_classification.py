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
    handle = Entrez.elink(dbfrom="nuccore", db="taxonomy", id=uid)
    record = Entrez.read(handle)

    link = record[0]["LinkSetDb"][0]["Link"][0]["Id"]

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

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",'--seq_id', type=str, help="sequence ncbi id")
    parser.add_argument("-p", '--protein_seq', action="store_true", help="Protein sequence (default=False)")

    args = parser.parse_args()


    print sequence_id2scientific_classification(args.seq_id, protein=args.protein_seq)


