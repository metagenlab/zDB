#!/usr/bin/env python

from Bio import Entrez


Entrez.email = "trestan.pillonel@unil.ch"

def pmid2abstract_info(pmid):
    handle = Entrez.esummary(db="pubmed", id="%s" % pmid, retmode="xml")
    records = Entrez.parse(handle)
    for record in records:
        print record['Title']
        print record['FullJournalName']
        print record['SO']

def doi2abstract_info():
    pass

def sequence_accession2abstract_info():
    pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i",'--pmid_id', default=False, type=str, help="pubmed id (mpid)")

    args = parser.parse_args()
    pmid2abstract_info(args.pmid_id)