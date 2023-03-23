#!/usr/bin/env python

from Bio import Entrez


Entrez.email = "trestan.pillonel@unil.ch"

def pmid2abstract_info(pmid):
    from Bio import Medline
    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline",
                           retmode="text")
    record = Medline.read(handle)

    print record
    pmid_data = {}
    pmid_data["title"] = record.get("TI", "?")
    pmid_data["authors"] = record.get("AU", "?")
    pmid_data["source"] = record.get("SO", "?")
    pmid_data["abstract"] = record.get("AB", "?")
    pmid_data["pmid"] = pmid

    return pmid_data

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