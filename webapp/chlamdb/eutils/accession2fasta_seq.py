#!/usr/bin/env python

from Bio import Entrez


Entrez.email = "trestan.pillonel@unil.ch"


def gi(ncbi_term, database="nuccore", retmax=20):

    handle = Entrez.esearch(db=database, term=ncbi_term, retmax=retmax)
    record = Entrez.read(handle)

    return record["IdList"]


def accession2record(ncbi_id_list, db="protein"):
    '''
    :param genbank/refseq protein accession
    :return: dictionnary protein_id2description
    '''
    from Bio import SeqIO
    import re

    gi_list = gi(','.join(ncbi_id_list), db, retmax=len(ncbi_id_list))

    # print "gi_list", ','.join(ncbi_gi_list)
    # get link from genbank 2 refseq
    handle = Entrez.efetch(db=db, id=','.join(
        gi_list), rettype="gb", retmode="text")

    records = [i for i in SeqIO.parse(handle, "genbank")]

    return records


if __name__ == '__main__':
    import argparse
    import sys
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--seq_id_genbank', default=False,
                        type=str, help="genbank2refseq", nargs="+")
    parser.add_argument("-d", '--ncbi_database', default="nucleotide",
                        type=str, help="database to search (protein/nucleotide/...)")
    parser.add_argument("-o", '--out_name', default="out.fa",
                        type=str, help="out name")
    parser.add_argument("-a", '--header_name', default=False, type=str,
                        help="header (default False, uses from fetched data)")
    args = parser.parse_args()

    rec = accession2record(args.seq_id_genbank, args.ncbi_database)

    if args.header_name:
        rec[0].id = args.header_name
    SeqIO.write(rec, args.out_name, "fasta")
