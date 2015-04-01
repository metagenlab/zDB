#!/usr/bin/env python

from Bio import Entrez


Entrez.email = "trestan.pillonel@unil.ch"

def gi(ncbi_term, database="nuccore"):

    handle = Entrez.esearch(db=database, term=ncbi_term)
    record= Entrez.read(handle)

    return record["IdList"]


def genbank2refseq(ncbi_id):
    '''
    :param genbank id
    :return: refseq_id
    '''

    # get link from genbank 2 refseq
    handle = Entrez.elink(dbfrom="nuccore", db="nuccore", id=gi(ncbi_id), term="srcdb+refseq[prop]")
    record = Entrez.read(handle)

    try:
        # if no result, return None
        refseq_id = record[0]['LinkSetDb'][0]['Link'][0]['Id'] #['IdList']
    except IndexError:
        return None
    ref_seq_accession = Entrez.efetch(db="nucleotide", id=refseq_id, rettype="acc")

    return ref_seq_accession.read()


def refseq2genbank(ncbi_id):
    '''
    :param refseq_id
    :return: genbank id
    '''

    handle = Entrez.elink(dbfrom="nuccore", db="nuccore", id=gi(ncbi_id), term="srcdb+ddbj/embl/genbank[prop]")
    record = Entrez.read(handle)
    try:
        genbank_id = record[0]['LinkSetDb'][0]['Link'][0]['Id'] #['IdList']
    except IndexError:
        return None
    genbank_accession = Entrez.efetch(db="nucleotide", id=genbank_id, rettype="acc")

    return genbank_accession.read()


def identify_id(ncbi_id, database = "nuccore"):

    '''
    :param refseq_id
    :return:
    '''
    import SeqIO

    my_gi = gi(ncbi_id, database)

    handle = Entrez.efetch(db=database, id=my_gi, rettype="gb")

    list(SeqIO.parse(handle, "genbank"))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",'--seq_id_genbank', default=False, type=str, help="genbank2refseq")
    parser.add_argument("-y",'--seq_id_refseq', default=False, type=str, help="refseq2genbank")

    args = parser.parse_args()

    if args.seq_id_genbank:
        print genbank2refseq(args.seq_id_genbank)
    if args.seq_id_refseq:
        print refseq2genbank(args.seq_id_refseq)
