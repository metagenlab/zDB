#!/usr/bin/env python

from Bio import Entrez


Entrez.email = "trestan.pillonel@unil.ch"


def genbank2refseq(ncbi_id):
    '''
    :param genbank id
    :return: refseq_id
    '''
    handle = Entrez.esearch(db="nuccore", term=ncbi_id)
    record = Entrez.read(handle)

    # get gi
    genome_id = record["IdList"]

    # get link from genbank 2 refseq
    handle = Entrez.elink(dbfrom="nuccore", db="nuccore", id=genome_id, term="srcdb+refseq[prop]")
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
    handle = Entrez.esearch(db="nuccore", term=ncbi_id)
    record = Entrez.read(handle)

    # get genome overview id
    genome_id = record["IdList"]

    handle = Entrez.elink(dbfrom="nuccore", db="nuccore", id=genome_id, term="srcdb+ddbj/embl/genbank[prop]")
    record = Entrez.read(handle)
    try:
        genbank_id = record[0]['LinkSetDb'][0]['Link'][0]['Id'] #['IdList']
    except IndexError:
        return None
    genbank_accession = Entrez.efetch(db="nucleotide", id=genbank_id, rettype="acc")

    return genbank_accession.read()



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