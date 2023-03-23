#!/usr/bin/env python

from Bio import Entrez


Entrez.email = "trestan.pillonel@unil.ch"

def gi(ncbi_term, database="nuccore"):

    handle = Entrez.esearch(db=database, term=ncbi_term)
    record= Entrez.read(handle)

    return record["IdList"]


def genbank2refseq(ncbi_id, database="nuccore"):
    '''
    :param genbank id
    :return: refseq_id
    '''

    # get link from genbank 2 refseq
    handle = Entrez.elink(dbfrom="nuccore", db=database, id=gi(ncbi_id), term="srcdb+refseq[prop]")
    record = Entrez.read(handle)

    try:
        # if no result, return None
        refseq_id = record[0]['LinkSetDb'][0]['Link'][0]['Id'] #['IdList']
    except IndexError:
        return None
    ref_seq_accession = Entrez.efetch(db=database, id=refseq_id, rettype="acc")

    return ref_seq_accession.read()


def refseq2genbank(ncbi_id, database="nuccore"):
    '''
    :param refseq_id
    :return: genbank id
    '''

    try:
        int(ncbi_id)
        print 'ncbi gi!'
        handle = Entrez.elink(dbfrom=database, db=database, id=ncbi_id, term="srcdb+ddbj/embl/genbank[prop]")
    except:
        handle = Entrez.elink(dbfrom=database, db=database, id=gi(ncbi_id), term="srcdb+ddbj/embl/genbank[prop]")
    record = Entrez.read(handle)
    print record
    try:
        genbank_id = record[0]['LinkSetDb'][0]['Link'][0]['Id'] #['IdList']
    except IndexError:
        return None
    genbank_accession = Entrez.efetch(db=database, id=genbank_id, rettype="acc")

    return genbank_accession.read()


def identify_id(ncbi_id, database = "nuccore"):

    '''
    :param refseq_id
    :return:
    '''
    from Bio import SeqIO

    my_gi = gi(ncbi_id, database)

    handle = Entrez.efetch(db=database, id=my_gi, rettype="gb")

    return list(SeqIO.parse(handle, "genbank"))


if __name__ == '__main__':
    import argparse
    import sys
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-g",'--seq_id_genbank', default=False, type=str, help="genbank2refseq")
    parser.add_argument("-r",'--seq_id_refseq', default=False, type=str, help="refseq2genbank")
    parser.add_argument("-i",'--info', default=False, type=str, help="Get search result from input identifier")
    parser.add_argument("-d",'--database', default="nuccore", type=str, help="Target database (genome, protein, nucccore, bioproject,...")

    args = parser.parse_args()

    if args.seq_id_genbank:
        sys.stdout.write(genbank2refseq(args.seq_id_genbank, args.database))
    if args.seq_id_refseq:
        sys.stdout.write(refseq2genbank(args.seq_id_refseq, args.database))
    if args.info:
        SeqIO.write(identify_id(args.seq_id_refseq, args.database), sys.stdout, "genbank")
