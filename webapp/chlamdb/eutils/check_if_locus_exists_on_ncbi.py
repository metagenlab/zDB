#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

def locus_on_ncbi(locus_tag):

    from Bio import Entrez
    from Bio import SeqIO
    from dateutil.parser import parse

    Entrez.email = "trestan.pillonel@unil.ch"

    handle = Entrez.esearch(db='nucleotide', term=locus_tag)

    record1 = Entrez.read(handle)
    print record1['IdList']

    if int(record1['Count']) == 0:
        return False
    else:
        handle2 = Entrez.efetch(
            db='nucleotide', id=record1['IdList'][0], rettype="gb", retmode="text")
        record2 = SeqIO.read(handle2, "genbank")
        dt = parse(record2.annotations['date'])

        return dt
