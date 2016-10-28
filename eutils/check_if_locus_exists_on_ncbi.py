#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

def locus_on_ncbi(locus_tag):

    from Bio import Entrez

    Entrez.email = "trestan.pillonel@unil.ch"


    handle = Entrez.esearch(db='nucleotide', term=locus_tag)

    record1 = Entrez.read(handle)
    print record1['Count']
    if int(record1['Count']) == 0:
        return False
    else:
        return True

