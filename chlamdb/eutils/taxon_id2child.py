#!/usr/bin/env python

from Bio import Entrez
Entrez.email = "trestan.pillonel@unil.ch"

def taxon_id2child(taxon_id):
    '''
    get all child taxon_id from one taxon_id
    :return:
    '''

    s_term = 'txid%s[Organism:exp]'
    handle = Entrez.esearch(db="taxonomy", term=s_term % taxon_id, retmax=100000)
    record = Entrez.read(handle)

    return record['IdList']

#taxid_list =  taxon_id2child(813)
#print '813' in taxid_list