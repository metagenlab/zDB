#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import accession2taxon_id


def protein_id2synonyms(protein_id):

    from Bio import Entrez
    Entrez.email = "trestan.pillonel@unil.ch"

    ncbi_id_list = accession2taxon_id.gi(protein_id, 'protein')

    handle = Entrez.elink(dbfrom='protein', db='protein', id=','.join(
        ncbi_id_list), linkname="proprotein_cdd_seqr")
    record = [i for i in Entrez.read(handle)]

    # handle = Entrez.efetch(db='protein', id=','.join(ncbi_id_list), retmax=len(ncbi_id_list))

    # record = Entrez.parse(handle)

    print record
    all_links = [i['Id'] for i in record[0]['LinkSetDb'][1]['Link']]
    # for link in all_links:
    handle = Entrez.esummary(db='protein', id=','.join(
        all_links), retmax=len(all_links))

    record = Entrez.parse(handle)
    for one_rec in record:
        print one_rec


protein_id2synonyms('ADI37996')
