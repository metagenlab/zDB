#!/usr/bin/env python

from Bio import Entrez, SeqIO
import eutils
Entrez.email = "trestan.pillonel@unil.ch"

#handle_sequences = Entrez.elink(dbfrom="nuccore", db="bioproject", id=470217946, term="srcdb+ddbj/embl/genbank[prop]")
#record_sequences =  Entrez.read(handle_sequences)
handle_bioproject = Entrez.elink(dbfrom="nuccore", db="bioproject", id=470217946)
record_bioproject = Entrez.read(handle_bioproject)

print record_bioproject
bioproject_link = record_bioproject[0]["LinkSetDb"][0]["Link"][0]["Id"]
print "buiolink", bioproject_link
handle_sequences = Entrez.elink(dbfrom="bioproject", db="nuccore", id=bioproject_link)
record_sequences =  Entrez.read(handle_sequences)

refseq_links = [link["Id"] for link in record_sequences[0]["LinkSetDb"][0]["Link"]]

print "all", refseq_links

handle_refseq = Entrez.elink(dbfrom="nuccore", db="nuccore", id=refseq_links[0], term="srcdb+ddbj/embl/genbank[prop]")
record_refseq = Entrez.read(handle_refseq)
print "refseq2genbank", record_refseq

