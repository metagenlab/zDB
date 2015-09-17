#!/usr/bin/env python

from Bio import Entrez, SeqIO

Entrez.email = "trestan.pillonel@unil.ch"
import gbk2ffn
import gbk2faa
import urllib2
import genbank2refseq
from Bio.Alphabet import generic_nucleotide
from Bio.Alphabet import generic_dna

def get_genomic_data(ncbi_accession, genbank2refseq_id=False):
    import time

    if genbank2refseq_id:
        genbank2refseq_id = genbank2refseq.genbank2refseq(ncbi_accession)


    i = 0
    handle = None
    while not handle:
        print i
        if i == 10:
            print 'reached max iteration number, %s could not be downloaded' % ncbi_accession
            return
        try:
            handle = Entrez.efetch(db="nucleotide", id=ncbi_accession, rettype="gb", retmode="text")
        except (urllib2.URLError, urllib2.HTTPError) as e:
            print 'url error, trying again...'
            time.sleep(1)
            i+=1

    seq_records = list(SeqIO.parse(handle, "genbank"))
    for record in seq_records:
        print "writing record %s..." % record.name

        if record.seq.count("N") == len(record.seq):
            print 'no_sequences in gbk, getting fasta record %s' % record.name
            one_handle = Entrez.efetch(db="nucleotide", id=record.name, rettype="fasta", retmode="text")
            fasta_record = list(SeqIO.parse(one_handle, "fasta"))[0]
            #fasta_record.seq.alphabet = generic_nucleotide
            record.seq = fasta_record.seq
            record.seq.alphabet = generic_dna
        print record

        #if "wgs" in record.annotations:
        #    print "WGS, not writing record %s" % record.name, record.description
        #    handle.close()
        #    break
        #try:
        SeqIO.write(record, "%s.gbk" % record.name, "genbank")
        #except:
        #    print "problem writing genbank"
        try:
            SeqIO.write(record, "%s.fna" % record.name, "fasta")
        except:
            print "problem writing fasta"
        #try:
        gbk2faa.gbk2faa([record], "%s.faa" % record.name)
        #except:
        #    print "problem writing faa"
        try:
            gbk2ffn.gbk2ffn([record], "%s.ffn" % record.name)
        except:
            print "problem writing ffn"
    handle.close()




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-a",'--accession',type=str,help="NCBI accession")

    args = parser.parse_args()
    get_genomic_data(args.accession)



