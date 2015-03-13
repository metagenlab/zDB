#!/usr/bin/env python

from Bio import Entrez, SeqIO
Entrez.email = "trestan.pillonel@unil.ch"
import gbk2ffn



def get_genomic_data(ncbi_accession):

    handle = Entrez.efetch(db="nucleotide", id=ncbi_accession, rettype="gb", retmode="text")
    seq_records = list(SeqIO.parse(handle, "genbank"))
    for record in seq_records:
        print "writing record %s..." % record.name

        if "wgs" in record.annotations:
            print "WGS, not writing record %s" % record.name, record.description
            handle.close()
            break
        try:
            SeqIO.write(record, "%s.gbk" % record.name, "genbank")
        except:
            print "problem writing genbank"
        try:
            SeqIO.write(record, "%s.fna" % record.name, "fasta")
        except:
            print "problem writing fasta"
        try:
            gbk2faa.gbk2faa([record], "%s.faa" % record.name)
        except:
            print "problem writing faa"
        #try:
        gbk2ffn.gbk2ffn([record], "%s.ffn" % record.name)
        #except:
        #    print "problem writing ffn"
    handle.close()




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-a",'--accession',type=str,help="NCBI accession")

    args = parser.parse_args()
    get_genomic_data(args.accession)



