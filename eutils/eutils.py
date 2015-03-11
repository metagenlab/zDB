#!/usr/bin/env python

from Bio import Entrez, SeqIO
Entrez.email = "trestan.pillonel@unil.ch"
import gbk2faa



def gbk2ffn(seq_record, outname):
    output_handle = open(outname, "w")
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            #print seq_feature
            assert len(seq_feature.qualifiers['translation'])==1
            # gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
            try:
                output_handle.write(">gi|%s|ref|%s| %s [%s]\n%s\n" % (
                        seq_feature.qualifiers["db_xref"][0].split(":")[1],
                        seq_record.id,
                        seq_feature.qualifiers["note"][0],
                        seq_record.description,
                        seq_feature.extract(seq_record.seq)))
            except:
                output_handle.write(">gi|%s|ref|%s| [%s]\n%s\n" % (
                        seq_feature.qualifiers["db_xref"][0].split(":")[1],
                        seq_record.id,
                        seq_record.description,
                        seq_feature.extract(seq_record.seq)))



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
        gbk2faa.gbk2ffn([record], "%s.ffn" % record.name)
        #except:
        #    print "problem writing ffn"
    handle.close()




if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-a",'--accession',type=str,help="NCBI accession")

    args = parser.parse_args()
    get_genomic_data(args.accession)



