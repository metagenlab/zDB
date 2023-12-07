#!/usr/bin/env python

from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_nucleotide
import genbank2refseq
import urllib2
import gbk2faa
import gbk2ffn
from Bio import Entrez, SeqIO

Entrez.email = "trestan.pillonel@unil.ch"


def chunks(l, n):
    "return sublists of l of minimum length n (work subdivision for the subprocesing module"
    return [l[i:i+n] for i in range(0, len(l), n)]


def get_protein_sequences(crispdata):

    import time
    from Bio.SeqRecord import SeqRecord

    mychunks = chunks(crispdata.keys(), 10)

    for i, chunk in enumerate(mychunks):
        print "%s/%s" % (i, len(mychunks))
        table = 1

        print chunk
        i = 0
        handle = None
        while not handle:
            print i
            if i == 10:
                print 'reached max iteration number, %s could not be downloaded' % chunk
                break
            try:
                handle = Entrez.efetch(db="nucleotide", id=','.join(
                    chunk), rettype="fasta", retmode="text")
            except (urllib2.URLError, urllib2.HTTPError) as e:
                print e
                print 'url error, trying again...'
                time.sleep(1)
                i += 1
        if not handle:
            continue

        seq_records = list(SeqIO.parse(handle, "fasta"))

        for seq_record in seq_records:
            # print seq_record
            ncbi_accession = seq_record.id.split('|')[-2].split('.')[0]
            print 'accession', ncbi_accession
            for gene_data in crispdata[ncbi_accession]:

                localisation = gene_data[0].split('..')
                strand = gene_data[1]
                start = int(localisation[0])-1
                end = int(localisation[1])
                gene_name = ncbi_accession + '_' + gene_data[2]

                # print seq_record[start:end].seq

                if strand == '-':
                    print ncbi_accession, 'strand -', gene_name, start, end, table, seq_record.description
                    if "Mycoplasma" in seq_record.description:
                        print 'Myco!'
                        table = 4
                    seq = seq_record[start:end].seq.reverse_complement(
                    ).translate(table=table)
                    if '*' in seq[0:-1] and table == 1:
                        print 'problem with translation, trying to identitfy translation table'
                        handle = Entrez.efetch(
                            db="nucleotide", id=ncbi_accession, rettype="gbwithparts", retmode="text")
                        seq_record = SeqIO.read(handle, "genbank")
                        for feature in seq_record.features:
                            if feature.type == 'CDS':
                                table = int(
                                    feature.qualifiers["transl_table"][0])
                        seq = seq_record[start:end].seq.reverse_complement(
                        ).translate(table=table)
                        print ncbi_accession, 'strand -', gene_name, start, end, table, seq_record.description
                        # print seq

                    if '*' in seq[0:-1]:
                        print 'changing start and stop: %s %s' % (start+1, end+1)
                        seq = seq_record[(
                            start+1):(end+1)].seq.reverse_complement().translate(table=table)
                        # print seq
                    if '*' in seq[0:-1]:
                        print 'changing again start and stop: %s %s' % (start, end+1)
                        seq = seq_record[start:(
                            end+1)].seq.reverse_complement().translate(table=table)
                        # print seq
                    if '*' in seq[0:-1]:
                        print 'changing again start and stop: %s %s' % (start+2, end+1)
                        seq = seq_record[(
                            start+2):(end+1)].seq.reverse_complement().translate(table=table)
                        # print seq

                else:
                    if "Mycoplasma" in seq_record.description:
                        print 'Myco!'
                        table = 4
                    print ncbi_accession, 'strand +', gene_name, start, end, table, seq_record.description
                    seq = seq_record[start:end].seq.translate(table)

                    if '*' in seq[0:-1] and table == 1:
                        print 'problem with translation, trying to identitfy translation table'
                        handle = Entrez.efetch(
                            db="nucleotide", id=ncbi_accession, rettype="gbwithparts", retmode="text")
                        seq_record = SeqIO.read(handle, "genbank")
                        for feature in seq_record.features:
                            if feature.type == 'CDS':
                                table = int(
                                    feature.qualifiers["transl_table"][0])
                        seq = seq_record[start:end].seq.reverse_complement(
                        ).translate(table=table)
                        print ncbi_accession, 'strand -', gene_name, start, end, table, seq_record.description
                        # print seq

                    if '*' in seq[0:-1]:
                        print 'changing start and stop: %s %s' % (start+1, end+1)
                        seq = seq_record[(start+1):(end+1)
                                         ].seq.translate(table=table)
                        # print seq
                    if '*' in seq[0:-1]:
                        print 'changing again start and stop: %s %s' % (start, end+1)
                        # print seq
                        seq = seq_record[start:(
                            end+1)].seq.translate(table=table)
                    if '*' in seq[0:-1]:
                        print 'changing again start and stop: %s %s' % (start+2, end+1)
                        seq = seq_record[(start+2):(end+1)
                                         ].seq.translate(table=table)

                    if '*' in seq[0:-1]:
                        print 'changing again start and stop: %s %s' % (start+3, end+2)
                        seq = seq_record[(start+3):(end+2)
                                         ].seq.translate(table=table)

                if '*' in seq:
                    seq = seq[0:-1]
                if '*' in seq:
                    print '############  could not translate sequence ##########', seq_record.description, start, end
                    break

                newrecord = SeqRecord(
                    seq, id=gene_name, description=seq_record.description)

                with open("all_seqs.fa", 'a') as f:
                    SeqIO.write(newrecord, f, 'fasta')

                handle.close()


def parse_crisprfile(crisp_table):
    accession2localisation = {}
    with open(crisp_table, 'r') as f:
        for line in f:
            data = line.rstrip().split('\t')
            if not '..' in data[1]:
                continue
            else:
                if data[0] not in accession2localisation:
                    accession2localisation[data[0]] = [data[1:4]]
                else:
                    accession2localisation[data[0]].append(data[1:4])
    return accession2localisation


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", '--accession', type=str, help="NCBI accession")
    parser.add_argument("-t", '--crisp_table', type=str, help="crisp_table")

    args = parser.parse_args()
    get_protein_sequences(parse_crisprfile(args.crisp_table))
    # get_genomic_data(args.accession)
