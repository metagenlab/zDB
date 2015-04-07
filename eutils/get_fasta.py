#!/usr/bin/env python

# get a set of nucleotide or protein sequences from ncbi using eutils
# Author: Claire Bertelli (claire.bertelli[@]gmail.com)
# Date: 04.2015
# ---------------------------------------------------------------------------

import argparse
from Bio import Entrez, SeqIO
import sys
from pandas import read_table

Entrez.email = "claire.bertelli@chuv.ch"

def download_seq_entrez(accession, multifasta, seq_type):
    # retrieves information about entry in ncbi using Entrez


    # if the user requests individual fasta file per accession
    if multifasta=="":
        for one_accession in accession:
            output_handle = open("%s.fasta" %one_accession, "w")
            handle = Entrez.efetch(db=seq_type, id=one_accession, rettype="fasta")
            output_handle.write(handle.read())
            output_handle.close()

     # if the user requests a single multifasta file
    else:
        output_handle = open("%s.fasta" %multifasta, "w")
        handle = Entrez.efetch(db=seq_type, id=accession, rettype="fasta")
        output_handle.write(handle.read())

        output_handle.close()


# reads table and retrieve GI in the column specified by the user
def table_to_gi(file, column):
    df = read_table(file, header=None, sep='\t')
    gi_list=list(df[column-1])
    gi_list=[str(i)for i in gi_list]
    return gi_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-a",'--accession',type=str,  help="file with accession numbers")
    parser.add_argument("-t",'--seq_type',type=str, default="protein", help="Type of sequence to be downloaded: protein, gene, nuccore. Default: protein")
    parser.add_argument("-m",'--multifasta',type=str, default="", help="For multiple accession numbers, enable to download all sequences in a multifasta file. Please provide a prefix for the output.")
    parser.add_argument("-c",'--column',type=int, default="", help="Input is a file containing a tabulated table containing the GI of sequences to retrieve in column c. Please specify column number.")

    args = parser.parse_args()

    if args.column!="":
        accession=table_to_gi(args.accession, args.column)
    else:
        accession=list(args.accession)

    sys.stdout.write("getting protein sequences %s ... \n" % accession)
    download_seq_entrez(accession, args.multifasta, args.seq_type)