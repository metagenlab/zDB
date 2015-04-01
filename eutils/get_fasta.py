#!/usr/bin/env python

# get a set of nucleotide or protein sequences from ncbi using eutils
# TODO replace prints by sys.stdout/err.write
# Author: Claire Bertelli (claire.bertelli[@]gmail.com)
# Date: 04.2015
# ---------------------------------------------------------------------------


from Bio import Entrez, SeqIO
import sys

Entrez.email = "claire.bertelli@unil.ch"

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


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-a",'--accession',type=str, nargs='+',  help="accession numbers")
    parser.add_argument("-t",'--seq_type',type=str, default="protein", help="Type of sequence to be downloaded: protein, gene, nuccore. Default: protein")
    parser.add_argument("-m",'--multifasta',type=str, default="", help="For multiple accession numbers, enable to download all sequences in a multifasta file. Please provide a prefix for the output.")


    args = parser.parse_args()

    if args.seq_type=="protein":
        sys.stdout.write("getting protein sequences %s ... \n" % args.accession)
        download_seq_entrez(args.accession, args.multifasta, args.seq_type)