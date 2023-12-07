#!/usr/bin/env python

# get a set of nucleotide or protein sequences from ncbi using Entrez
# Author: Claire Bertelli (claire.bertelli[@]gmail.com)
# Date: 04.2015
# ---------------------------------------------------------------------------

import argparse
from Bio import Entrez
import sys
from pandas import read_table
import re
from shell_command import shell_command

Entrez.email = "claire.bertelli@chuv.ch"


def download_seq_entrez(accession, multifasta, seq_type, append):
    # retrieves information about entry in ncbi using Entrez

    # if the user requests individual fasta file per accession
    if multifasta == "":
        for one_accession in accession:
            output_name = re.sub(
                "([a-zA-Z_0-9]+)\.([a-zA-Z]+)", "\\1.fasta", one_accession)
            output_handle = open(output_name, "w")
            handle = Entrez.efetch(
                db=seq_type, id=one_accession, rettype="fasta")
            output_handle.write(handle.read())
            output_handle.close()

     # if the user requests a single multifasta file
    else:
        # if the user required to add the blast file
        if append != "":
            output_handle = open("%s_only.fasta" % multifasta, "w")
            handle = Entrez.efetch(db=seq_type, id=accession, rettype="fasta")
            output_handle.write(handle.read())
            output_handle.close()

            cmd = "cat %s %s_only.fasta > %s.fasta" % (
                append, multifasta, multifasta)
            print cmd
            shell_command(cmd)
            cmd = "rm %s_only.fasta" % multifasta
            print cmd
            shell_command(cmd)

        else:
            output_handle = open("%s.fasta" % multifasta, "w")
            handle = Entrez.efetch(db=seq_type, id=accession, rettype="fasta")
            output_handle.write(handle.read())
            output_handle.close()


# reads table and retrieve GI in the column specified by the user
def table_to_gi(file, column):
    df = read_table(file, header=None, sep='\t')
    gi_list = list(df[column-1])
    gi_list = [str(i)for i in gi_list]
    return gi_list


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-a", '--accession', type=str,
                        help="File with accession numbers")
    parser.add_argument("-t", '--seq_type', type=str, default="protein",
                        help="Type of sequence to be downloaded: protein, gene, nuccore. Default: protein")
    parser.add_argument("-m", '--multifasta', type=str, default="",
                        help="For multiple accession numbers, enable to download all sequences in a multifasta file. Please provide a prefix for the output.")
    parser.add_argument("-c", '--column', type=int, default="",
                        help="Input is a file containing a tabulated table containing the GI of sequences to retrieve in column c. Please specify column number.")
    parser.add_argument("-p", '--append', type=str, default="",
                        help="If required, name of the blast query file that should be appended to the output fasta file")

    args = parser.parse_args()

    if args.column != "":
        accession = table_to_gi(args.accession, args.column)
    else:
        accession = list(args.accession)

    sys.stdout.write("getting protein sequences %s ... \n" % accession)
    download_seq_entrez(accession, args.multifasta, args.seq_type, args.append)
