#!/usr/bin/env python

# blast a fasta file against a user selected database over internet (NCBI blast) or a local database
# Author: Claire Bertelli (claire.bertelli[@]gmail.com)
# Date: 04.2015
# ---------------------------------------------------------------------------

import argparse
from Bio import Entrez
from shell_command import shell_command
import re

Entrez.email = "claire.bertelli@chuv.ch"


def blast_fasta(blast_flavor, input, blastdb, evalue, nb_hit, local):
    '''
    :parameters for BLAST
    :return: blast results with taxonomic information
    scientific names and kingdom will only be retrieved if the taxid database has been installed locally ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
    '''
    output_name = re.sub(r"([a-zA-Z_0-9]+)\.([a-zA-Z]+)", "blast\\1.tab", input)
    output_format = "6 qgi qacc sgi sacc sscinames sskingdoms staxids evalue "\
                    "nident pident positive gaps length qstart qend qcovs "\
                    "sstart send sstrand stitle"
    if local:
        cmd = "%s -task %s -query %s -out %s -db %s -evalue %s -max_target_seqs %s -outfmt '%s'" % (
            blast_flavor, blast_flavor, input, output_name, blastdb, evalue, nb_hit, output_format)
        print cmd
        out, err, code = shell_command(cmd)
        if code != 0:
            print err
    else:
        cmd = "%s -remote -task %s -query %s -out %s -db %s -evalue %s -max_target_seqs %s -outfmt '%s'" % (
            blast_flavor, blast_flavor, input, output_name, blastdb, evalue, nb_hit, output_format)
        out, err, code = shell_command(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", '--flavor', type=str,
        help="BLAST flavor: blastp for protein-protein, blastn for "
             "nucleotide-nucleotide, blastx for translated nucleotide query "
             "against protein db, tblastn for protein query against translated db")
    parser.add_argument("-i", '--input', type=str,
                        help="Fasta or multifasta sequence")
    parser.add_argument("-d", '--blastdb', type=str,
                        help="BLAST database to use")
    parser.add_argument("-e", '--evalue', default=0.00001, type=str,
                        help="E-value cutoff to report results (default=10-5)")
    parser.add_argument("-n", '--nb_hit', default=10, type=int,
                        help="Number of hits to keep (default=10)")
    parser.add_argument("-l", '--local', action="store_true",
                        help="Run BLAST locally")

    args = parser.parse_args()

    print "Performing blast..."
    blast_fasta(args.flavor, args.input, args.blastdb,
                args.evalue, args.nb_hit, local=args.local)
    print "BLAST result now available"
