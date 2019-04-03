#!/usr/bin/env python

from Bio import Entrez


Entrez.email = "trestan.pillonel@unil.ch"

def cog_id2sequences(cog_id):

    import MySQLdb
    import os
    sqlpsw = os.environ['SQLPSW']
    from tempfile import NamedTemporaryFile
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="COG") # name of the data base
    cursor = conn.cursor()

    sql = 'select protein_id from cog_2014 where COG_id="%s";' % cog_id
    cursor.execute(sql, )
    prot_id_list = [str(i[0]) for i in cursor.fetchall()]
    sys.stdout.write('N sequences for cog %s: %s' % (cog_id, len(prot_id_list)))
    handle = Entrez.efetch(db="nucleotide", id=','.join(prot_id_list), rettype="fasta", retmode="text")
    seq_records = list(SeqIO.parse(handle, "fasta"))
    return seq_records


if __name__ == '__main__':
    import argparse
    import sys
    from Bio import SeqIO
    import shell_command
    parser = argparse.ArgumentParser()
    parser.add_argument("-c",'--cog_id', type=str, help="cog id")

    args = parser.parse_args()
    records = cog_id2sequences(args.cog_id)


    with open('%s.fa' % args.cog_id, 'w') as f:
        SeqIO.write(records,f, 'fasta')

    mafft_cmd = 'mafft --auto --maxiterate 1000 %s.fa > %s_mafft.fa' % (args.cog_id, args.cog_id)
    hmm_cmd = 'hmmbuild %s.hmm %s_mafft.fa' % (args.cog_id, args.cog_id)
    a, b, c = shell_command.shell_command(mafft_cmd)
    if c == 0:
        a, b, c = shell_command.shell_command(hmm_cmd)
        if c == 0:
            sys.stdout.write('\n%s ... done\n' % args.cog_id)
        else:
            sys.stdout.write('\n%s\n%s' % (a, b))
    else:
        sys.stdout.write('\n%s\n%s' % (a, b))

