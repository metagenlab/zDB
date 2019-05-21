#! /usr/bin/env python


import sys
from Bio import GenBank
from Bio import Entrez
from Bio import SeqIO
from BioSQL import BioSeqDatabase
from manipulate_biosqldb import load_db
from manipulate_biosqldb import query_yes_no








def insert_blast_table(input_file, table_name):

    from biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db('chlamydia_12_15')
    # query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
    sql_profiles_table = 'CREATE TABLE IF NOT EXISTS temp_tables.%s (' \
                         ' query_id varchar(100), ' \
                         ' subject_id varchar(100), ' \
                         ' identity float, ' \
                         ' alignment_length int, ' \
                         ' mismatches int, ' \
                         ' gap_opens int, ' \
                         ' q_start int, ' \
                         ' q_end int, ' \
                         ' s_start int, ' \
                         ' s_end int, ' \
                         ' evalue float, ' \
                         ' bit_score float)' % (table_name)
    try:
        print sql_profiles_table
        server.adaptor.execute(sql_profiles_table)
    except:
        print 'problem creating the sql table'

    import shell_command
    import os
    wd = os.getcwd()
    path = os.path.join(wd, input_file)

    sqlpsw = os.environ['SQLPSW']
    cmd = 'mysql -uroot -p%s biosqldb -e \'LOAD DATA LOCAL INFILE "%s" INTO TABLE temp_tables.%s;\'' % (sqlpsw,path, table_name)
    print cmd
    a, b, c = shell_command.shell_command(cmd)
    print a
    print b
    print c






if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--table_name', type=str, help="sql table name", required=True)
    parser.add_argument("-i", '--input_table', type=str, help="file name", required=True)

    args = parser.parse_args()

    insert_blast_table(args.input_table, args.table_name)


