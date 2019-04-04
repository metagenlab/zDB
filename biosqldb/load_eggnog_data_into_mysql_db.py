#! /usr/bin/env python


import sys
from Bio import GenBank
from Bio import Entrez
from Bio import SeqIO
from BioSQL import BioSeqDatabase
from manipulate_biosqldb import load_db
from manipulate_biosqldb import query_yes_no

def import_eggnog(eggnog_file, biodb):
    from biosqldb import manipulate_biosqldb
    import biosql_own_sql_tables
    import ete2
    import re
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table IF NOT EXISTS COG.eggnog_hits_%s (query_seq varchar(100), ' \
          ' best_hit_eggNOG_ortholog varchar(400),' \
          ' best_hit_evalue FLOAT,' \
          ' best_hit_score FLOAT,' \
          ' predicted_name varchar(100),' \
          ' index query_seq(query_seq));' % biodb
    print sql
    server.adaptor.execute(sql,)
    with open(eggnog_file, 'r') as f:
        for n, one_hit in enumerate(f):
            if n == 0:
                continue
            else:
                data = one_hit.rstrip().split('\t')
                if data[1] == '-':
                    continue
                else:
                    sql = 'insert into  COG.eggnog_hits_%s values ("%s", "%s", %s, %s, "%s");' % (biodb, data[0], data[1], data[2], data[3], data[4])
                    print sql
                    server.adaptor.execute(sql,)
        server.commit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--eggnog_file', type=str, help="eggnog_file (http://beta-eggnogdb.embl.de/#/app/emapper)")
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)

    args = parser.parse_args()

    import_eggnog(args.eggnog_file, args.db_name)


