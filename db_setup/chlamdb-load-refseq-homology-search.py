#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# Load blast/plast/diamond results into biosql databas
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: april 2019
# ---------------------------------------------------------------------------


def load_blastnr_file_into_db(locus_tag2taxon_id,
                              locus_tag2seqfeature_id,
                              protein_id2seqfeature_id,
                              locus_tag2bioentry_id,
                              biodb,
                              hash2locus_list,
                              linear_taxonomy,
                              diamond_refseq,
                              refseq_taxonomy):

    '''
    Load tabulated blast results into sql table blastnr_`db_name`
    Ab unique identifier (primary sql key) is attributed to each blast hsp

    :param seqfeature_id2locus_tag: dictionnary with whole data for `biodb` biodatabase
    :param locus_tag2seqfeature_id: dictionnary with whole data for `biodb` biodatabase
    :param protein_id2seqfeature_id: dictionnary with whole data for `biodb` biodatabase
    :param db_name: name of the biodatabase
    :param input_blast_files: all input tabulated blast files
    :return: None
    '''


    import time
    import sqlite3
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    mysql_conn = server.adaptor.conn
    mysql_cursor = server.adaptor.cursor

    sqlite3_conn = sqlite3.connect(diamond_refseq)
    sqlite3_cursor = sqlite3_conn.cursor()

    sql1 = 'attach "%s" as linear_taxonomy' % linear_taxonomy
    sql2 = 'attach "%s" as refseq_taxonomy' % refseq_taxonomy

    sqlite3_cursor.execute(sql1)
    sqlite3_cursor.execute(sql2)

    sql3 = 'select t1.*,t2.*,t3.superkingdom from diamond_refseq t1 ' \
           ' inner join refseq_taxonomy_refseq_hits t2 on t1.sseqid=t2.accession ' \
           'inner join linear_taxonomy_ncbi_taxonomy t3 on t2.taxid=t3.tax_id'

    n = 0
    for row in sqlite3_cursor.execute(sql3):

        '''
        0 hit_count
        1 qseqid
        2 sseqid
        3 pident
        4 length
        5 mismatch
        6 gapopen
        7 qstart
        8 qend
        9 sstart
        10 send
        11 evalue
        12 bitscore
        13 accession
        14 taxid
        15 description
        16 length
        17 superkingdom
        '''

        if n % 10000 == 0:
            print(time.ctime() + ': %s...' % n)
            mysql_conn.commit()

        hit_n = row[0]
        query_hash = row[1]
        evalue = row[11]
        percent_identity = float(row[3])
        gaps = int(row[6])
        align_length = int(row[4])
        query_start = int(row[7])
        query_end = int(row[8])
        subject_start = int(row[9])
        subject_end = int(row[10])
        bit_score = float(row[12])
        subject_accession = row[2]
        subject_taxon_id = row[14]
        subject_kingdom = row[17]
        hit_length = row[16]

        subject_description_data = re.findall("(.+) \[(.*)\]", row[15])
        try:
            subject_scientific_names = subject_description_data[0][1]
            subject_title = subject_description_data[0][0]
        except:
            subject_scientific_names = "?"
            subject_title = row[15]


        '''
        1 query_taxon_id int
        2 query_bioentry_id INT
        3 seqfeature_id INT
        4 hit_number int
        5 subject_gi int
        6 subject_accession varchar(200)
        7 subject_kingdom varchar(200)
        8 subject_scientific_name TEXT(2000000)
        9 subject_taxid INT
        10 subject_title VARCHAR(2000)
        11 evalue varchar(200)
        12 bit_score float
        12 percent_identity float
        14 gaps int
        15 length int
        16 query_start int
        17 query_end int
        18 subject_start int
        19 subject_end
        '''

        locus_list = hash2locus_list[query_hash]

        n += 1
        for locus_tag in locus_list:
            seqfeature_id = locus_tag2seqfeature_id[locus_tag]

            values = '(%s,%s,%s,%s,"%s","%s","%s",%s,"%s","%s",%s,%s,%s,%s,%s,%s,%s,%s, %s);' % (locus_tag2taxon_id[locus_tag],
                                                                                                 locus_tag2bioentry_id[locus_tag],
                                                                                                 seqfeature_id,
                                                                                                 hit_n,
                                                                                                 subject_accession,
                                                                                                 subject_kingdom,
                                                                                                 subject_scientific_names,
                                                                                                 subject_taxon_id,
                                                                                                 subject_title,
                                                                                                 evalue,
                                                                                                 bit_score,
                                                                                                 percent_identity,
                                                                                                 gaps,
                                                                                                 align_length,
                                                                                                 query_start,
                                                                                                 query_end,
                                                                                                 subject_start,
                                                                                                 subject_end,
                                                                                                 hit_length
                                                                                                 )

            sql = 'insert into blastnr_%s values %s' % (biodb,
                                                        values)
            mysql_cursor.execute(sql,)


def create_sql_plastnr_tables(db_name):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor
    
    sql_plast = 'CREATE TABLE IF NOT EXISTS blastnr_%s (query_taxon_id INT,' \
                ' query_bioentry_id INT,' \
                ' seqfeature_id INT,' \
                ' hit_number int,' \
                ' subject_accession varchar(200),' \
                ' subject_kingdom varchar(200),' \
                ' subject_scientific_name TEXT(2000000), ' \
                ' subject_taxid INT,' \
                ' subject_title VARCHAR(2000),' \
                ' evalue varchar(200),' \
                ' bit_score float,' \
                ' percent_identity float,' \
                ' gaps int,' \
                ' length int,' \
                ' query_start int,' \
                ' query_end int,' \
                ' subject_start int,' \
                ' subject_end int,' \
                ' subject_length int,' \
                ' INDEX query_taxon_id (query_taxon_id),' \
                ' INDEX query_bioentry_id (query_bioentry_id),' \
                ' INDEX hit_number (hit_number),' \
                ' INDEX seqfeature_id (seqfeature_id),' \
                ' INDEX subject_taxid(subject_taxid))' % (db_name)

    try:

        cursor.execute(sql_plast)
        conn.commit()
    except:
        print(sql_plast)
        print('not created')


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    import sys
    import os
    import json
    import numpy
    import re
    import time
    from datetime import datetime
    from multiprocessing import Process
    import  chlamdb_setup_utils

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--mysql_database', type=str, help="Biosql biodatabase name")
    parser.add_argument("-t", '--create_tables', action='store_true', help="Create SQL tables")
    parser.add_argument("-f", '--filter_n_hits', action='store_true', help="filter_n_hits (max 100 hits/locus)")
    parser.add_argument("-u", '--hash2locus_tag', type=str, help="Tab separated file with correspondance between sequence hashes and locus tags")
    parser.add_argument("-lt", '--linear_taxonomy', type=str, help="linear_taxonomy.db")
    parser.add_argument("-rd", '--diamond_refseq', type=str, help="diamond_refseq.db")
    parser.add_argument("-rt", '--refseq_taxonomy', type=str, help="refseq_taxonomy.db")

    args = parser.parse_args()

    mysql_host = 'localhost'
    mysql_user = 'root'

    mysql_pwd = os.environ['SQLPSW']
    mysql_db = 'blastnr'

    db_name = args.mysql_database

    biodb = args.mysql_database

    server, db = manipulate_biosqldb.load_db(biodb)

    hash2locus_list = chlamdb_setup_utils.get_hash2locus_list(args.hash2locus_tag)

    sys.stdout.write("creating locus_tag2seqfeature_id")
    locus_tag2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)

    sys.stdout.write("creating protein_id2seqfeature_id")
    protein_id2seqfeature_id = manipulate_biosqldb.protein_id2seqfeature_id_dict(server, biodb)

    create_sql_plastnr_tables(db_name)

    print('get locus2taxon_id')
    sql = 'select locus_tag, taxon_id from orthology_detail'
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    print('get locus2bioentry')
    sql2 = 'select locus_tag,bioentry_id from biodatabase t1 ' \
           ' inner join bioentry as t2 on t1.biodatabase_id=t2.biodatabase_id' \
           ' inner join orthology_detail t3 on t2.accession=t3.accession where t1.name="%s"' % (db_name)

    locus_tag2bioentry_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

    load_blastnr_file_into_db(locus_tag2taxon_id,
                              locus_tag2seqfeature_id,
                              protein_id2seqfeature_id,
                              locus_tag2bioentry_id,
                              db_name,
                              hash2locus_list,
                              args.linear_taxonomy,
                              args.diamond_refseq,
                              args.refseq_taxonomy)
