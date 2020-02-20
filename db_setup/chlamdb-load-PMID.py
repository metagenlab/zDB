#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2019
# ---------------------------------------------------------------------------


def load_PMID(PMID_db_path, hash2locus_list, db_name):

    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.biosqldb import biosql_own_sql_tables
    import sqlite3 
    
    sqlite_conn = sqlite3.connect(PMID_db_path)
    sqlite_cursor = sqlite_conn.cursor()
    
    server, db = manipulate_biosqldb.load_db(db_name)
    
    sql1 = 'create table string_seqfeature_id2pmid (seqfeature_id INTEGER, pmid INTEGER);'
    sql2 = 'create table if not exists string_pmid2data (pmid INTEGER, title TEXT, authors TEXT, source TEXT, abstract TEXT);'
    
    server.adaptor.execute(sql1,)
    server.adaptor.execute(sql2,)
    
    # get locus_tag2seqfeaute_id dictionnary
    sql = 'select locus_tag, seqfeature_id from annotation_seqfeature_id2locus'
    locus_tag2seqfeaure_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    # retrieve existing pmid data (if any)
    sql = 'select pmid from string.pmid2data '
    pmid_already_in_db = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql,)]
    
    sql = 'select * from hash2pmid'
    hash_and_pmid = sqlite_cursor.execute(sql,).fetchall()
    
    sql = 'select * from pmid2data'
    pmid2data = manipulate_biosqldb.to_dict(sqlite_cursor.execute(sql,).fetchall())
    
    template1 = 'insert into string_pmid2data values (%s, %s, %s, %s, %s)'
    template2 = 'insert into string_seqfeature_id2pmid values(%s, %s)'
    for n, row in enumerate(hash_and_pmid):
        hash, pmid = row 
        locus_list = hash2locus_list[hash]
        if pmid not in pmid_already_in_db:
            try:
                server.adaptor.execute(template1, [pmid] + list(pmid2data[str(pmid)]))
                pmid_already_in_db.append(pmid)
            except KeyError:
                print("Problem with pmid %s" % pmid)
        for locus_tag in locus_list:
            server.adaptor.execute(template2,[locus_tag2seqfeaure_id[locus_tag], pmid])
    server.commit()

    sql1 = 'create index sf on string_seqfeature_id2pmid (seqfeature_id)'
    sql2 = 'create index pm1 on string_seqfeature_id2pmid (pmid)'
    sql3 = 'create index pm2 on string_pmid2data (pmid)'
    
    server.adaptor.execute(sql1,)
    server.adaptor.execute(sql2,)
    server.adaptor.execute(sql3,)

if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    import chlamdb_setup_utils

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--PMID_sqlite_db', type=str, help="input PMID sqlite3 db file")
    parser.add_argument("-d", '--database_name', type=str, help="database name")
    parser.add_argument("-c", '--corresp_table', type=str, help="hash to locus correspondance table")

    args = parser.parse_args()

    hash2locus_list = chlamdb_setup_utils.get_hash2locus_list(args.corresp_table)

    load_PMID(args.PMID_sqlite_db,
              hash2locus_list,
              args.database_name)
