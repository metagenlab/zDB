#!/usr/bin/env python

import sqlite3 
from chlamdb.biosqldb import manipulate_biosqldb

def load_paperblast_snipperts(paperblast_sqlite, db_name):
    paperblast_conn = sqlite3.connect(paperblast_sqlite)
    paperblast_cursor = paperblast_conn.cursor()

    server, db = manipulate_biosqldb.load_db(db_name)

    sql1 = 'create table if not exists string.snippets (snipped_id INT AUTO_INCREMENT PRIMARY KEY, pmid INT, accession varchar(200), query_term varchar(200), snippet TEXT)'

    server.adaptor.execute(sql1,)

    sql_paperblast_hits = 'select distinct accession,pmid,id from string.paperblast_entry t1 ' \
                          ' inner join string.paperblast2pmid t2 on t1.id=t2.paperblast_id;'

    match_Snippet = 0
    match_GenePaper = 0
    match_both = 0
    no_match = 0

    data = server.adaptor.execute_and_fetchall(sql_paperblast_hits,)

    insert_template = 'insert into string.snippets (pmid,accession,query_term,snippet) values (%s,%s,%s,%s)'
    
    snipped_id_insert_template = 'update string.paperblast2pmid set snippet_id=%s where paperblast_id=%s and PMID=%s'

    for n, row in enumerate(data):
        snippet_id = False
        if n % 1000 == 0:
            server.adaptor.commit()
            print(n)
        sql = 'select pmId, geneId, queryTerm,snippet from Snippet where geneId="%s" and pmId="%s"' % (row[0], row[1])

        sql2 = 'select pmId, geneId, queryTerm from GenePaper where geneId="%s" and pmId="%s"' % (row[0], row[1])

        snippet = paperblast_cursor.execute(sql,).fetchall()

        if len(snippet) != 0:
            pmid = snippet[0][0]
            accession = snippet[0][1]
            query_term = snippet[0][2]
            snippet = snippet[0][3]

            server.adaptor.execute(insert_template, 
                                  [pmid, accession, query_term, snippet])
            snippet_id = server.adaptor.last_id("string.snippets")
        else:
            GenePaper = paperblast_cursor.execute(sql2,).fetchall()
            if len(GenePaper) > 0:
                pmid = GenePaper[0][0]
                accession = GenePaper[0][1]
                query_term = GenePaper[0][2]
                snippet = '-'

                server.adaptor.execute(insert_template, 
                                    [pmid, accession, query_term, snippet])
                
                snippet_id = server.adaptor.last_id("string.snippets")

        if snippet_id:
            server.adaptor.execute(snipped_id_insert_template,
                                  [snippet_id, row[2], pmid])

    server.adaptor.commit()

if __name__ == '__main__':
    import argparse
    import chlamdb_setup_utils
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)
    parser.add_argument("-p", '--paperblast_db_name', type=str, help="paperblast db_name", required=True)

    args = parser.parse_args()

    load_paperblast_snipperts(args.paperblast_db_name, 
                              args.db_name)
    