#!/usr/bin/env python


def load_hash2locus(hash2locus_list, 
                    biodb):
    import MySQLdb
    import os
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    sql = 'select locus_tag, seqfeature_id from annotation_seqfeature_id2locus'
    cursor.execute(sql,)
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    sql = 'create table annotation_hash2seqfeature_id (hash varchar(300), seqfeature_id INTEGER)'
    cursor.execute(sql,)

    for hash in hash2locus_list:
        locus_list = hash2locus_list[hash]
        for locus in locus_list:
            cursor.execute('insert into annotation_hash2seqfeature_id values ("%s", "%s")' % (hash, 
                                                                                              locus_tag2seqfeature_id[locus]))
    conn.commit()
    sql1 = 'create index h1 ON annotation_hash2seqfeature_id(hash)'
    sql2 = 'create index h2 ON annotation_hash2seqfeature_id(seqfeature_id)'
    cursor.execute(sql1,)
    cursor.execute(sql2,)
    conn.commit()



if __name__ == '__main__':
    import argparse
    import chlamdb_setup_utils
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", '--hash2locus_tag', type=str, help="Tab separated file with correspondance between sequence hashes and locus tags")
    parser.add_argument("-d", '--database_name', type=str, help="database name", default=False)

    args = parser.parse_args()

    hash2locus_list = chlamdb_setup_utils.get_hash2locus_list(args.hash2locus_tag)

    load_hash2locus(hash2locus_list, args.database_name)
