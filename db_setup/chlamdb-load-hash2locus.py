#!/usr/bin/env python


def load_hash2locus(hash2locus_list, biodb):
    import MySQLdb
    import os
    from chlamdb.biosqldb import manipulate_biosqldb
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                           user="root", # your username
                           passwd=sqlpsw) # name of the data base
    cursor = conn.cursor()

    sql = 'select locus_tag, seqfeature_id from annotation.seqfeature_id2locus_%s' % biodb
    cursor.execute(sql,)
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    sql = 'create table annotation.hash2seqfeature_id_%s (hash varchar(300), seqfeature_id INTEGER)' % biodb
    cursor.execute(sql,)

    for hash in hash2locus_list:
        locus_list = hash2locus_list[hash]
        for locus in locus_list:
            cursor.execute('insert into annotation.hash2seqfeature_id_%s values ("%s", "%s")' % (biodb, hash, locus_tag2seqfeature_id[locus]))
    conn.commit()
    sql1 = 'create index h1 ON annotation.hash2seqfeature_id_%s(hash)' % biodb
    sql2 = 'create index h2 ON annotation.hash2seqfeature_id_%s(seqfeature_id)' % biodb
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
