#!/usr/bin/env python


def get_interpro_entry_tables(interpro_release,
                              biodb):
    '''
    :param interpro_release: par exemple: 60.0
    :return:
    '''
    import urllib
    from chlamdb.biosqldb import manipulate_biosqldb
    
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    sql = 'CREATE table if not exists interpro_entry (interpro_id INTEGER PRIMARY KEY, name varchar(400), description TEXT)'

    cursor.execute(sql,)
    conn.commit()
    link = 'ftp://ftp.ebi.ac.uk/pub/databases/interpro/%s/entry.list' % interpro_release
    print(link)
    entry_list = urllib.request.urlopen(link).read().decode('utf-8').split("\n")
    for i, line in enumerate(entry_list):
        print(i)
        if not 'IPR' in line:
            continue
        # ENTRY_AC	ENTRY_TYPE	ENTRY_NAME
        # IPR000126	Active_site	Serine proteases, V8 family, serine active site
        
        data = line.rstrip().split("\t")
        accession = data[0]
        description = data[2]
        sql = 'insert into interpro_entry (name, description) values ("%s", "%s")' % (accession, description)
        cursor.execute(sql,)
    conn.commit()


def get_interpro_parent2child_table(interpro_release):
    link = ''


def get_interpro2go_table(biodb):

    '''

    :param interpro_release: par exemple: 60.0
    :return:
    '''
    import urllib2
    from chlamdb.biosqldb import manipulate_biosqldb
    
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor


    sql = 'CREATE table interpro2gene_ontology (interpro_id INT, go_id INT, index interpro_id(interpro_id), index go_id(go_id))'
    #cursor.execute(sql,)
    #conn.commit()

    sql = 'select name,interpro_id from entry'
    cursor.execute(sql,)
    interpro_name2interpro_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    sql = 'select acc,id from gene_ontology.term where acc like "GO%"'
    cursor.execute(sql,)
    go_name2go_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    link = 'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro2go'
    print(link)
    req = urllib2.Request(link)
    entry_list = urllib2.urlopen(req)
    echec = 0
    for i, line in enumerate(entry_list):
        if i % 1000 == 0:
            print(i)
        if line[0] == '!':
            continue
        interpro_name = line.rstrip().split(':')[1].split(' ')[0]
        go_name = line.rstrip().split(' ; ')[1]
        try:
            sql = 'insert into interpro2gene_ontology values(%s, %s)' % (interpro_name2interpro_id[interpro_name], go_name2go_id[go_name])
            cursor.execute(sql,)
            conn.commit()
        except:
            echec+=1
    print (echec, i)


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--db_name', type=str, help='DB name', default='Biodb name')
    parser.add_argument('-v', '--interpro_version', type=str, help='Interpro version')

    args = parser.parse_args()

    get_interpro_entry_tables(args.interpro_version, args.db_name)
    #get_interpro2go_table(args.db_name)