#!/usr/bin/python

def get_pathway_ko_association_table():
    import os
    import MySQLdb
    import manipulate_biosqldb
    import urllib2
    from Bio.KEGG.KGML import KGML_parser
    import re

    sqlpsw = os.environ['SQLPSW']

    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="enzyme") # name of the data base
    cursor = conn.cursor()

    sql = 'create table enzyme.pathway2ortholog_associations (pathway_id INT, node_id INT, ko_id varchar(200), ' \
          ' index pathway_id(pathway_id), index node_id(node_id), index ko_id(ko_id));'
    #cursor.execute(sql,)
    #conn.commit()

    sql2 = 'select pathway_name,pathway_id from enzyme.kegg_pathway'
    cursor.execute(sql2,)

    pathway2pathway_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    for pathway in pathway2pathway_id:
        print pathway

        url_template = 'http://rest.kegg.jp/get/%s/kgml' % re.sub('map', 'ko', pathway)
        print url_template
        try:
            f = urllib2.urlopen(url_template)
        except:
            continue
        from Bio.Graphics import KGML_vis


        pathway_KGML = KGML_parser.read(f.read())

        # Loop over the orthologs in the pathway, and change the
        # background colour
        orthologs = [e for e in pathway_KGML.orthologs]
        for o in orthologs:
            #print dir(o)
            #print o.id
            ko_temp_list = list(set([i.rstrip() for i in o.name.split('ko:')]))
            ko_temp_list = filter(None, ko_temp_list)
            for ko in ko_temp_list:
                sql = 'insert into enzyme.pathway2ortholog_associations values(%s, %s, "%s")' % (pathway2pathway_id[pathway],
                                                                                               o.id,
                                                                                               ko)
                cursor.execute(sql,)
        conn.commit()

get_pathway_ko_association_table()
