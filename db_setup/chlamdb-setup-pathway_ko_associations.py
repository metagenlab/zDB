#!/usr/bin/env python

def get_pathway_ko_association_table(biodb):
    from chlamdb.biosqldb import manipulate_biosqldb
    import urllib
    from Bio.KEGG.KGML import KGML_parser
    import re
    from chlamdb.biosqldb import manipulate_biosqldb
    
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    sql = 'create table enzyme_pathway2ortholog_associations (pathway_id INT, node_id INT, ko_id varchar(200), ' \
          ' index pathway_id(pathway_id), index node_id(node_id), index ko_id(ko_id));'
    cursor.execute(sql,)
    conn.commit()

    sql2 = 'select pathway_name,pathway_id from enzyme_kegg_pathway'
    cursor.execute(sql2,)

    pathway2pathway_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    for pathway in pathway2pathway_id:
        print(pathway)

        url_template = 'http://rest.kegg.jp/get/%s/kgml' % re.sub('map', 'ko', pathway)
        print(url_template)
        try:
            f = urllib.request.urlopen(url_template)
        except:
            continue
        from Bio.Graphics import KGML_vis

        pathway_KGML = KGML_parser.read(f.read().decode("UTF-8"))

        # Loop over the orthologs in the pathway, and change the
        # background colour
        orthologs = [e for e in pathway_KGML.orthologs]
        for o in orthologs:
            ko_temp_list = list(set([i.rstrip() for i in o.name.split('ko:')]))
            ko_temp_list = filter(None, ko_temp_list)
            for ko in ko_temp_list:
                sql = 'insert into enzyme_pathway2ortholog_associations values(%s, %s, "%s")' % (pathway2pathway_id[pathway],
                                                                                                 o.id,
                                                                                                 ko)
                cursor.execute(sql,)
        conn.commit()

get_pathway_ko_association_table(biodb)
