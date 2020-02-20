#!/usr/bin/python

from chlamdb.biosqldb import manipulate_biosqldb

def locus2inference_table(biodb):

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'CREATE TABLE locus_tag2uniprot_hit (locus_tag varchar(400),' \
          ' uniprot_id varchar(400), index locus_tag(locus_tag))' % biodb

    server.adaptor.execute(sql,)

    locus2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)

    for locus in locus2seqfeature_id:
        sql = 'select value from seqfeature_qualifier_value where seqfeature_id=%s and value like "%%%%UniProtKB%%%%"' % (locus2seqfeature_id[locus])
        try:
            data = server.adaptor.execute_and_fetchall(sql,)[0][0]
            sql2 = 'insert into locus_tag2uniprot_hit values ("%s", "%s")' %(biodb, locus, data.split(':')[2])
            try:
                server.adaptor.execute(sql2,)
                server.commit()
            except:
                print sql2
        except:
            pass
locus2inference_table("chlamydia_04_16")