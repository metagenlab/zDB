#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


def create_locus_tag2seqfeature_table(biodb):

    import manipulate_biosqldb



    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'CREATE TABLE IF NOT EXISTS custom_tables.locus2seqfeature_id_%s (locus_tag varchar(400), seqfeature_id INT,' \
          ' index locus_tag (locus_tag))' % biodb
    server.adaptor.execute(sql)
    server.commit()

    locus2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)

    for locus in locus2seqfeature_id:
        sql = 'insert into custom_tables.locus2seqfeature_id_%s values ("%s", %s)' % (biodb, locus, locus2seqfeature_id[locus])
        server.adaptor.execute(sql)
        server.commit()



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--biodb', type=str, help="biodatabase name")


    args = parser.parse_args()
    #door_accession2door_operon_table(accession=1101)
    create_locus_tag2seqfeature_table(args.biodb)