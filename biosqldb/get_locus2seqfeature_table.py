#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


def create_locus_tag2seqfeature_table(biodb):

    import manipulate_biosqldb



    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'CREATE TABLE IF NOT EXISTS custom_tables.locus2seqfeature_id_%s (locus_tag varchar(400), seqfeature_id INT, taxon_id INT,' \
          ' orthogroup varchar(400), index orthogroup(orthogroup), ' \
          ' index locus_tag (locus_tag), index seqfeature_id(seqfeature_id), index taxon_id (taxon_id))' % biodb

    server.adaptor.execute(sql)
    server.commit()

    locus2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)

    sql = 'select locus_tag, taxon_id from orthology_detail_%s' % biodb

    locus2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql2 = 'select locus_tag, orthogroup from orthology_detail_%s' % biodb

    locus2orthogroup = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

    for locus in locus2seqfeature_id:
        try:
            sql = 'insert into custom_tables.locus2seqfeature_id_%s values ("%s", %s, %s, "%s")' % (biodb, locus,
                                                                                              locus2seqfeature_id[locus],
                                                                                              locus2taxon_id[locus],
                                                                                              locus2orthogroup[locus])
            server.adaptor.execute(sql)
        except:
            # pseudogenes
            sql = 'insert into custom_tables.locus2seqfeature_id_%s values ("%s", %s, %s, %s)' % (biodb, locus,
                                                                                              locus2seqfeature_id[locus],
                                                                                              "NULL",
                                                                                              "NULL")
            server.adaptor.execute(sql)
        server.commit()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--biodb', type=str, help="biodatabase name")


    args = parser.parse_args()
    #door_accession2door_operon_table(accession=1101)
    create_locus_tag2seqfeature_table(args.biodb)
