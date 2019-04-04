#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


def create_locus_tag2seqfeature_table(biodb, locus2seqfeature_id=False, locus2taxon_id=False):

    from biosqldb import manipulate_biosqldb



    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'CREATE TABLE IF NOT EXISTS custom_tables.locus2seqfeature_id_%s (locus_tag varchar(400), ' \
          ' seqfeature_id INT, ' \
          ' taxon_id INT,' \
          ' index locus_tag (locus_tag), ' \
          ' index seqfeature_id(seqfeature_id), ' \
          ' index taxon_id (taxon_id))' % biodb

    server.adaptor.execute(sql)
    server.commit()
    if not locus2seqfeature_id:
        locus2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)
    if not locus2taxon_id:
        locus2taxon_id = manipulate_biosqldb.locus_tag2genome_taxon_id(server, biodb)

    for locus in locus2seqfeature_id:
        try:
            sql = 'insert into custom_tables.locus2seqfeature_id_%s values ("%s", %s, %s)' % (biodb, locus,
                                                                                              locus2seqfeature_id[locus],
                                                                                              locus2taxon_id[locus])
            server.adaptor.execute(sql)
        except:
            # pseudogenes
            sql = 'insert into custom_tables.locus2seqfeature_id_%s values ("%s", %s, %s)' % (biodb,
                                                                                              locus,
                                                                                              locus2seqfeature_id[locus],
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
