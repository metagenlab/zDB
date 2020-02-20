#! /usr/bin/env python

def create_locus2old_locus_table(biodb):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)
    sql = 'create table custom_tables_seqfeature_id2old_locus_tag (seqfeature_id INT, old_locus_tag varchar(300), index seqfeature_id(seqfeature_id))' % biodb
    server.adaptor.execute(sql,)
    #locus2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)

    sql = 'select t2.seqfeature_id,t3.value from bioentry t1 ' \
          ' inner join seqfeature t2 on t1.bioentry_id=t2.bioentry_id ' \
          ' inner join seqfeature_qualifier_value as t3 on t2.seqfeature_id=t3.seqfeature_id ' \
          ' inner join term as t4 on t3.term_id=t4.term_id ' \
          ' inner join biodatabase as t5 on t5.biodatabase_id=t1.biodatabase_id ' \
          ' inner join term as t6 on t2.type_term_id=t6.term_id where t6.name="CDS" and t5.name="%s" and t4.name="old_locus_tag"' % biodb
    seq_feature_id2old_locus_tag = server.adaptor.execute_and_fetchall(sql,)
    for row in seq_feature_id2old_locus_tag:
        sql = 'insert into custom_tables_seqfeature_id2old_locus_tag values(%s,"%s")' % (biodb,
                                                                                            row[0],
                                                                                            row[1])
        server.adaptor.execute(sql,)
        server.commit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)


    args = parser.parse_args()

    create_locus2old_locus_table(args.db_name)