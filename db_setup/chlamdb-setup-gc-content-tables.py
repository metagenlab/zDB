#! /usr/bin/env python


def biodb2cds_gc(biodb):
    from chlamdb.biosqldb import manipulate_biosqldb
    from Bio.SeqUtils import GC123

    server, db = manipulate_biosqldb.load_db(biodb)

    sql1 = 'select distinct accession from orthology_detail'
    sql2 = 'select locus_tag, taxon_id from orthology_detail'
    sql3 = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id'

    accession_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql1,)]

    locus2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
    locus2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3,))

    sql_head = 'create table IF NOT EXISTS custom_tables_gc_content (taxon_id INT, ' \
          ' seqfeature_id INT,' \
          ' seq_length INT, gc_percent FLOAT, gc_1 FLOAT, gc_2 FLOAT, gc_3 FLOAT, ' \
               ' INDEX seqfeature_id(seqfeature_id), index taxon_id(taxon_id))'

    server.adaptor.execute(sql_head,)

    count_all=0
    for accession in accession_list:
        print (accession)
        record = db.lookup(accession=accession)
        seq = record.seq
        for n, feature in enumerate(record.features):
            if feature.type == 'CDS' and not 'pseudo' in feature.qualifiers and not 'pseudogene' in feature.qualifiers and 'translation' in feature.qualifiers:
                count_all+=1
                dna_sequence = feature.extract(seq)
                locus = feature.qualifiers['locus_tag'][0]

                gc, gc1, gc2, gc3 = GC123(str(dna_sequence))
                sql = 'insert into  custom_tables_gc_content values (%s, %s, %s, %s, %s, %s, %s);' % (locus2taxon_id[locus],
                                                                                                      locus2seqfeature_id[locus],
                                                                                                      len(dna_sequence),
                                                                                                      round(gc,2),
                                                                                                      round(gc1),
                                                                                                      round(gc2),
                                                                                                      round(gc3))
                server.adaptor.execute(sql,)
        server.commit()

if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)

    args = parser.parse_args()

    biodb2cds_gc(args.db_name)
    
    manipulate_biosqldb.update_config_table(args.database_name, "GC_statistics")
