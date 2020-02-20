#! /usr/bin/env python


def index_silix_families(family_list, biodb, silix_cutoff):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    sql = 'create table comparative_tables.silix_families_%s_%s (silix_id INT primary key, silix_name)' % (biodb, silix_cutoff)




def import_silix(silix_output, biodb, silix_cutoff):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)
    import re
    # columns
    # 0 taxon_id
    # 1 seqfeature_id
    # 2 silix_family_id

    sql = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % biodb
    locus2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    sql = 'select locus_tag, taxon_id from orthology_detail' % biodb
    locus2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    sql = 'create table IF NOT EXISTS comparative_tables_silix_%s (taxon_id INT, ' \
          ' seqfeature_id INT, ' \
          ' silix_id INT,' \
          ' INDEX seqfeature_id (seqfeature_id), index taxon_id(taxon_id), index silix_id(silix_id));' % (biodb, silix_cutoff)
    server.adaptor.execute(sql,)
    server.commit()
    with open(silix_output, 'r') as f:
        row_list = [i for i in f]
        family_list = list(set([i.rstrip().split('\t')[0] for i in row_list]))

        for n, row in enumerate(row_list):
            if n == 0:
                continue
            data = row.rstrip().split('\t')

            family_index = family_list.index(data[0])
            locus_tag = data[1]
            taxon_id = locus2taxon_id[locus_tag]
            seqfeature_id = locus2seqfeature_id[locus_tag]
            print family_index, locus_tag

            sql = 'insert into comparative_tables_silix_%s values (%s,  %s, %s);' % (biodb,
                                                                 silix_cutoff,
                                                                 taxon_id,
                                                                 seqfeature_id,
                                                                 family_index)

            server.adaptor.execute(sql,)

        server.commit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", '--silix_output', type=str, help="silix_output")
    parser.add_argument("-d", '--biodb', type=str, help="biodb")
    parser.add_argument("-c", '--silix_cutoff', type=str, help="silix identity cutoff", required=True)

    args = parser.parse_args()

    import_silix(args.silix_output, args.biodb, args.silix_cutoff)


