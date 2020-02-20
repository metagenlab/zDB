#! /usr/bin/env python


def import_taxonomy(annot_file, biodb):
    from chlamdb.biosqldb import manipulate_biosqldb

    '''
    import manually established taxonomy

    0 accession
    1 phylum_id
    2 phylum
    3 subphylum_id
    4 subphylum
    5 order_id
    6 order
    7 family_id
    8 family
    9 genus_id
    10 genus
    11 species_id
    12 species_description
    13 description
    '''

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select accession, taxon_id from bioentry t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id' \
          ' where t2.name="%s"' % biodb

    accession2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'create table IF NOT EXISTS custom_tables_taxonomy (taxon_id INT, ' \
          ' phylum_id INT,' \
          ' phylum_description TEXT,' \
          ' subphylum_id INT,' \
          ' subphylum_description TEXT,' \
          ' order_id INT,' \
          ' order_description TEXT,' \
          ' family_id INT,' \
          ' family_description TEXT,' \
          ' genus_id INT,' \
          ' genus_description TEXT,' \
          ' species_id INT,' \
          ' species_description TEXT);' % biodb
    print sql
    server.adaptor.execute(sql,)
    with open(annot_file, 'r') as f:
        for n, one_row in enumerate(f):
            if n==0:
                continue
            data = one_row.rstrip().split('\t')
            try:
                taxon_id = int(data[0])
            except:
                try:
                    taxon_id = accession2taxon_id[data[0]]
                except:
                    IOError('invalid accssion/taxon_id')
            sql = 'insert into  custom_tables_taxonomy values (%s, %s, "%s", %s,"%s", %s, "%s", %s, "%s", %s, "%s", %s, "%s");' % (biodb,
                                                                                                                      taxon_id,
                                                                                                                      data[1],
                                                                                                                      data[2],
                                                                                                                      data[3],
                                                                                                                      data[4],
                                                                                                                      data[5],
                                                                                                                      data[6],
                                                                                                                        data[7],
                                                                                                                        data[8],
                                                                                                                        data[9],
                                                                                                                        data[10],
                                                                                                                        data[11],
                                                                                                                        data[12])
            print sql
            server.adaptor.execute(sql,)
        server.commit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--taxonomy_file', type=str, help="taxonomy tab file (see format in fct)")
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)

    args = parser.parse_args()

    import_taxonomy(args.taxonomy_file, args.db_name)
