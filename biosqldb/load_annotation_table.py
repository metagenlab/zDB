#! /usr/bin/env python


def import_annot(annot_file, biodb):
    import manipulate_biosqldb

    '''
    import annotation file of the form as tabulated file with 4 columns:
    1. category
    2. gene name
    3. locus tag
    4. description
    '''

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table IF NOT EXISTS custom_tables.annot_table_%s (category varchar(100), ' \
          ' gene varchar(400),' \
          ' locus_tag varchar(400),' \
          ' description TEXT,' \
          ' index locus_tag(locus_tag));' % biodb
    print sql
    server.adaptor.execute(sql,)
    with open(annot_file, 'r') as f:
        for n, one_row in enumerate(f):
            data = one_row.rstrip().split('\t')
            sql = 'insert into  custom_tables.annot_table_%s values ("%s", "%s", "%s", "%s");' % (biodb,
                                                                                              data[0],
                                                                                              data[1],
                                                                                              data[2],
                                                                                              data[3])
            print sql
            server.adaptor.execute(sql,)
        server.commit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--annot_file', type=str, help="eggnog_file (http://beta-eggnogdb.embl.de/#/app/emapper)")
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)

    args = parser.parse_args()

    import_annot(args.annot_file, args.db_name)
