#! /usr/bin/env python


def import_annot(annot_file, biodb):
    from chlamdb.biosqldb import manipulate_biosqldb
    from datetime import datetime

    '''
    import annotation file of the form as tabulated file with 4 columns:
    1. category
    2. gene name
    3. locus tag
    4. description
    '''

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table IF NOT EXISTS custom_tables_annot_table (category varchar(100), ' \
          ' gene varchar(400),' \
          ' locus_tag varchar(400),' \
          ' description TEXT,' \
          ' reference TEXT,' \
          ' date varchar(400),' \
          ' index locus_tag(locus_tag));' % biodb
    print sql
    server.adaptor.execute(sql,)
    with open(annot_file, 'r') as f:
        for n, one_row in enumerate(f):
            now = datetime.now()
            str_date = "%s-%s-%s" % (now.year, now.month, now.day)

            data = one_row.rstrip().split('\t')
            print data
            print len(data)
            if len(data) == 5:
                sql = 'insert into  custom_tables_annot_table (category, gene, locus_tag, description,' \
                      ' reference, date) values ("%s", "%s", "%s", "%s", "%s", "%s");' % (biodb,
                                                                                                  data[0],
                                                                                                  data[1],
                                                                                                  data[2],
                                                                                                  data[3],
                                                                                                  data[4],
                                                                                                  str_date)
                #print sql
                server.adaptor.execute(sql,)
        server.commit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--annot_file', type=str, help="eggnog_file (http://beta-eggnogdb.embl.de/#/app/emapper)")
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)

    args = parser.parse_args()

    import_annot(args.annot_file, args.db_name)
