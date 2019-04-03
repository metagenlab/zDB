#! /usr/bin/env python


def import_checkm(checkm_file, biodb):
    import manipulate_biosqldb
    import re

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select accession, taxon_id from bioentry t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id' \
          ' where t2.name="%s"' % biodb

    accession2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

#   Bin Id             Marker lineage   # genomes   # markers   # marker sets   0     1    2    3   4   5+   Completeness   Contamination   Strain heterogeneity
    sql = 'create table IF NOT EXISTS custom_tables.checkm_%s (taxon_id INT, ' \
          ' n_missing INT,' \
          ' n_1 INT,' \
          ' n_2 INT,' \
          ' n_3 INT,' \
          ' n_4 INT,' \
          ' n_5_plus INT,' \
          ' n_total INT,' \
          ' completeness FLOAT,' \
          ' contamination FLOAT,' \
          ' heterogeneity FLOAT);' % biodb
    print sql
    server.adaptor.execute(sql,)
    with open(checkm_file, 'r') as f:
        for n, one_row in enumerate(f):
            data = re.split('\s+', one_row.rstrip())
            print data, len(data)
            if len(data) != 1 and data[1] != 'Bin':
                accession = data[1]
                n_missing = data[6]
                n_1 = data[7]
                n_2 = data[8]
                n_3 = data[9]
                n_4 = data[10]
                n_5_plus = data[11]
                completeness = data[12]
                contamination = data[13]
                heterogeneity = data[14]
                n_sum = int(n_2) + int(n_3)+ int(n_4) + int(n_5_plus)
                print n_sum

                sql = 'insert into  custom_tables.checkm_%s values (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);' % (biodb,
                                                                                  accession2taxon_id[accession],
                                                                                  n_missing,
                                                                                  n_1,
                                                                                  n_2,
                                                                                  n_3,
                                                                                  n_4,
                                                                                  n_5_plus,
                                                                                  n_sum,
                                                                                  completeness,
                                                                                  contamination,
                                                                                  heterogeneity)
                print sql
                server.adaptor.execute(sql,)
        server.commit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--checkm_table', type=str, help="check table")
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)

    args = parser.parse_args()

    import_checkm(args.checkm_table, args.db_name)
