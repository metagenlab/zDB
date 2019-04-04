#! /usr/bin/env python


def import_16S_identity(identity_table, biodb, sqlite3=False):
    from biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb, sqlite=sqlite3)

    sql = 'select accession, taxon_id from bioentry t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id' \
          ' where t2.name="%s"' % biodb
    accession2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))


    sql = 'create table IF NOT EXISTS comparative_tables.identity_16S_%s (taxon_1 INT, ' \
          ' taxon_2 INT,' \
          ' identity float);' % biodb
    print sql
    server.adaptor.execute(sql,)
    pair_ok = []
    with open(identity_table, 'r') as f:
        for n, one_row in enumerate(f):
            if n == 0:
                accession_list = one_row.rstrip().split('\t')[1:]
                print 'accessions', accession_list
            else:
                data = one_row.rstrip().split('\t')
                accession_1 = data[0]
                for n, identity in enumerate(data[1:]):
                    accession_2 = accession_list[n]
                    print accession_1, accession_2, accession2taxon_id[accession_1],accession2taxon_id[accession_2],identity
                    if [accession2taxon_id[accession_1],accession2taxon_id[accession_2]] not in pair_ok:
                        pair_ok.append([accession2taxon_id[accession_1],accession2taxon_id[accession_2]])

                        sql = 'insert into  comparative_tables.identity_16S_%s values (%s, %s, %s);' % (biodb,
                                                                                                accession2taxon_id[accession_1],
                                                                                                accession2taxon_id[accession_2],
                                                                                                identity)
                        print sql
                        server.adaptor.execute(sql,)
        server.commit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--identity_table', type=str, help="identity table")
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)
    parser.add_argument("-s", '--sqlite', type=str, help="sqlite database path", default=False)

    args = parser.parse_args()

    import_16S_identity(args.identity_table, args.db_name, args.sqlite)
