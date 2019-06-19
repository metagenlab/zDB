#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2019
# ---------------------------------------------------------------------------


def locus2ec_table(hash2ec_dico, biodatabase, hash2locus_list):

    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.biosqldb import biosql_own_sql_tables

    server, db = manipulate_biosqldb.load_db(biodatabase)

    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme.seqfeature_id2ec_%s (enzyme_id INT AUTO_INCREMENT PRIMARY KEY,' \
           ' seqfeature_id INT,' \
           ' ec_id INT,' \
           ' FOREIGN KEY(ec_id) REFERENCES enzymes(enzyme_id));' % (biodatabase)

    server.adaptor.execute_and_fetchall(sql2,)

    sql = 'select locus_tag, seqfeature_id from annotation.seqfeature_id2locus_%s' % biodatabase

    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    for hash in hash2ec_dico:
        for ec_data in hash2ec_dico[hash]:
            for locus in hash2locus_list[hash]:
                sql = 'select enzyme_id from enzyme.enzymes where ec="%s"' % ec_data[0]
                ec_id = server.adaptor.execute_and_fetchall(sql,)[0][0]
                seqfeature_id = locus_tag2seqfeature_id[locus]
                sql = 'insert into enzyme.seqfeature_id2ec_%s (seqfeature_id, ec_id) values (%s, %s)' % (biodatabase,
                                                                                                         seqfeature_id,
                                                                                                         ec_id)
                server.adaptor.execute(sql,)
    server.commit()


def hash2EC(priam_file):

    locus_tag2EC_dico = {}
    with open(priam_file, "r") as f:
        lines = [i.rstrip() for i in f]
        for line in lines:
            if len(line) == 0:
                continue
            elif line[0] == '#':
                continue
            elif line[0] == '>':
                locus_tag = line.split(' ')[0][1:]
                if locus_tag not in locus_tag2EC_dico:
                    locus_tag2EC_dico[locus_tag] = []
            else:
                # valid EC
                locus_tag2EC_dico[locus_tag].append(line.split('\t'))
    return locus_tag2EC_dico


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    import chlamdb_setup_utils

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_priam_files', type=str, help="input interpro csv file", nargs='+', default=False)
    parser.add_argument("-d", '--database_name', type=str, help="database name")
    parser.add_argument("-c", '--corresp_table', type=str, help="hash to locus correspondance table")

    args = parser.parse_args()

    hash2locus_list = chlamdb_setup_utils.get_hash2locus_list(args.corresp_table)

    hash2ec = {}

    for priam_file in args.input_priam_files:
        hash2ec.update(hash2EC(priam_file))

    locus2ec_table(hash2ec,
                   args.database_name,
                   hash2locus_list)
