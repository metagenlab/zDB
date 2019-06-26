#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2019
# ---------------------------------------------------------------------------



def locus2ko_table(hash2ko_dico,
                   biodatabase,
                   ko_accession2ko_id,
                   hash2locus_list):

    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodatabase)

    sql2 = 'select locus_tag, seqfeature_id from annotation.seqfeature_id2locus_%s' % biodatabase

    locus2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2))

    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme.seqfeature_id2ko_%s (seqfeature_id INT,' \
           ' ko_id INT, ' \
           ' thrshld FLOAT, ' \
           ' score FLOAT, ' \
           ' evalue FLOAT, ' \
           ' index ko_id (ko_id),' \
           ' index seqid (seqfeature_id));' % (biodatabase)

    server.adaptor.execute_and_fetchall(sql2,)

    for hash in hash2ko_dico:
        for locus_tag in hash2locus_list[hash]:
            ko = hash2ko_dico[hash]["KO"]
            thrshld = hash2ko_dico[hash]["thrshld"]
            score = hash2ko_dico[hash]["score"]
            evalue = hash2ko_dico[hash]["E-value"]
            ko_id = ko_accession2ko_id[ko]
            seqfeature_id = locus2seqfeature_id[locus_tag]

            sql = 'insert into enzyme.seqfeature_id2ko_%s (seqfeature_id, ko_id, thrshld, score, evalue) values (%s, %s, %s, %s, %s)' % (biodatabase,
                                                                                                                                         seqfeature_id,
                                                                                                                                         ko_id,
                                                                                                                                         thrshld,
                                                                                                                                         score,
                                                                                                                                         evalue)


            server.adaptor.execute(sql,)
    server.commit()


def locus2ko_table_legacy(biodatabase, hash2ko, hash2locus_list):
    # create legacy table locus2ko
    # TODO: remove all depenancies to this table
    # content:

    '''
    taxon_id, locus_tag, orthogroup, ko_id (not numerical id but KO*)
    '''

    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodatabase)

    sql = 'select locus_tag, taxon_id from orthology_detail_%s' % biodatabase
    sql2 = 'select locus_tag, orthogroup from orthology_detail_%s' % biodatabase

    locus2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql))
    locus2orthogroup = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2))

    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme.locus2ko_%s (taxon_id INT,'\
           ' locus_tag VARCHAR(20),' \
           ' orthogroup varchar(20),' \
           ' ko_id VARCHAR(20), index taxon_id (taxon_id), index ko_id (ko_id));' % (biodatabase)

    server.adaptor.execute_and_fetchall(sql2,)
    for hash in hash2ko:
        for locus in hash2locus_list[hash]:
            ko = hash2ko[hash]["KO"]

            sql = 'insert into enzyme.locus2ko_%s (taxon_id, locus_tag, orthogroup, ko_id) values ("%s", "%s", "%s", "%s")' % (biodatabase,
                                                                                                                               locus2taxon_id[locus],
                                                                                                                               locus,
                                                                                                                               locus2orthogroup[locus],
                                                                                                                               ko)

            server.adaptor.execute(sql,)
        server.commit()


def parse_kofamscan_output(result_file_list):
    hash2ko = {}
    for result_file in result_file_list:
        # # gene name            KO     thrshld  score   E-value KO definition
        with open(result_file, 'r') as f:
            for line in f:
                if line.startswith("#"):
                    continue
                if line.startswith("*"):
                    data = line[2:].strip().split()
                    hash2ko[data[0]] = {}
                    hash2ko[data[0]]["KO"] = data[1]
                    hash2ko[data[0]]["thrshld"] = data[2]
                    hash2ko[data[0]]["score"] = data[3]
                    hash2ko[data[0]]["E-value"] = data[4]
    return hash2ko


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    import chlamdb_setup_utils

    parser = argparse.ArgumentParser()
    parser.add_argument("-k", '--ko_table_list', type=str, help="input blastGhost file", nargs='+')
    parser.add_argument("-d", '--database_name', type=str, help="database name")
    parser.add_argument("-c", '--corresp_table', type=str, help="hash to locus correspondance table")
    parser.add_argument("-l", '--legacy', action='store_true', help="Create legacy table(s)")

    args = parser.parse_args()

    hash2locus_list = chlamdb_setup_utils.get_hash2locus_list(args.corresp_table)

    server, db = manipulate_biosqldb.load_db(args.database_name)
    sql = 'select ko_accession, ko_id from enzyme.ko_annotation'
    ko_accession2ko_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    hash2ko = parse_kofamscan_output(args.ko_table_list)

    locus2ko_table(hash2ko,
                   args.database_name,
                   ko_accession2ko_id,
                   hash2locus_list)

    if args.legacy:
        locus2ko_table_legacy(args.database_name, hash2ko, hash2locus_list)
