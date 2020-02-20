#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2019
# ---------------------------------------------------------------------------


def load_BPBAac_table(table_file, 
                      biodb, 
                      hash2locus_list):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)
    
    sql = 'create database if not exists effectors;'
    server.adaptor.execute(sql,)

    sql = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id'
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag, taxon_id from orthology_detail'
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'CREATE table effectors_predicted_BPBAac (seqfeature_id INT primary key,taxon_id INT,SVM_value FLOAT,effectors INT,' \
          ' INDEX taxon_id (taxon_id))'

    server.adaptor.execute(sql,)

    with open(table_file, 'r') as f:
        for row in f:
            data = row.rstrip().split(',')
            if data[0] == 'Protein':
                continue
            else:
                hash = data[0]
                if data[2] == 'NO':
                    effector = 0
                    continue
                elif data[2] == 'YES':
                    effector = 1
                else:
                    print ('unknown result', data)
                for locus in hash2locus_list[hash]:
                    sql = 'insert into effectors_predicted_BPBAac values (%s,%s,%s,%s)' % (locus_tag2seqfeature_id[locus],
                                                                                           locus_tag2taxon_id[locus],
                                                                                           data[1],
                                                                                           effector)
                    server.adaptor.execute(sql,)
                server.commit()


def load_T3_MM_table(table_file, 
                     biodb, 
                     hash2locus_list):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create database if not exists effectors;'
    server.adaptor.execute(sql,)

    sql = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id'
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag, taxon_id from orthology_detail'
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'CREATE table effectors_predicted_T3MM (seqfeature_id INT primary key,taxon_id INT,value FLOAT,effectors INT,' \
            ' probability FLOAT, INDEX taxon_id (taxon_id))'

    server.adaptor.execute(sql,)

    with open(table_file, 'r') as f:
        for row in f:
            data = row.rstrip().split(',')
            if data[0] == 'seqName':
                continue
            else:
                hash = data[0]
                if data[2] == 'NO':
                    effector = 0
                    continue
                elif data[2] == 'YES':
                    effector = 1
                else:
                    print ('unknown result', data)
                for locus in hash2locus_list[hash]:
                    
                    # seqName,value,T3SE,probability
                    sql = 'insert into effectors_predicted_T3MM values (%s,%s,%s,%s, %s)' % (locus_tag2seqfeature_id[locus],
                                                                                             locus_tag2taxon_id[locus],
                                                                                             data[1],
                                                                                             effector,
                                                                                             data[3])
                    server.adaptor.execute(sql,)
                server.commit()


def load_effectiveT3_table(table_file, 
                           biodb,
                           hash2locus_list):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create database if not exists effectors;'
    server.adaptor.execute(sql,)

    sql = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id'
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag, taxon_id from orthology_detail'
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'CREATE table IF NOT EXISTS effectors_predicted_effectiveT3 (seqfeature_id INT primary key,taxon_id INT,score FLOAT,effectors INT,' \
            ' INDEX taxon_id (taxon_id))'

    server.adaptor.execute(sql,)

    with open(table_file, 'r') as f:
        for i, row in enumerate(f):
            if row[0] == '#':
                continue
            data = row.rstrip().split(';')
            if data[0] == 'Id':
                continue
            hash = data[0]
            if data[3] == 'false':
                effector = 0
                continue
            elif data[3] == 'true':
                effector = 1
            else:
                print ('unknown result', data)

            for locus in hash2locus_list[hash]:
                # CRC-02125D79C6DFBC1D; RhT_01099 Rhabdochlamydia helvetica T3358;1.000000000000000;true
                sql = 'insert into effectors_predicted_effectiveT3 values (%s,%s,%s,%s)' % (locus_tag2seqfeature_id[locus],
                                                                                            locus_tag2taxon_id[locus],
                                                                                            data[2],
                                                                                            effector)
                server.adaptor.execute(sql,)
            server.commit()
                

def load_DeepT3_table(table_file, 
                      biodb,
                      hash2locus_list):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create database if not exists effectors;'
    server.adaptor.execute(sql,)

    sql = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id'
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag, taxon_id from orthology_detail'
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'CREATE table IF NOT EXISTS effectors_predicted_DeepT3 (seqfeature_id INT primary key,taxon_id INT,effectors INT,' \
            ' INDEX taxon_id (taxon_id))'

    server.adaptor.execute(sql,)

    with open(table_file, 'r') as f:
        for i, row in enumerate(f):
            data = row.rstrip().split('\t')
            if i < 3:
                continue
            else:
                hash = data[0]
                if data[1] == 'non-T3SE':
                    effector = 0
                    continue
                elif data[1] == 'T3SE':
                    effector = 1
                else:
                    print ('unknown result', data)
                for locus in hash2locus_list[hash]:
                    # CRC-20ACA5744563D026	non-T3SE
                    sql = 'insert into effectors_predicted_DeepT3 values (%s,%s,%s)' % (locus_tag2seqfeature_id[locus],
                                                                                        locus_tag2taxon_id[locus],
                                                                                        effector)
                    server.adaptor.execute(sql,)
                server.commit()
                                

if __name__ == '__main__':
    import argparse
    import chlamdb_setup_utils
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--table_file', type=str, help="effector prediction table", required=True)
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)
    parser.add_argument("-b", '--BPBAac', action="store_true", help="BPBAac table")
    parser.add_argument("-t3", '--T3_MM', action="store_true", help="T3_MM table")
    parser.add_argument("-eff", '--effectiveT3', action="store_true", help="effectiveT3 table")
    parser.add_argument("-dt", '--DeepT3', action="store_true", help="DeepT3 table")
    parser.add_argument("-c", '--corresp_table',  type=str, help="hash to locus correspondance table")


    args = parser.parse_args()
    
    hash2locus_list = chlamdb_setup_utils.get_hash2locus_list(args.corresp_table)

    if args.BPBAac:
        load_BPBAac_table(args.table_file, 
                          args.db_name,
                          hash2locus_list)
    elif args.T3_MM:
        load_T3_MM_table(args.table_file, 
                         args.db_name,
                         hash2locus_list)
    elif args.effectiveT3:
        load_effectiveT3_table(args.table_file, 
                               args.db_name,
                               hash2locus_list)
    elif args.DeepT3:
        load_DeepT3_table(args.table_file, 
                          args.db_name,
                          hash2locus_list)
