#! /usr/bin/env python



def load_BPBAac_table(table_file, biodb='chlamydia_04_16'):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % biodb
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag, taxon_id from orthology_detail' % biodb
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'CREATE table effectors_predicted_BPBAac (seqfeature_id INT primary key,taxon_id INT,SVM_value FLOAT,effectors INT,' \
          ' INDEX taxon_id (taxon_id))' % biodb

    server.adaptor.execute(sql,)

    with open(table_file, 'r') as f:
        for row in f:
            data = row.rstrip().split(',')
            if data[0] == 'Protein':
                continue
            else:

                if data[2] == 'NO':
                    effector = 0
                    continue
                elif data[2] == 'YES':
                    effector = 1
                else:
                    print 'unknown result', data
                sql = 'insert into effectors_predicted_BPBAac values (%s,%s,%s,%s)' % (biodb,
                                                                                          locus_tag2seqfeature_id[data[0]],
                                                                                          locus_tag2taxon_id[data[0]],
                                                                                          data[1],
                                                                                          effector)
                server.adaptor.execute(sql,)
                server.commit()


def load_T3_MM_table(table_file, biodb='chlamydia_04_16'):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % biodb
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag, taxon_id from orthology_detail' % biodb
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'CREATE table effectors_predicted_T3MM (seqfeature_id INT primary key,taxon_id INT,value FLOAT,effectors INT,' \
          ' probability FLOAT, INDEX taxon_id (taxon_id))' % biodb

    server.adaptor.execute(sql,)

    with open(table_file, 'r') as f:
        for row in f:
            data = row.rstrip().split(',')
            if data[0] == 'seqName':
                continue
            else:

                if data[2] == 'NO':
                    effector = 0
                    continue
                elif data[2] == 'YES':
                    effector = 1
                else:
                    print 'unknown result', data

                # seqName,value,T3SE,probability
                sql = 'insert into effectors_predicted_T3MM values (%s,%s,%s,%s, %s)' % (biodb,
                                                                                            locus_tag2seqfeature_id[data[0]],
                                                                                            locus_tag2taxon_id[data[0]],
                                                                                            data[1],
                                                                                            effector,
                                                                                            data[3])
                server.adaptor.execute(sql,)
                server.commit()


def load_effectiveT3_table(table_file, biodb='chlamydia_04_16'):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % biodb
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag, taxon_id from orthology_detail' % biodb
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'CREATE table IF NOT EXISTS effectors_predicted_effectiveT3 (seqfeature_id INT primary key,taxon_id INT,score FLOAT,effectors INT,' \
          ' INDEX taxon_id (taxon_id))' % biodb

    server.adaptor.execute(sql,)

    with open(table_file, 'r') as f:
        for i, row in enumerate(f):
            data = row.rstrip().split(';')
            if i < 3:
                continue
            else:

                if data[3] == 'false':
                    effector = 0
                    continue
                elif data[3] == 'true':
                    effector = 1
                else:
                    print 'unknown result', data

                # Id; Description; Score; is Effective
                sql = 'insert into effectors_predicted_effectiveT3 values (%s,%s,%s,%s)' % (biodb,
                                                                                            locus_tag2seqfeature_id[data[0]],
                                                                                            locus_tag2taxon_id[data[0]],
                                                                                            data[2],
                                                                                            effector)
                server.adaptor.execute(sql,)
                server.commit()


def load_T4SEpre_bpbAac_table(table_file, biodb='chlamydia_04_16'):
    from chlamdb.biosqldb import manipulate_biosqldb
    import re

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % biodb
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag, taxon_id from orthology_detail' % biodb
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'CREATE table IF NOT EXISTS effectors_predicted_T4SEpre_bpbAac (seqfeature_id INT primary key,taxon_id INT,SVM_value FLOAT,effectors INT,' \
          ' INDEX taxon_id (taxon_id))' % biodb

    server.adaptor.execute(sql,)

    with open(table_file, 'r') as f:
        for i, row in enumerate(f):
            data = row.rstrip().split(',')
            if i < 2:
                continue
            else:

                if data[2] == 'NO':
                    effector = 0
                    continue
                elif data[2] == 'YES':
                    effector = 1
                else:
                    print 'unknown result', data

                #Protein,SVM-Value,T4S protein or not
                locus = re.sub('"', "", data[0])
                sql = 'insert into effectors_predicted_T4SEpre_bpbAac values (%s,%s,%s,%s)' % (biodb,
                                                                                            locus_tag2seqfeature_id[locus],
                                                                                            locus_tag2taxon_id[locus],
                                                                                            data[1],
                                                                                            effector)
                server.adaptor.execute(sql,)
                server.commit()


def load_T4SEpre_psAac_table(table_file, biodb='chlamydia_04_16'):
    from chlamdb.biosqldb import manipulate_biosqldb
    import re

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % biodb
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag, taxon_id from orthology_detail' % biodb
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'CREATE table IF NOT EXISTS effectors_predicted_T4SEpre_psAac (seqfeature_id INT primary key,taxon_id INT,SVM_value FLOAT,effectors INT,' \
          ' INDEX taxon_id (taxon_id))' % biodb

    server.adaptor.execute(sql,)

    with open(table_file, 'r') as f:
        for i, row in enumerate(f):
            data = row.rstrip().split(',')
            if i < 2:
                continue
            else:

                if data[2] == 'NO':
                    effector = 0
                    continue
                elif data[2] == 'YES':
                    effector = 1
                else:
                    print 'unknown result', data

                #Protein,SVM-Value,T4S protein or not
                #locus = re.sub('"', "", data[0])
                sql = 'insert into effectors_predicted_T4SEpre_psAac values (%s,%s,%s,%s)' % (biodb,
                                                                                            locus_tag2seqfeature_id[data[0]],
                                                                                            locus_tag2taxon_id[data[0]],
                                                                                            data[1],
                                                                                            effector)
                server.adaptor.execute(sql,)
                server.commit()


def load_chaperones_table(table_file, biodb='chlamydia_04_16'):
    from chlamdb.biosqldb import manipulate_biosqldb
    import re

    '''
    ATTENTION possible d'avoir plusieurs sequences cibles identifiees dans une meme sequence: dubplicate seqfeature dans la table
    '''

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % biodb
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag, taxon_id from orthology_detail' % biodb
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'CREATE table IF NOT EXISTS effectors_predicted_chaperones (seqfeature_id INT,taxon_id INT,motif varchar(400), info TEXT,' \
          ' INDEX taxon_id (taxon_id), index seqfeature_id (seqfeature_id))' % biodb

    server.adaptor.execute(sql,)

    with open(table_file, 'r') as f:
        for i, row in enumerate(f):
            data = row.rstrip().split('\t')
            if data[1] == 'Protein does not contain a CCBD within the first 150 amino acids':
                continue
            else:
                motif = data[1].split('(')[1].split(')')[0]
                info = data[1].split(') ')[1]

            #Protein,SVM-Value,T4S protein or not
            #locus = re.sub('"', "", data[0])
            sql = 'insert into effectors_predicted_chaperones values (%s,%s,"%s", "%s")' % (biodb,
                                                                                        locus_tag2seqfeature_id[data[0]],
                                                                                        locus_tag2taxon_id[data[0]],
                                                                                        motif,
                                                                                        info)
            try:
                server.adaptor.execute(sql,)
                server.commit()
            except:
                print data
                import sys
                sys.exit()


def load_ELD_table(table_file, biodb='chlamydia_04_16'):
    from chlamdb.biosqldb import manipulate_biosqldb
    import re

    '''
    ATTENTION possible d'avoir plusieurs domaines identifiees dans une meme sequence: dubplicate seqfeature dans la table
    '''

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % biodb
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag, taxon_id from orthology_detail' % biodb
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'CREATE table IF NOT EXISTS effectors_predicted_ELD (seqfeature_id INT,taxon_id INT,pfam_domain varchar(400), description TEXT, score INT,' \
          ' INDEX taxon_id (taxon_id), index seqfeature_id (seqfeature_id))' % biodb

    server.adaptor.execute(sql,)

    with open(table_file, 'r') as f:
        for i, row in enumerate(f):
            data = row.rstrip().split('\t')

            #Protein,SVM-Value,T4S protein or not
            #locus = re.sub('"', "", data[0])
            sql = 'insert into effectors_predicted_ELD values (%s,%s,"%s", "%s", %s)' % (biodb,
                                                                                        locus_tag2seqfeature_id[data[0]],
                                                                                        locus_tag2taxon_id[data[0]],
                                                                                        data[1],
                                                                                        data[2],
                                                                                        data[3])
            try:
                server.adaptor.execute(sql,)
                server.commit()
            except:
                print data
                import sys
                sys.exit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--table_file', type=str, help="effector prediction table", required=True)
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)
    parser.add_argument("-b", '--BPBAac', action="store_true", help="BPBAac table")
    parser.add_argument("-t3", '--T3_MM', action="store_true", help="T3_MM table")
    parser.add_argument("-eff", '--effectiveT3', action="store_true", help="effectiveT3 table")
    parser.add_argument("-T4b", '--T4SEpre_bpbAac', action="store_true", help="T4SEpre_bpbAac table")
    parser.add_argument("-T4p", '--T4SEpre_psAac', action="store_true", help="T4SEpre_psAac table ")
    parser.add_argument("-c", '--chaperones', action="store_true", help="CHAPERONES table ")
    parser.add_argument("-e", '--ELD', action="store_true", help="ELD table ")


    args = parser.parse_args()

    if args.BPBAac:
        load_BPBAac_table(args.table_file, args.db_name)
    elif args.T3_MM:
        load_T3_MM_table(args.table_file, args.db_name)
    elif args.effectiveT3:
        load_effectiveT3_table(args.table_file, args.db_name)
    elif args.T4SEpre_bpbAac:
        load_T4SEpre_bpbAac_table(args.table_file, args.db_name)
    elif args.T4SEpre_psAac:
        load_T4SEpre_psAac_table(args.table_file, args.db_name)
    elif args.chaperones:
        load_chaperones_table(args.table_file, args.db_name)
    elif args.ELD:
        load_ELD_table(args.table_file, args.db_name)