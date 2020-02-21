#! /usr/bin/env python

def import_rnaseq(rnaseq_table, biodb, taxon_id):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)
    import re
    # columns
    # 0 old locus tag (used in this study)
    # 1 gene name
    # 2 new NCBI locus tag
    # 3 temporal class
    # 4 Pfam domains, families, repeats (bit score >=25)
    # 5 EffectiveT3
    # 6 detected in EB proteome
    # 7 function (ConsPred)
    # 8 eggNOG description
    # 9 comment
    # 10 2 hpi_1
    # 11 2 hpi_2
    # 12 2 hpi_3
    # 13 48 hpi_1
    # 14 48 hpi_2
    # 15 48 hpi_3
    # 16 96 hpi_1
    # 17 96 hpi_2
    # 18 96 hpi_3
    # 19 extracellular_1
    # 20 extracellular_2
    # 21 extracellular_3

    sql = 'select old_locus_tag, operon_id from custom_tables_DOOR2_operons' % biodb

    old_locus2operon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % biodb

    locus2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'create table IF NOT EXISTS rnaseq_%s (seqfeature_id INT, ' \
          ' old_locus_tag varchar(400), ' \
          ' operon_id INT,' \
          ' gene_name TEXT,' \
          ' temporal_class varchar(400),' \
          ' EB_proteome INT,' \
          ' function_conspred TEXT,' \
          ' eggNOG TEXT,' \
          ' comment TEXT,' \
          ' hpi_2_1 FLOAT,' \
          ' hpi_2_2 FLOAT,' \
          ' hpi_2_3 FLOAT,' \
          ' hpi_48_1 FLOAT,' \
          ' hpi_48_2 FLOAT,' \
          ' hpi_48_3 FLOAT,' \
          ' hpi_96_1 FLOAT,' \
          ' hpi_96_2 FLOAT,' \
          ' hpi_96_3 FLOAT,' \
          ' extracellular_1 FLOAT,' \
          ' extracellular_2 FLOAT,' \
          ' extracellular_3 FLOAT,' \
          ' INDEX seqfeature_id (seqfeature_id), index temporal_class(temporal_class));' % (biodb,
                                                     taxon_id)
    server.adaptor.execute(sql,)
    server.commit()
    with open(rnaseq_table, 'r') as f:

        for n, row in enumerate(f):
            if n == 0:
                continue
            data = row.rstrip().split('\t')
            try:
                seqfeature_id = locus2seqfeature_id[data[2]]
            except KeyError:
                seqfeature_id = 'NULL'
            operon_id = old_locus2operon_id[data[0]]
            gene_name = data[1]
            old_locus = data[0]
            temporal_class = re.sub(" ","_",data[3])
            if data[6] == 'Y':
                proteome = 1
            else:
                proteome = 0
            conspred = data[7]
            eggNOG = data[8]
            comment = data[9]
            hpi_2_1 = data[10]
            hpi_2_2 = data[11]
            hpi_2_3 = data[12]
            hpi_48_1 = data[13]
            hpi_48_2 = data[14]
            hpi_48_3 = data[15]
            hpi_96_1 = data[16]
            hpi_96_2 = data[17]
            hpi_96_3 = data[18]
            extracellular_1 = data[19]
            extracellular_2 = data[20]
            extracellular_3 = data[21]
            if hpi_2_1 == '-':
                hpi_2_1 = 'NULL'
            if hpi_2_2 == '-':
                hpi_2_2 = 'NULL'
            if hpi_2_3 == '-':
                hpi_2_3 = 'NULL'
            if hpi_48_1 == '-':
                hpi_48_1 = 'NULL'
            if hpi_48_2 == '-':
                hpi_48_2 = 'NULL'
            if hpi_48_3 == '-':
                hpi_48_3 = 'NULL'
            if hpi_96_1 == '-':
                hpi_96_1 = 'NULL'
            if hpi_96_2 == '-':
                hpi_96_2 = 'NULL'
            if hpi_96_3 == '-':
                hpi_96_3 = 'NULL'
            if extracellular_1 == '-':
                extracellular_1 = 'NULL'
            if extracellular_2 == '-':
                extracellular_2 = 'NULL'
            if extracellular_3 == '-':
                extracellular_3 = 'NULL'


            sql = 'insert into rnaseq_%s values (%s, "%s",%s, "%s", "%s", %s, "%s", "%s", "%s", %s, %s' \
                  ' , %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);' % (biodb,
                                                                 taxon_id,
                                                                 seqfeature_id,
                                                                 old_locus,
                                                                 operon_id,
                                                                 gene_name,
                                                                 temporal_class,
                                                                 proteome,
                                                                 conspred,
                                                                 eggNOG,
                                                                 comment,
                                                                 hpi_2_1,
                                                                 hpi_2_2,
                                                                 hpi_2_3,
                                                                 hpi_48_1,
                                                                 hpi_48_2,
                                                                 hpi_48_3,
                                                                 hpi_96_1,
                                                                 hpi_96_2,
                                                                 hpi_96_3,
                                                                 extracellular_1,
                                                                 extracellular_2,
                                                                 extracellular_3)
            try:
                server.adaptor.execute(sql,)
            except:
                print sql
                import sys
                sys.exit()
        server.commit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", '--rnaseq', type=str, help="rnaseq_file")
    parser.add_argument("-t", '--taxon_id', type=str, help="taxon_id")
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)

    args = parser.parse_args()

    import_rnaseq(args.rnaseq, args.db_name, args.taxon_id)


