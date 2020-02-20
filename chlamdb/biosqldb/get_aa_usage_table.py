#! /usr/bin/env python


def biodb2aa_usage(biodb):
    from chlamdb.biosqldb import manipulate_biosqldb
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    '''
    import annotation file of the form as tabulated file with 4 columns:
    1. category
    2. gene name
    3. locus tag
    4. description
    '''

    server, db = manipulate_biosqldb.load_db(biodb)

    sql1 = 'select locus_tag, translation from orthology_detail' % biodb
    sql2 = 'select locus_tag, taxon_id from orthology_detail' % biodb
    sql3 = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % biodb

    locus2translation = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql1,))
    locus2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
    locus2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3,))

    all_aa = []
    for locus in locus2translation:
        for aa in list(locus2translation[locus]):
            if aa not in all_aa:
                all_aa.append(aa)
    print len(all_aa), all_aa


    sql_head = 'create table IF NOT EXISTS custom_tables_aa_usage_count (taxon_id INT, ' \
          ' seqfeature_id INT,' \
          ' seq_length INT,' % biodb
    for aa in all_aa:
        sql_head+=' %s FLOAT,' % aa
    sql_head+=' INDEX taxon_id(taxon_id),' \
              ' INDEX seqfeature_id(seqfeature_id));'

    print sql_head

    server.adaptor.execute(sql_head,)

    for n, locus in enumerate(locus2translation):
        seq = locus2translation[locus]
        analysed_seq = ProteinAnalysis(seq)
        aa_percent = analysed_seq.get_amino_acids_percent()

        aa_list = [str(i) for i in aa_percent.keys()]
        columns = ','.join(aa_list)
        values = ','.join([str(aa_percent[i]) for i in aa_list])


        sql = 'insert into  custom_tables_aa_usage_count (taxon_id, seqfeature_id, seq_length, %s' \
              ' ) values (%s, %s, %s, %s);' % (biodb,
                                               columns,
                                               locus2taxon_id[locus],
                                               locus2seqfeature_id[locus],
                                               len(seq),
                                               values)
        print sql
        server.adaptor.execute(sql,)
        server.commit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)

    args = parser.parse_args()

    biodb2aa_usage(args.db_name)
