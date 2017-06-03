#! /usr/bin/env python


def biodb2aa_usage(biodb):
    import manipulate_biosqldb
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    from Bio.SeqUtils import CodonUsage
    import biosql_own_sql_tables

    CodonsDict = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
   'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
   'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
   'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
   'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
   'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
   'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
   'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
   'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
   'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
   'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
   'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
   'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}

    server, db = manipulate_biosqldb.load_db(biodb)

    sql1 = 'select locus_tag, accession from orthology_detail_%s' % biodb
    sql2 = 'select locus_tag, taxon_id from orthology_detail_%s' % biodb
    sql3 = 'select locus_tag, seqfeature_id from custom_tables.locus2seqfeature_id_%s' % biodb
    sql4 = 'select locus_tag, start from orthology_detail_%s' % biodb
    sql5 = 'select locus_tag, stop from orthology_detail_%s' % biodb

    locus2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
    locus2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3,))
    locus2accession = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql1,))
    locus2start = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql4,))
    locus2end = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql5,))

    sql_head = 'create table IF NOT EXISTS custom_tables.codons_usage_%s (taxon_id INT, ' \
          ' seqfeature_id INT,' \
          ' seq_length INT,' % biodb
    for aa in CodonsDict.keys():
        sql_head+=' %s FLOAT,' % aa
    sql_head+=' INDEX taxon_id(taxon_id),' \
              ' INDEX seqfeature_id(seqfeature_id));'

    print sql_head

    server.adaptor.execute(sql_head,)

    for n, locus in enumerate(locus2taxon_id):
        seq_start = int(locus2start[locus])
        seq_end = int(locus2end[locus])
        genome_accession = locus2accession[locus]
        leng = (seq_end-seq_start)
        seq = manipulate_biosqldb.location2sequence(server, genome_accession, biodb, seq_start, leng)

        codon_usage = CodonUsage(seq)
        aa_percent = analysed_seq.get_amino_acids_percent()

        aa_list = [str(i) for i in aa_percent.keys()]
        columns = ','.join(aa_list)
        values = ','.join([str(aa_percent[i]) for i in aa_list])


        sql = 'insert into  custom_tables.aa_usage_count_%s (taxon_id, seqfeature_id, seq_length, %s' \
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
