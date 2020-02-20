#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

def get_pepstats(locus2seq):

    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    import re

    locus2pepstats = {}
    for locus in locus2seq:
            seq = re.sub('X|J|Z|B','',locus2seq[locus])
            analysed_seq = ProteinAnalysis(str(seq))
            locus2pepstats[locus] = analysed_seq
    return locus2pepstats


def all_aa_possibles(locus2sequence):

    unique_list = []
    for locus in locus2sequence:
        seq = locus2sequence[locus]
        for aa in set(list(str(seq))):
            if aa not in unique_list:
                unique_list.append(aa)
    return unique_list


def pepstats2sql(biodb, locus2pepstats, aa_list):

    from chlamdb.biosqldb import manipulate_biosqldb
    import re
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select locus_tag, taxon_id from orthology_detail'

    locus2taxon = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql1 = 'CREATE table custom_tables_locus2pepstats (locus_tag varchar(300),' \
           ' taxon_id INT,' \
           ' mol_weight float,' \
           ' isoelectric_point float,' \
           ' aromaticity float,' \
           ' instability_index float,' \
           ' fraction_helix float,' \
           ' fraction_turn float,' \
           ' fraction_sheet float,' \
           ' index taxon_id(taxon_id),' \
           ' index locus_tag(locus_tag))'



    columns = 'locus_tag varchar(300),taxon_id INT,' + ' float,'.join(aa_list) + ' float,'
    sql2 = 'CREATE table custom_tables_locus2aa_freq (%s index taxon_id(taxon_id),index locus_tag(locus_tag))' % (columns)
    sql3 = 'CREATE table custom_tables_locus2aa_counts (%s index taxon_id(taxon_id),index locus_tag(locus_tag))' % (columns)



    server.adaptor.execute(sql1,)
    server.adaptor.execute(sql2,)
    server.adaptor.execute(sql3,)

    for locus in locus2pepstats:
        mol_weight = locus2pepstats[locus].molecular_weight()
        isoelectric_point = locus2pepstats[locus].isoelectric_point()
        aa_counts = locus2pepstats[locus].count_amino_acids()
        aa_percent = locus2pepstats[locus].get_amino_acids_percent()
        aromaticity = locus2pepstats[locus].aromaticity()
        instability_index = locus2pepstats[locus].instability_index()
        secondary_structure_fraction = locus2pepstats[locus].secondary_structure_fraction() # The list contains 3 values: [Helix, Turn, Sheet].

        sql = 'insert into custom_tables_locus2pepstats values ("%s",%s,%s,%s,%s,%s,%s,%s,%s)' % (locus,
                                                                                                  locus2taxon[locus],
                                                                                                  mol_weight,
                                                                                                  isoelectric_point,
                                                                                                  aromaticity,
                                                                                                  instability_index,
                                                                                                  secondary_structure_fraction[0],
                                                                                                  secondary_structure_fraction[1],
                                                                                                  secondary_structure_fraction[2])

        sql_freq = 'insert into custom_tables_locus2aa_freq values ("%s",%s,' % (locus, 
                                                                                 locus2taxon[locus])
        sql_counts = 'insert into custom_tables_locus2aa_counts values ("%s",%s,' % (locus, 
                                                                                     locus2taxon[locus])
        for aa in aa_list:
            try:
                sql_freq+= '%s,' % aa_percent[aa]
                sql_counts+= '%s,' % aa_counts[aa]
            except:
                sql_freq+= '0,'
                sql_counts+= '0,'
        sql_freq = sql_freq[0:-1] + ')'
        sql_counts = sql_counts[0:-1] + ')'

        server.adaptor.execute(sql,)
        server.adaptor.execute(sql_freq,)
        server.adaptor.execute(sql_counts,)
        server.commit()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_fasta', type=str, help="input fasta file")
    parser.add_argument("-d", '--biodb', type=str, help="biodatabase name")


    args = parser.parse_args()
    if args.input_fasta:
        with open(args.input_fasta) as f:
            from Bio import SeqIO
            record = list(SeqIO.parse(f, 'fasta'))
            locus2seq = {}
            for i in record:
                locus2seq[i.name] = i.seq
            aa_unique_list = all_aa_possibles(locus2seq)
            stats = get_pepstats(locus2seq)
            #pepstats2sql(args.biodb,
            #             stats,
            #             aa_unique_list)
    if args.biodb:
        from chlamdb.biosqldb import manipulate_biosqldb
        server, db = manipulate_biosqldb.load_db(args.biodb)

        sql = 'select locus_tag, translation from orthology_detail'

        locus2seq = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

        aa_unique_list = all_aa_possibles(locus2seq)
        stats = get_pepstats(locus2seq)
        pepstats2sql(args.biodb,
                     stats,
                     aa_unique_list)
