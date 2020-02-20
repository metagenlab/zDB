#! /usr/bin/env python


def biodb2aa_usage(biodb):
    from chlamdb.biosqldb import manipulate_biosqldb
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    from Bio.SeqUtils import CodonUsage
    import biosql_own_sql_tables
    import copy

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

    sql1 = 'select distinct accession from orthology_detail' % biodb
    sql2 = 'select locus_tag, taxon_id from orthology_detail' % biodb
    sql3 = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % biodb

    accession_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql1,)]

    locus2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
    locus2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3,))

    sql_head = 'create table IF NOT EXISTS custom_tables.codon_usage_%s (taxon_id INT, ' \
          ' seqfeature_id INT,' \
          ' seq_length INT,' % biodb
    for aa in CodonsDict.keys():
        sql_head+=' %s FLOAT,' % aa
    sql_head+=' INDEX taxon_id(taxon_id),' \
              ' INDEX seqfeature_id(seqfeature_id));'

    server.adaptor.execute(sql_head,)

    sql_head = 'create table IF NOT EXISTS custom_tables_codon_usage_percent (taxon_id INT, ' \
          ' seqfeature_id INT,' \
          ' seq_length INT,' % biodb
    for aa in CodonsDict.keys():
        sql_head+=' %s FLOAT,' % aa
    sql_head+=' INDEX taxon_id(taxon_id),' \
              ' INDEX seqfeature_id(seqfeature_id));'

    server.adaptor.execute(sql_head,)

    codon_list = [str(i) for i in CodonsDict.keys()]
    count_all=0
    for accession in accession_list:
        record = db.lookup(accession=accession)
        seq = record.seq
        for n, feature in enumerate(record.features):
            if feature.type == 'CDS' and not 'pseudo' in feature.qualifiers:
                count_all+=1
                print count_all
                dna_sequence = feature.extract(seq)
                locus = feature.qualifiers['locus_tag'][0]
                codon_count = copy.copy(CodonsDict)
                for i in range(0, len(dna_sequence), 3):

                    codon = dna_sequence[i:i + 3]
                    #if codon in codon_count:
                    try:
                        codon_count[codon] += 1
                    except KeyError:
                        print 'unknown codon'
                    #else:
                    #    raise TypeError("illegal codon %s in gene: %s" % (codon, feature))

                columns = ','.join(codon_list)
                values = ','.join([str(codon_count[i]) for i in codon_list])

                sql = 'insert into  custom_tables.codon_usage_%s (taxon_id, seqfeature_id, seq_length, %s' \
                      ' ) values (%s, %s, %s, %s);' % (biodb,
                                                       columns,
                                                       locus2taxon_id[locus],
                                                       locus2seqfeature_id[locus],
                                                       len(dna_sequence),
                                                       values)
                print sql
                server.adaptor.execute(sql,)
                n_codons = float(len(dna_sequence)/3)
                values = ','.join([str(codon_count[i]/n_codons) for i in codon_list])
                sql = 'insert into  custom_tables_codon_usage_percent (taxon_id, seqfeature_id, seq_length, %s' \
                      ' ) values (%s, %s, %s, %s);' % (biodb,
                                                       columns,
                                                       locus2taxon_id[locus],
                                                       locus2seqfeature_id[locus],
                                                       len(dna_sequence),
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
