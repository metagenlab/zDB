#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

def parse_seqid_file(orthofinder_seqid_file):
    locus_tag2orthofinder_id = {}
    with open(orthofinder_seqid_file, 'r') as f:
        for line in f:
            data = line.split(': ')
            seq_id = data[0]
            locus_tag = data[1].split(' ')[0]
            locus_tag2orthofinder_id[seq_id] = locus_tag
    return locus_tag2orthofinder_id

def parse_orthofinder_blast_files(biodb,
                                  orthofinder_blast_file_list,
                                  locus_tag2orthofinder_id,
                                  evalue_cutoff=0.01,
                                  identity_cutoff=0):

    import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select locus_tag, taxon_id from orthology_detail_%s' % biodb

    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))



    locus2taxon2best_hit_id = {}
    for file in orthofinder_blast_file_list:
        with open(file, 'r') as f:
            for line in f:
                #9_1	60_294	56.647	346	139	3	4	338	5	350	3.33e-141	405
                data = line.split('\t')
                seq_1_id = data[0]
                seq_2_id = data[1]
                locus_tag_1 = locus_tag2orthofinder_id[seq_1_id]
                locus_tag_2 = locus_tag2orthofinder_id[seq_2_id]
                identity = data[2]
                evalue = data[10]
                bitscore = data[11]
                taxon_id_seq_2 = locus_tag2taxon_id[locus_tag_2]
                if locus_tag_1 not in locus2taxon2best_hit_id:
                    locus2taxon2best_hit_id[locus_tag_1] = {}
                    locus2taxon2best_hit_id[locus_tag_1][taxon_id_seq_2] = [locus_tag_2, identity, evalue, bitscore]
                else:
                    # not the best hit, skip
                    if taxon_id_seq_2 in locus2taxon2best_hit_id[locus_tag_1]:
                        continue
                    # best hit
                    else:
                        locus2taxon2best_hit_id[locus_tag_1][taxon_id_seq_2] = [locus_tag_2, identity, evalue, bitscore]
    return locus2taxon2best_hit_id

def get_reciproval_BBH_table(biodb, locus2taxon2best_hit_id):

    # taxon_id, taxon_2, seqfeature_id_1, seqfeature_id_2, blast_identity_a_vs_b, blast_identity_b_vs_a, mean_blast_identity, msa_identity, orthogroup

    import manipulate_biosqldb
    import numpy
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table comparative_tables.reciprocal_BBH_%s (taxon_1 INT, taxon_2 INT, seqfeature_id_1 INT, seqfeature_id_2 INT,' \
          ' blast_evalue_a_vs_b FLOAT, blast_evalue_b_vs_a FLOAT,' \
          ' blast_score_a_vs_b FLOAT, blast_score_b_vs_a FLOAT, blast_identity_a_vs_b FLOAT, blast_identity_b_vs_a FLOAT, ' \
          ' mean_blast_identity FLOAT, msa_identity FLOAT, orthogroup_1 varchar(400), orthogroup_2 varchar(400),' \
          ' index taxon_1(taxon_1), INDEX taxon_2 (taxon_2), INDEX seqfeature_id_1(seqfeature_id_1),' \
          ' INDEX seqfeature_id_2 (seqfeature_id_2))' % biodb
    server.adaptor.execute(sql,)
    server.commit()

    sql = 'select locus_tag, taxon_id from orthology_detail_%s' % biodb

    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag, orthogroup from orthology_detail_%s' % biodb

    locus_tag2orthogroup = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag, seqfeature_id from custom_tables.locus2seqfeature_id_%s' % biodb

    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    for n1, locus_1 in enumerate(locus2taxon2best_hit_id):
        taxon_1 = locus_tag2taxon_id[locus_1]
        for taxon_2 in locus2taxon2best_hit_id[locus_1]:
            locus_2 = locus2taxon2best_hit_id[locus_1][taxon_2][0]
            try:
                if locus2taxon2best_hit_id[locus_2][taxon_1][0] == locus_1:

                        identity_a_vs_b = locus2taxon2best_hit_id[locus_1][taxon_2][1]
                        identity_b_vs_a = locus2taxon2best_hit_id[locus_2][taxon_1][1]

                        evalue_a_vs_b = locus2taxon2best_hit_id[locus_1][taxon_2][2]
                        evalue_b_vs_a = locus2taxon2best_hit_id[locus_2][taxon_1][2]

                        score_a_vs_b = locus2taxon2best_hit_id[locus_1][taxon_2][3]
                        score_b_vs_a = locus2taxon2best_hit_id[locus_2][taxon_1][3]

                        group_1 = locus_tag2orthogroup[locus_1]
                        group_2 = locus_tag2orthogroup[locus_2]
                        seqfeature_id_1 = locus_tag2seqfeature_id[locus_1]
                        seqfeature_id_2 = locus_tag2seqfeature_id[locus_2]
                        if group_1 != group_2:
                            msa_identity = "NULL"
                        else:
                            sql = 'select identity from orth_%s.%s where locus_a in ("%s", "%s") ' \
                                  ' and locus_b in ("%s","%s") and locus_a!=locus_b;' % (biodb,
                                                                                         group_1,
                                                                                         locus_1,
                                                                                         locus_2,
                                                                                         locus_1,
                                                                                         locus_2)

                            msa_identity = server.adaptor.execute_and_fetchall(sql,)[0][0]

                        sql = 'insert into comparative_tables.reciprocal_BBH_%s values (%s,%s,%s,%s,%s,%s,%s,%s' \
                              ',%s,%s,%s,%s,"%s","%s")' % (biodb,
                                                      taxon_1,
                                                      taxon_2,
                                                      seqfeature_id_1,
                                                      seqfeature_id_2,
                                                      evalue_a_vs_b,
                                                      evalue_b_vs_a,
                                                      score_a_vs_b,
                                                      score_b_vs_a,
                                                      identity_a_vs_b,
                                                      identity_b_vs_a,
                                                      numpy.mean( [float(identity_a_vs_b), float(identity_b_vs_a)]),
                                                      msa_identity,
                                                      group_1,
                                                      group_2)
                        try:
                            server.adaptor.execute_and_fetchall(sql,)
                        except:
                            print sql
                        server.commit()


            except:
                print 'FAIL!'
                continue


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s",'--seq_id_file',type=str,help="SequenceIDs.txt")
    parser.add_argument("-b",'--blast_files', nargs="+", help="blast_file")
    parser.add_argument("-d",'--database', help="bioadatabase name")

    args = parser.parse_args()
    seqid2locus_tag = parse_seqid_file(args.seq_id_file)
    locus2taxon2best_hit_id = parse_orthofinder_blast_files(args.database,
                                                            args.blast_files,
                                                            seqid2locus_tag)
    get_reciproval_BBH_table(args.database,
                             locus2taxon2best_hit_id)