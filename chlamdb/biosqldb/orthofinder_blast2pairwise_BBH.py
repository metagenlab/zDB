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
                                  identity_cutoff=0,
                                  sqlite3=False):

    from chlamdb.biosqldb import manipulate_biosqldb
    if not sqlite3:
        server, db = manipulate_biosqldb.load_db(biodb)
    else:
        server, db = manipulate_biosqldb.load_db(biodb, sqlite3)

    sql = 'select locus_tag, taxon_id from annotation_seqfeature_id2locus' % biodb

    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))



    locus2taxon2best_hit_id = {}
    for n, file in enumerate(orthofinder_blast_file_list):
        print(n, len(orthofinder_blast_file_list))
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
                query_align_length = int(data[7]) - int(data[6])
                hit_align_length = int(data[9]) - int(data[8])
                taxon_id_seq_2 = locus_tag2taxon_id[locus_tag_2]
                if locus_tag_1 not in locus2taxon2best_hit_id:
                    locus2taxon2best_hit_id[locus_tag_1] = {}
                    locus2taxon2best_hit_id[locus_tag_1][taxon_id_seq_2] = [locus_tag_2,
                                                                            identity,
                                                                            evalue,
                                                                            bitscore,
                                                                            query_align_length,
                                                                            hit_align_length]
                else:
                    # not the best hit, skip
                    if taxon_id_seq_2 in locus2taxon2best_hit_id[locus_tag_1]:
                        continue
                    # best hit
                    else:
                        locus2taxon2best_hit_id[locus_tag_1][taxon_id_seq_2] = [locus_tag_2,
                                                                                identity,
                                                                                evalue,
                                                                                bitscore,
                                                                                query_align_length,
                                                                                hit_align_length]
    return locus2taxon2best_hit_id

def get_parwise_genome_median_identity_table(biodb, sqlite3=False):

    from chlamdb.biosqldb import manipulate_biosqldb
    import numpy
    import rpy2.robjects.numpy2ri

    import numpy as np
    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri as numpy2ri
    from rpy2.robjects.packages import importr
    import rpy2.rlike.container as rlc
    from rpy2.robjects import pandas2ri
    pandas2ri.activate()

    server, db = manipulate_biosqldb.load_db(biodb, sqlite3)

    sql = 'select taxon_id from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id ' \
          ' where t1.name="%s" group by taxon_id' % biodb
    taxon_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    sql = 'CREATE table if not exists comparative_tables_reciprocal_BBH_average_identity (taxon_1 INT, taxon_2 INT, ' \
          ' average_identity FLOAT, median_identity FLOAT, n_pairs INT)' % biodb
    server.adaptor.execute(sql,)
    server.commit()
    for n, taxon_1 in enumerate(taxon_list):
        print ("%s / %s" % (n, len(taxon_list)))
        for taxon_2 in taxon_list[n+1:len(taxon_list)]:
            sql = 'select blast_identity_a_vs_b from comparative_tables_reciprocal_BBH where taxon_1 =%s ' \
                  'and taxon_2=%s;' % (biodb,
                                taxon_1,
                                taxon_2)

            data = [float(i[0]) for i in server.adaptor.execute_and_fetchall(sql,)]
            median_id = numpy.median(data)
            mean_id = numpy.mean(data)

            data_array = np.asarray(data, dtype='float')
            robjects.r.assign('identity_values', numpy2ri.numpy2ri(data_array))
            robjects.r('''
                d <- density(identity_values)
                w <- which(d$y==max(d$y))
                max_density <- d$x[w]
            ''')
            max_density = robjects.r('max_density')

            sql = 'insert into comparative_tables_reciprocal_BBH_average_identity values(%s, %s, %s, %s, %s, %s)' % (biodb,
                                                                                                                        taxon_1,
                                                                                                                        taxon_2,
                                                                                                                        round(mean_id, 2),
                                                                                                                        round(median_id, 2),
                                                                                                                        round(float(max_density[0]), 2),
                                                                                                                        len(data)
                                                                                                                        )
            server.adaptor.execute(sql,)
        server.commit()

def get_reciproval_BBH_table(biodb, locus2taxon2best_hit_id, sqlite3=False):

    # taxon_id, taxon_2, seqfeature_id_1, seqfeature_id_2, blast_identity_a_vs_b, blast_identity_b_vs_a, mean_blast_identity, msa_identity, orthogroup

    from chlamdb.biosqldb import manipulate_biosqldb
    import numpy
    from datetime import datetime

    if not sqlite3:
        server, db = manipulate_biosqldb.load_db(biodb)
    else:
        server, db = manipulate_biosqldb.load_db(biodb,
                                                 sqlite=sqlite3)

    sql = 'create table if not exists comparative_tables_reciprocal_BBH (taxon_1 INT, taxon_2 INT, seqfeature_id_1 INT, seqfeature_id_2 INT,' \
          ' blast_evalue_a_vs_b FLOAT, blast_evalue_b_vs_a FLOAT,' \
          ' blast_score_a_vs_b FLOAT, blast_score_b_vs_a FLOAT, blast_identity_a_vs_b FLOAT, blast_identity_b_vs_a FLOAT, ' \
          ' mean_blast_identity FLOAT, msa_identity FLOAT, orthogroup_1 varchar(400), orthogroup_2 varchar(400))' % biodb

    server.adaptor.execute(sql,)
    server.commit()

    sql = 'select locus_tag, taxon_id from annotation_seqfeature_id2locus' % biodb

    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag, orthogroup from orthology_detail' % biodb

    try:
        locus_tag2orthogroup = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    except:
        locus_tag2orthogroup = False
    sql = 'select locus_tag, seqfeature_id from annotation_seqfeature_id2locus' % biodb

    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    if not sqlite3:
        sql = 'select seqfeature_id,char_length(translation) from annotation_seqfeature_id2CDS_annotation' % biodb
    else:
        # sqlite use LENGTH to return the number of characters (return the size in bite in MySQL)
        sql = 'select seqfeature_id,LENGTH(translation) from annotation_seqfeature_id2CDS_annotation' % biodb

    seqfeature_id2sequence_length = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    for n1, locus_1 in enumerate(locus2taxon2best_hit_id):
        if n1 % 50 == 0:
            server.commit()
            print ("%s / %s -- %s " % (n1, len(locus2taxon2best_hit_id),str(datetime.now())))
        taxon_1 = locus_tag2taxon_id[locus_1]

        for taxon_2 in locus2taxon2best_hit_id[locus_1]:
            if taxon_1 == taxon_2:
                # self blast, do not store results
                continue

            ok = False
            locus_2 = locus2taxon2best_hit_id[locus_1][taxon_2][0]
            #try:
            '''
            0 locus_tag_2,
            1 identity,
            2 evalue,
            3 bitscore,
            4 query_align_length,
            5 hit_align_length
            '''

            try:
                ll = locus2taxon2best_hit_id[locus_2][taxon_1][0]
            except:
                continue
            if ll != locus_1:
                # non reciprocal BBH hit
                continue
            else:
                # reciprocal BBH
                identity_a_vs_b = locus2taxon2best_hit_id[locus_1][taxon_2][1]
                identity_b_vs_a = locus2taxon2best_hit_id[locus_2][taxon_1][1]

                evalue_a_vs_b = locus2taxon2best_hit_id[locus_1][taxon_2][2]
                evalue_b_vs_a = locus2taxon2best_hit_id[locus_2][taxon_1][2]

                if (float(evalue_a_vs_b) > 0.00001) or (float(evalue_b_vs_a) > 0.00001):
                    # skip if one has too low evalue
                    continue
                query_align_length_a_vs_b = locus2taxon2best_hit_id[locus_1][taxon_2][4]
                hit_align_length_a_vs_b = locus2taxon2best_hit_id[locus_1][taxon_2][5]
                query_align_length_b_vs_a = locus2taxon2best_hit_id[locus_2][taxon_1][4]
                hit_align_length_b_vs_a = locus2taxon2best_hit_id[locus_2][taxon_1][5]

                seqfeature_id_1 = locus_tag2seqfeature_id[locus_1]
                seqfeature_id_2 = locus_tag2seqfeature_id[locus_2]

                locus_1_length = seqfeature_id2sequence_length[str(seqfeature_id_1)]
                locus_2_length = seqfeature_id2sequence_length[str(seqfeature_id_2)]

                query_cov_a_vs_b = query_align_length_a_vs_b/float(locus_1_length)
                hit_cov_a_vs_b = hit_align_length_a_vs_b/float(locus_2_length)
                query_cov_b_vs_a = query_align_length_b_vs_a/float(locus_2_length)
                hit_cov_b_vs_a = hit_align_length_b_vs_a/float(locus_1_length)

                # keep only if > 50% query and hit coverage
                if (hit_cov_b_vs_a < 0.5) or (hit_cov_a_vs_b < 0.5) or (query_cov_b_vs_a < 0.5) or (query_cov_a_vs_b < 0.5):
                    continue

                score_a_vs_b = locus2taxon2best_hit_id[locus_1][taxon_2][3]
                score_b_vs_a = locus2taxon2best_hit_id[locus_2][taxon_1][3]

                if locus_tag2orthogroup:
                    group_1 = locus_tag2orthogroup[locus_1]
                    group_2 = locus_tag2orthogroup[locus_2]

                    if group_1 != group_2:
                        msa_identity = "NULL"
                    else:
                        sql = 'select identity from orthology_identity where locus_a in ("%s", "%s") ' \
                                ' and locus_b in ("%s","%s") and locus_a!=locus_b;' % (locus_1,
                                                                                       locus_2,
                                                                                       locus_1,
                                                                                       locus_2)

                        msa_identity = server.adaptor.execute_and_fetchall(sql,)[0][0]
                else:
                    msa_identity = 0
                    group_1 = '-'
                    group_2 = '-'

                sql = 'insert into comparative_tables_reciprocal_BBH values (%s,%s,%s,%s,%s,%s,%s,%s' \
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
                    server.adaptor.execute(sql,)
                except:
                    print ('FAIL')
                    print (sql)



            #except:
            #    print 'FAIL!'
            #    print sql
            #    continue
    server.commit()
    sql2 = 'CREATE INDEX taxid1 ON comparative_tables_reciprocal_BBH (taxon_1);' % (biodb, biodb)
    sql3 = 'CREATE INDEX taxid2 ON comparative_tables_reciprocal_BBH (taxon_2);' % (biodb, biodb)
    sql4 = 'CREATE INDEX sfid1 ON comparative_tables_reciprocal_BBH (seqfeature_id_1);' % (biodb, biodb)
    sql5 = 'CREATE INDEX sfid2 ON reciprocal_BBH_%s (seqfeature_id_2);' % (biodb, biodb)
    server.adaptor.execute(sql2,)
    server.adaptor.execute(sql3,)
    server.adaptor.execute(sql4,)
    server.adaptor.execute(sql5,)
    server.commit()


def median_RBBH2species(biodb):

    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table if not exists species_%s (taxon_id INT, species_id INT)' % biodb
    server.adaptor.execute(sql,)
    server.commit()

    sql_taxon = 'select taxon_id from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                ' where t1.name="%s" and t2.description not like "%%%%plasmid%%%%"' % biodb
    taxon_id_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql_taxon,)]

    sql2 = 'select taxon_1,taxon_2,median_identity from comparative_tables_reciprocal_BBH_average_identity;' % biodb

    taxon2taxon2identity = {}
    for row in server.adaptor.execute_and_fetchall(sql2,):
        if row[0] not in taxon2taxon2identity:
            taxon2taxon2identity[row[0]] = {}
            taxon2taxon2identity[row[0]][row[1]] = row[2]
        else:
            taxon2taxon2identity[row[0]][row[1]] = row[2]
        if row[1] not in taxon2taxon2identity:
            taxon2taxon2identity[row[1]] = {}
            taxon2taxon2identity[row[1]][row[0]] = row[2]
        else:
            taxon2taxon2identity[row[1]][row[0]] = row[2]
    #print taxon2taxon2identity
    species_index = 0
    taxons_classified = []
    for taxon_1 in taxon_id_list:
        if taxon_1 in taxons_classified:
            continue
        species_list = [taxon_1]
        for taxon_2 in taxon_id_list:
            if taxon_1 == taxon_2:
                continue
            else:
                if float(taxon2taxon2identity[taxon_1][taxon_2]) >=98:
                    species_list.append(taxon_2)
        for taxon in species_list:
            taxons_classified.append(taxon)
            sql = 'insert into species_%s values(%s,%s)' % (biodb, taxon, species_index)
            server.adaptor.execute(sql,)
        server.commit()
        species_index+=1


def seqfeature_id2n_species_chlamydiae_only(biodb, chlamydiae_taxon_list):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql0 = 'create table if not exists custom_tables.seqfeature_id2n_species_chlamydiae_%s (seqfeature_id INT primary KEY, n_species INT, INDEX n_species(n_species))' % biodb

    server.adaptor.execute(sql0,)

    sql1 ='select locus_tag, orthogroup from orthology_detail' % biodb

    orthogroup2locus_list = {}
    for row in server.adaptor.execute_and_fetchall(sql1,):
        if row[1] not in orthogroup2locus_list:
            orthogroup2locus_list[row[1]] = [row[0]]
        else:
            orthogroup2locus_list[row[1]].append(row[0])

    sql2 = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % biodb
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

    sql3 = 'select taxon_id, species_id from species_%s' % biodb

    taxon_id2species_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3,))

    sql4 = 'select locus_tag, taxon_id from orthology_detail' % biodb

    locus2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql4,))

    for group in orthogroup2locus_list:
        species = []
        locus_keep = []
        for locus in orthogroup2locus_list[group]:
            taxon_id = locus2taxon_id[locus]
            if taxon_id in chlamydiae_taxon_list:
                locus_keep.append(locus)
                species_id = taxon_id2species_id[str(taxon_id)]
                if species_id not in species:
                    species.append(species_id)
        for one_locus in locus_keep:
            seqfeature_id = locus_tag2seqfeature_id[one_locus]
            sql = 'insert into custom_tables.seqfeature_id2n_species_chlamydiae_%s values (%s, %s)' % (biodb, seqfeature_id, len(species))
            server.adaptor.execute(sql,)
        server.commit()

def seqfeature_id2n_species(biodb):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql0 = 'create table if not exists custom_tables_seqfeature_id2n_species (seqfeature_id INT primary KEY, n_species INT, INDEX n_species(n_species))' % biodb

    server.adaptor.execute(sql0,)

    sql1 ='select locus_tag, orthogroup from orthology_detail' % biodb

    orthogroup2locus_list = {}
    for row in server.adaptor.execute_and_fetchall(sql1,):
        if row[1] not in orthogroup2locus_list:
            orthogroup2locus_list[row[1]] = [row[0]]
        else:
            orthogroup2locus_list[row[1]].append(row[0])
    sql2 = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % biodb
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

    sql3 = 'select taxon_id, species_id from species_%s' % biodb

    taxon_id2species_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3,))

    sql4 = 'select locus_tag, taxon_id from orthology_detail' % biodb

    locus2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql4,))

    for group in orthogroup2locus_list:
        species = list(set([taxon_id2species_id[str(locus2taxon_id[i])] for i in orthogroup2locus_list[group]]))
        print (group, 'n species:', len(species))
        for one_locus in orthogroup2locus_list[group]:
            seqfeature_id = locus_tag2seqfeature_id[one_locus]
            sql = 'insert into custom_tables_seqfeature_id2n_species values (%s, %s)' % (biodb, seqfeature_id, len(species))
            server.adaptor.execute(sql,)
        server.commit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-s",'--seq_id_file',type=str,help="SequenceIDs.txt")
    parser.add_argument("-b",'--blast_files', nargs="+", help="blast_file")
    parser.add_argument("-d",'--database', help="bioadatabase name")

    args = parser.parse_args()

    '''
    taxid_chlam_list = [314,
                        886707,
                        804807,
                            48,
                        1143376,
                            283,
                            55,
                        1279839,
                        1279496,
                        1279822,
                            46,
                        1279815,
                            49,
                            87925,
                            52,
                        1137444,
                            67,
                        1172028,
                        1069693,
                        1172027,
                            307,
                            59,
                            60,
                        1279767,
                        1035343,
                            313,
                            62,
                            293,
                        1069694,
                        1279497,
                        1279774,
                            64,
                            66]

    seqfeature_id2n_species_chlamydiae_only(args.database, taxid_chlam_list)
    '''
    print("Parsing SeqID file...")
    seqid2locus_tag = parse_seqid_file(args.seq_id_file)
    print ('Parsing blast files...')
    locus2taxon2best_hit_id = parse_orthofinder_blast_files(args.database,
                                                            args.blast_files,
                                                            seqid2locus_tag)
    print('Loading data...')
    get_reciproval_BBH_table(args.database,
                             locus2taxon2best_hit_id)

    print ('Get pairwise median id table...')
    get_parwise_genome_median_identity_table(args.database)

    print("Median_RBBH2species...")
    median_RBBH2species(args.database)

    print("seqfeature_id2n_species...")
    seqfeature_id2n_species(args.database)
