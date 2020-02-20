#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

def get_best_non_top_phylum_hit(biodb):

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select taxon_id, count(*) as n, best_hit_phylum from blastnr_blastnr_majority_phylum' \
          ' group by taxon_id,best_hit_phylum order by taxon_id,n DESC;' % biodb

    data = server.adaptor.execute_and_fetchall(sql,)

    taxon2most_freq_phylum = {}
    for row in data:
        if row[0] not in taxon2most_freq_phylum:
            taxon2most_freq_phylum[row[0]] = row[2]

    sql = 'create table blastnr_blastnr_best_non_self_phylum (query_taxon_id INT, ' \
          ' seqfeature_id INT, ' \
          ' hit_number INT, ' \
          ' superkingdom varchar(200),' \
          ' subject_taxon_id INT,' \
          ' percent_identity FLOAT,' \
          ' subject_accession varchar(200),' \
          ' subject_title TEXT,' \
          ' index query_taxon_id(query_taxon_id),' \
          ' index seqfeature_id(seqfeature_id),' \
          ' index subject_taxon_id(subject_taxon_id),' \
          ' index superkingdom(superkingdom))' % (biodb)

    server.adaptor.execute(sql,)
    server.commit()

    # iter taxons
    for taxon in taxon2most_freq_phylum:

        print 'taxon', taxon
        # get all non top phylum hits
        phylum_filter = taxon2most_freq_phylum[taxon]

        sql = 'select t1.query_taxon_id, t1.seqfeature_id, t1.hit_number, t2.superkingdom, t1.subject_scientific_name, ' \
              ' t1.percent_identity, t2.taxon_id,kingdom,phylum,class,subject_title, subject_accession from blastnr_blastnr t1 ' \
              ' inner join blastnr_blastnr_taxonomy t2 on t1.subject_taxid=t2.taxon_id ' \
              ' where (t1.query_taxon_id=%s and t2.phylum != "%s") order by seqfeature_id' % (biodb,
                                                                                        taxon,
                                                                                        phylum_filter)

        data = server.adaptor.execute_and_fetchall(sql,)

        # keep best non top phylum hit for each feature
        for n, row in enumerate(data):

            query_taxon_id = row[0]
            seqfeature_id = row[1]
            hit_number = row[2]
            superkingdom = row[3]
            subject_taxon_id = row[6]
            percent_identity = row[5]
            subject_accession = row[10]
            subject_title = row[11]
            if n%10000 == 0:
                print "%s / %s" % (n, len(data))
            if n == 0:

                sql = 'insert into blastnr_blastnr_best_non_self_phylum values (%s, %s, %s, "%s", %s, %s, "%s", "%s")' % (biodb,
                                                                                                                 query_taxon_id,
                                                                                                                seqfeature_id,
                                                                                                                hit_number,
                                                                                                                superkingdom,
                                                                                                                subject_taxon_id,
                                                                                                                percent_identity,
                                                                                                                subject_accession,
                                                                                                                subject_title)
                print sql
                server.adaptor.execute(sql,)
                server.commit()

            else:
                # if new feature
                if row[1] != data[n-1][1]:
                    sql = 'insert into blastnr_blastnr_best_non_self_phylum values (%s, %s, %s, "%s", %s, %s, "%s", "%s")' % (biodb,
                                                                                                                     query_taxon_id,
                                                                                                                    seqfeature_id,
                                                                                                                    hit_number,
                                                                                                                    superkingdom,
                                                                                                                    subject_taxon_id,
                                                                                                                    percent_identity,
                                                                                                                    subject_accession,
                                                                                                                    subject_title)
                    server.adaptor.execute(sql,)
                    server.commit()


def best_blast_hit_majority_species(biodb):

    '''

    1. for each protein, get the species of the best hit, and calculate the proportion of the
    proteome it represent

    :param biodb:
    :param n_hits:
    :return:
    '''


    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select taxon_id from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id' \
          ' where t1.name="%s" group by taxon_id' % biodb

    taxon_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    sql = 'CREATE table blastnr.blastnr_best_hits_species_%s (taxon_id INT, hit_taxid INT, count INT, proportion FLOAT)' % biodb

    server.adaptor.execute(sql,)
    server.commit()

    sql = 'select taxon_id, count(*) as n from orthology_detail group by taxon_id' % (biodb)
    taxon_id2n_CDS = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    for one_taxon in taxon_list:
        print one_taxon

        sql2 = 'select A.subject_taxid, count(*) from blastnr_blastnr A ' \
               ' where A.query_taxon_id=%s and hit_number=1 group by A.subject_taxid;' % (biodb, one_taxon)

        taxon_id2count = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
        #print seqfeature2best_hit_phylum
        total_bbh = float(sum([int(i) for i in taxon_id2count.values()]))

        for n, taxon in enumerate(taxon_id2count):
                sql = 'insert into blastnr.blastnr_best_hits_species_%s values(%s,%s,%s, %s)' % (biodb,
                                                                                                 one_taxon,
                                                                                                 taxon,
                                                                                                 taxon_id2count[taxon],
                                                                                                 round((taxon_id2count[taxon]/total_bbh)*100,2))


                server.adaptor.execute(sql,)


        server.commit()


def majority_phylum(biodb, n_hits):

    '''

    1. for each protein, get the most frequent phylum among the top n_hits (typically 100)
    2. for each protein, get the phylum of the best hit

    :param biodb:
    :param n_hits:
    :return:
    '''


    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select taxon_id from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id' \
                   ' where t1.name="%s" group by taxon_id' % biodb

    taxon_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    print taxon_list
    sql = 'CREATE table blastnr_blastnr_majority_phylum (taxon_id INT, seqfeature_id INT, best_hit_phylum varchar(400), majority_phylum varchar(400))' % biodb

    server.adaptor.execute(sql,)
    server.commit()

    for one_taxon in taxon_list:
        print one_taxon
        sql = 'select seqfeature_id,phylum, count(*) as n from blastnr_blastnr A ' \
              ' inner join blastnr_blastnr_taxonomy B on A.subject_taxid=B.taxon_id ' \
              ' where A.query_taxon_id=%s and hit_number<=%s group by seqfeature_id,phylum order by seqfeature_id, n DESC;' % (biodb,
                                                                                                                       one_taxon,
                                                                                                                       n_hits)
        sql2 = 'select seqfeature_id,phylum from blastnr_blastnr A ' \
               ' inner join blastnr_blastnr_taxonomy B on A.subject_taxid=B.taxon_id  ' \
               ' where A.query_taxon_id=%s and hit_number=1;' % (biodb, one_taxon)

        print sql2
        data_majority = server.adaptor.execute_and_fetchall(sql,)
        seqfeature2best_hit_phylum = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
        #print seqfeature2best_hit_phylum
        seqfeature2phylum = {}
        for n, row in enumerate(data_majority):
            if row[0] not in seqfeature2phylum:
                try:
                    if row[0] == data_majority[n+1][0]:
                        if row[2] == data_majority[n+1][0]:
                            print 'equality!!!!!!!!!!!!!!!!!!!!'
                        seqfeature2phylum[row[0]] = row[1]
                        sql = 'insert into blastnr_blastnr_majority_phylum values(%s,%s,"%s", "%s")' % (biodb,
                                                                                                           one_taxon,
                                                                                                           row[0],
                                                                                                           seqfeature2best_hit_phylum[str(row[0])],
                                                                                                           row[1])
                except (IndexError, KeyError) as e:
                    if isinstance(e, IndexError):
                        continue
                    else:
                        print 'error!!!!!, missing BEST BLAST HIT!', sql

                server.adaptor.execute(sql,)


        server.commit()

def count_less_than_n_hits(biodb, cutoff=100):

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    # n 0 hits
    # n 1-99 hits
    # n 100 hits

    sql = 'CREATE TABLE blastnr_count_n_blast (taxon_id INT, total INT, n_no_hits INT, n_100_hits INT, n_less_100_hits INT)' % biodb

    server.adaptor.execute(sql)

    sql = 'select query_taxon_id, seqfeature_id, count(*) as n from blastnr_blastnr group by query_taxon_id, seqfeature_id;' % biodb
    taxon2n_100 = {}
    taxon2n_less_than_100 = {}

    data = server.adaptor.execute_and_fetchall(sql,)
    for row in data:
        if row[0] not in taxon2n_100:
            taxon2n_100[row[0]] = 0
            taxon2n_less_than_100[row[0]] = 0
        if row[2]>=cutoff:
            taxon2n_100[row[0]]+=1
        else:
            taxon2n_less_than_100[row[0]]+=1

    sql = 'select taxon_id, count(*) from orthology_detail group by taxon_id' % biodb

    taxon2n_proteins = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    for taxon in taxon2n_proteins:
        total = taxon2n_proteins[taxon]
        n_100 = taxon2n_100[int(taxon)]
        n_less_than_100 = taxon2n_less_than_100[int(taxon)]
        sql = 'insert into blastnr_count_n_blast values (%s,%s,%s,%s, %s)' % (biodb,
                                                                               taxon,
                                                                               total,
                                                                               total-n_100-n_less_than_100,
                                                                               n_100,
                                                                               n_less_than_100)
        server.adaptor.execute(sql,)
    server.commit()

def count_majority_phylum(biodb, hit_number=1):

    # n Chlamydiae (or any top ranking phylum)
    # n non (or any top ranking phylum)
    # n eukaryotes
    # n viruses
    # n archae
    # n undefined

    import biosql_own_sql_tables as bsql
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table blastnr.BBH_taxo_hit_number_%s_%s (taxon_id INT, n_top_phylum INT, ' \
          ' n_non_top_phylum INT, n_euk INT, n_arch INT, n_virus INT, undefined INT)' % (hit_number,biodb)
    print sql
    server.adaptor.execute(sql,)
    server.commit()

    sql = 'select taxon_id, count(*) as n, best_hit_phylum from blastnr_blastnr_majority_phylum' \
          ' group by taxon_id,best_hit_phylum order by taxon_id,n DESC;' % biodb
    print sql
    data = server.adaptor.execute_and_fetchall(sql,)
    taxon2most_freq_phylum = {}
    for row in data:
        if row[0] not in taxon2most_freq_phylum:
            taxon2most_freq_phylum[row[0]] = row[2]

    sql = 'select taxon_id from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id' \
          ' where t1.name="%s" group by taxon_id' % biodb

    taxon_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]


    for taxon_id in taxon_list:
        print 'taxon id, most freq!', taxon_id, taxon2most_freq_phylum[taxon_id]

        # for best hit and for majority rule
        sql_bacteria_chlamydiae = 'select t1.query_taxon_id, count(*) as n from blastnr_blastnr t1 inner join blastnr_blastnr_taxonomy t2 ' \
                                  ' on t1.subject_taxid=t2.taxon_id ' \
                                  ' where t1.hit_number=%s and t1.query_taxon_id=%s and t2.superkingdom="Bacteria" ' \
                                  ' and t2.phylum="%s" group by t1.query_taxon_id;' % (biodb,
                                                                                       hit_number,
                                                                                       taxon_id,
                                                                                       taxon2most_freq_phylum[taxon_id]
                                                                                       )

        sql_bacteria_non_chlamydia = 'select t1.query_taxon_id, count(*) as n from blastnr_blastnr t1 inner join blastnr_blastnr_taxonomy t2 ' \
                                  ' on t1.subject_taxid=t2.taxon_id ' \
                                  ' where t1.hit_number=%s and t1.query_taxon_id=%s and t2.superkingdom="Bacteria" ' \
                                  ' and t2.phylum!="%s" group by t1.query_taxon_id;' % (biodb,
                                                                                        hit_number,
                                                                                        taxon_id,
                                                                                        taxon2most_freq_phylum[taxon_id])

        sql_euk = 'select t1.query_taxon_id, count(*) as n from blastnr_blastnr t1 inner join blastnr_blastnr_taxonomy t2 ' \
                  ' on t1.subject_taxid=t2.taxon_id where t1.hit_number=%s and t1.query_taxon_id=%s' \
                  ' and t2.superkingdom="Eukaryota" group by t1.query_taxon_id;' % (biodb,
                                                                                    hit_number,
                                                                                    taxon_id)

        sql_archae = 'select t1.query_taxon_id, count(*) as n from blastnr_blastnr t1 inner join blastnr_blastnr_taxonomy t2 ' \
                  ' on t1.subject_taxid=t2.taxon_id where t1.hit_number=%s and t1.query_taxon_id=%s' \
                  ' and t2.superkingdom="Archaea" group by t1.query_taxon_id;' % (biodb,
                                                                                  hit_number,
                                                                                  taxon_id)

        sql_viruses = 'select t1.query_taxon_id, count(*) as n from blastnr_blastnr t1 inner join blastnr_blastnr_taxonomy t2 ' \
                  ' on t1.subject_taxid=t2.taxon_id where t1.hit_number=%s and t1.query_taxon_id=%s' \
                  ' and t2.superkingdom="Viruses" group by t1.query_taxon_id;' % (biodb,
                                                                                  hit_number,
                                                                                  taxon_id)

        sql_undefined = 'select t1.query_taxon_id, count(*) as n from blastnr_blastnr t1 inner join blastnr_blastnr_taxonomy t2 ' \
                  ' on t1.subject_taxid=t2.taxon_id where t1.hit_number=%s and t1.query_taxon_id=%s' \
                  ' and t2.superkingdom="-" group by t1.query_taxon_id;' % (biodb,
                                                                            hit_number,
                                                                            taxon_id)

        taxon2n_chlamydiae = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_bacteria_chlamydiae))
        taxon2n_non_chlamydiae = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_bacteria_non_chlamydia))
        taxon2n_euk = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_euk))
        taxon2n_arch = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_archae))
        taxon2n_virus = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_viruses))
        taxon2undefined = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_undefined))

        try:
            non_phylum = taxon2n_non_chlamydiae[str(taxon_id)]
        except KeyError:
            non_phylum = 0
        try:
            euk = taxon2n_euk[str(taxon_id)]
        except KeyError:
            euk= 0
        try:
            archae = taxon2n_arch[str(taxon_id)]
        except KeyError:
            archae = 0
        try:
            virus = taxon2n_virus[str(taxon_id)]
        except KeyError:
            virus = 0
        try:
            undef = taxon2undefined[str(taxon_id)]
        except KeyError:
            undef = 0


        sql = 'insert into blastnr.BBH_taxo_hit_number_%s_%s values (%s,%s,%s,%s,%s,%s,%s)' % (hit_number,
                                                                                               biodb,
                                                                                               taxon_id,
                                                                                               taxon2n_chlamydiae[str(taxon_id)],
                                                                                               non_phylum,
                                                                                               euk,
                                                                                               archae,
                                                                                               virus,
                                                                                               undef)
        server.adaptor.execute(sql,)
        server.commit()


def count_majority_excluding_self_species(biodb):
    # taxonomy of the best non identical hit

    # n Chlamydiae (or any top ranking phylum)
    # n non (or any top ranking phylum)
    # n eukaryotes
    # n viruses
    # n archae
    # n undefined

    '''


    select seqfeature_id,hit_number,subject_taxid from
    (select * from blastnr_2017_06_29b_motile_chlamydiae where percent_identity!=100 group by seqfeature_id) A;
    '''


    import biosql_own_sql_tables as bsql
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table blastnr.BBH_taxo_hit_non_identical_%s (taxon_id INT, n_top_phylum INT, ' \
          ' n_non_top_phylum INT, n_euk INT, n_arch INT, n_virus INT, undefined INT)' % (biodb)
    print sql
    server.adaptor.execute(sql,)
    server.commit()

    sql = 'select taxon_id, count(*) as n, best_hit_phylum from blastnr_blastnr_majority_phylum' \
          ' group by taxon_id,best_hit_phylum order by taxon_id,n DESC;' % biodb

    data = server.adaptor.execute_and_fetchall(sql,)
    taxon2most_freq_phylum = {}
    for row in data:
        if row[0] not in taxon2most_freq_phylum:
            taxon2most_freq_phylum[row[0]] = row[2]

    sql = 'select taxon_id from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id' \
          ' where t1.name="%s" group by taxon_id' % biodb

    taxon_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]


    for taxon_id in taxon_list:
        print 'taxon id, most freq!', taxon_id, taxon2most_freq_phylum[taxon_id]

        # for best hit and for majority rule
        sql_bacteria_chlamydiae = 'select query_taxon_id,count(*) from ' \
                                  ' (select query_taxon_id,seqfeature_id,hit_number,subject_taxid from ' \
                                  ' (select * from blastnr_blastnr ' \
                                  ' where percent_identity!=100 and query_taxon_id=%s ' \
                                  ' group by seqfeature_id) A )B ' \
                                  ' inner join blastnr_blastnr_taxonomy C on B.subject_taxid=C.taxon_id ' \
                                  ' where C.phylum="%s" group by query_taxon_id;' % (biodb,
                                                                                       taxon_id,
                                                                                       taxon2most_freq_phylum[taxon_id]
                                                                                       )

        sql_bacteria_non_chlamydia = 'select query_taxon_id,count(*) from ' \
                                  ' (select query_taxon_id,seqfeature_id,hit_number,subject_taxid from ' \
                                  ' (select * from blastnr_blastnr ' \
                                  ' where percent_identity!=100 and query_taxon_id=%s ' \
                                  ' group by seqfeature_id) A )B ' \
                                  ' inner join blastnr_blastnr_taxonomy C on B.subject_taxid=C.taxon_id ' \
                                  ' where C.phylum!="%s" and C.superkingdom="Bacteria" group by query_taxon_id;' % (biodb,
                                                                                        taxon_id,
                                                                                        taxon2most_freq_phylum[taxon_id])

        sql_euk = 'select query_taxon_id,count(*) from ' \
                                  ' (select query_taxon_id,seqfeature_id,hit_number,subject_taxid from ' \
                                  ' (select * from blastnr_blastnr ' \
                                  ' where percent_identity!=100 and query_taxon_id=%s ' \
                                  ' group by seqfeature_id) A )B ' \
                                  ' inner join blastnr_blastnr_taxonomy C on B.subject_taxid=C.taxon_id ' \
                                  ' where C.superkingdom="Eukaryota" group by query_taxon_id;' % (biodb,
                                                                                    taxon_id)

        sql_archae = 'select query_taxon_id,count(*) from ' \
                                  ' (select query_taxon_id,seqfeature_id,hit_number,subject_taxid from ' \
                                  ' (select * from blastnr_blastnr ' \
                                  ' where percent_identity!=100 and query_taxon_id=%s ' \
                                  ' group by seqfeature_id) A )B ' \
                                  ' inner join blastnr_blastnr_taxonomy C on B.subject_taxid=C.taxon_id ' \
                                  ' where C.superkingdom="Archaea" group by query_taxon_id;' % (biodb,
                                                                                    taxon_id)

        sql_viruses = 'select query_taxon_id,count(*) from ' \
                                  ' (select query_taxon_id,seqfeature_id,hit_number,subject_taxid from ' \
                                  ' (select * from blastnr_blastnr ' \
                                  ' where percent_identity!=100 and query_taxon_id=%s ' \
                                  ' group by seqfeature_id) A )B ' \
                                  ' inner join blastnr_blastnr_taxonomy C on B.subject_taxid=C.taxon_id ' \
                                  ' where C.superkingdom="Viruses" group by query_taxon_id;' % (biodb,
                                                                                  taxon_id)

        sql_undefined = 'select query_taxon_id,count(*) from ' \
                                  ' (select query_taxon_id,seqfeature_id,hit_number,subject_taxid from ' \
                                  ' (select * from blastnr_blastnr ' \
                                  ' where percent_identity!=100 and query_taxon_id=%s ' \
                                  ' group by seqfeature_id) A )B ' \
                                  ' inner join blastnr_blastnr_taxonomy C on B.subject_taxid=C.taxon_id ' \
                                  ' where C.superkingdom="-" group by query_taxon_id;' % (biodb,
                                                                                  taxon_id)

        taxon2n_chlamydiae = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_bacteria_chlamydiae))
        taxon2n_non_chlamydiae = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_bacteria_non_chlamydia))
        taxon2n_euk = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_euk))
        taxon2n_arch = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_archae))
        taxon2n_virus = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_viruses))
        taxon2undefined = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_undefined))

        try:
            non_phylum = taxon2n_non_chlamydiae[str(taxon_id)]
        except KeyError:
            non_phylum = 0
        try:
            euk = taxon2n_euk[str(taxon_id)]
        except KeyError:
            euk= 0
        try:
            archae = taxon2n_arch[str(taxon_id)]
        except KeyError:
            archae = 0
        try:
            virus = taxon2n_virus[str(taxon_id)]
        except KeyError:
            virus = 0
        try:
            undef = taxon2undefined[str(taxon_id)]
        except KeyError:
            undef = 0


        sql = 'insert into blastnr.BBH_taxo_hit_non_identical_%s values (%s,%s,%s,%s,%s,%s,%s)' % (                                                                                            biodb,
                                                                                               taxon_id,
                                                                                               taxon2n_chlamydiae[str(taxon_id)],
                                                                                               non_phylum,
                                                                                               euk,
                                                                                               archae,
                                                                                               virus,
                                                                                               undef)
        server.adaptor.execute(sql,)
        server.commit()

def count_majority_phylum_non_identical(biodb, accession2expluded_taxon_id):

    # taxonomy of the best non identical hit

    # n Chlamydiae (or any top ranking phylum)
    # n non (or any top ranking phylum)
    # n eukaryotes
    # n viruses
    # n archae
    # n undefined

    '''
 select seqfeature_id,hit_number,subject_taxid from
 (select * from blastnr_2017_06_29b_motile_chlamydiae where percent_identity!=100 group by seqfeature_id) A;
    '''

    import biosql_own_sql_tables as bsql
    import taxon_id2child
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table blastnr.BBH_taxo_hit_non_identical_%s (taxon_id INT, n_top_phylum INT, ' \
          ' n_non_top_phylum INT, n_euk INT, n_arch INT, n_virus INT, undefined INT)' % (biodb)
    print sql
    server.adaptor.execute(sql,)
    server.commit()

    sql = 'select taxon_id, accession from biosqldb.bioentry t1 inner join biosqldb.biodatabase t2' \
          ' on t1.biodatabase_id=t2.biodatabase_id where t2.name="%s" and t1.description not like "%%%%plasmid%%%%"' % biodb

    taxon_id2genome_accession = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql),)
    print taxon_id2genome_accession

    sql = 'select taxon_id, count(*) as n, best_hit_phylum from blastnr_blastnr_majority_phylum' \
          ' group by taxon_id,best_hit_phylum order by taxon_id,n DESC;' % biodb

    data = server.adaptor.execute_and_fetchall(sql,)
    taxon2most_freq_phylum = {}
    for row in data:
        if row[0] not in taxon2most_freq_phylum:
            taxon2most_freq_phylum[row[0]] = row[2]

    sql = 'select taxon_id from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id' \
          ' where t1.name="%s" group by taxon_id' % biodb

    taxon_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]


    for taxon_id in taxon_list:
        print 'taxon id, most freq!', taxon_id, taxon2most_freq_phylum[taxon_id]
        try:
            excluded_taxid =str( accession2expluded_taxon_id[taxon_id2genome_accession[str(taxon_id)]])
            excluded_taxon_list = taxon_id2child.taxon_id2child(excluded_taxid)
            print 'excluded taxon list:', excluded_taxon_list
            excluded_taxid_filter = 'and subject_taxid not in (%s)' % ','.join(excluded_taxon_list)
            print excluded_taxid_filter
        except:
            excluded_taxid_filter = ''


            # for best hit and for majority rule
        sql_bacteria_chlamydiae = 'select query_taxon_id,count(*) from ' \
                                  ' (select query_taxon_id,seqfeature_id,hit_number,subject_taxid from ' \
                                  ' (select * from blastnr_blastnr ' \
                                  ' where percent_identity!=0 and query_taxon_id=%s ' \
                                  ' %s group by seqfeature_id) A )B ' \
                                  ' inner join blastnr_blastnr_taxonomy C on B.subject_taxid=C.taxon_id ' \
                                  ' where C.phylum="%s" group by query_taxon_id;' % (biodb,
                                                                                     taxon_id,
                                                                                     excluded_taxid_filter,
                                                                                     taxon2most_freq_phylum[taxon_id]
                                                                                     )

        sql_bacteria_non_chlamydia = 'select query_taxon_id,count(*) from ' \
                                  ' (select query_taxon_id,seqfeature_id,hit_number,subject_taxid from ' \
                                  ' (select * from blastnr_blastnr ' \
                                  ' where percent_identity!=0 and query_taxon_id=%s ' \
                                  ' %s group by seqfeature_id) A )B ' \
                                  ' inner join blastnr_blastnr_taxonomy C on B.subject_taxid=C.taxon_id ' \
                                  ' where C.phylum!="%s" and C.superkingdom="Bacteria" group by query_taxon_id;' % (biodb,
                                                                                        taxon_id,
                                                                                        excluded_taxid_filter,
                                                                                        taxon2most_freq_phylum[taxon_id])

        sql_euk = 'select query_taxon_id,count(*) from ' \
                                  ' (select query_taxon_id,seqfeature_id,hit_number,subject_taxid from ' \
                                  ' (select * from blastnr_blastnr ' \
                                  ' where percent_identity!=0 and query_taxon_id=%s ' \
                                  ' %s group by seqfeature_id) A )B ' \
                                  ' inner join blastnr_blastnr_taxonomy C on B.subject_taxid=C.taxon_id ' \
                                  ' where C.superkingdom="Eukaryota" group by query_taxon_id;' % (biodb,
                                                                                    taxon_id,
                                                                                    excluded_taxid_filter,)

        sql_archae = 'select query_taxon_id,count(*) from ' \
                                  ' (select query_taxon_id,seqfeature_id,hit_number,subject_taxid from ' \
                                  ' (select * from blastnr_blastnr ' \
                                  ' where percent_identity!=0 and query_taxon_id=%s ' \
                                  ' %s group by seqfeature_id) A )B ' \
                                  ' inner join blastnr_blastnr_taxonomy C on B.subject_taxid=C.taxon_id ' \
                                  ' where C.superkingdom="Archaea" group by query_taxon_id;' % (biodb,
                                                                                                taxon_id,
                                                                                                excluded_taxid_filter,)

        sql_viruses = 'select query_taxon_id,count(*) from ' \
                                  ' (select query_taxon_id,seqfeature_id,hit_number,subject_taxid from ' \
                                  ' (select * from blastnr_blastnr ' \
                                  ' where percent_identity!=0 and query_taxon_id=%s ' \
                                  ' %s group by seqfeature_id) A )B ' \
                                  ' inner join blastnr_blastnr_taxonomy C on B.subject_taxid=C.taxon_id ' \
                                  ' where C.superkingdom="Viruses" group by query_taxon_id;' % (biodb,
                                                                                                taxon_id,
                                                                                                excluded_taxid_filter)

        sql_undefined = 'select query_taxon_id,count(*) from ' \
                                  ' (select query_taxon_id,seqfeature_id,hit_number,subject_taxid from ' \
                                  ' (select * from blastnr_blastnr ' \
                                  ' where percent_identity!=0 and query_taxon_id=%s ' \
                                  ' %s group by seqfeature_id) A )B ' \
                                  ' inner join blastnr_blastnr_taxonomy C on B.subject_taxid=C.taxon_id ' \
                                  ' where C.superkingdom="-" group by query_taxon_id;' % (biodb,
                                                                                          taxon_id,
                                                                                          excluded_taxid_filter)

        print 'sql_bacteria_chlamydiae', sql_bacteria_chlamydiae
        taxon2n_chlamydiae = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_bacteria_chlamydiae))
        print 'sql_bacteria_non_chlamydia',sql_bacteria_non_chlamydia
        taxon2n_non_chlamydiae = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_bacteria_non_chlamydia))
        print 'sql_euk',sql_euk
        taxon2n_euk = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_euk))
        print 'sql_archae',sql_archae
        taxon2n_arch = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_archae))
        print 'sql_viruses', sql_viruses
        taxon2n_virus = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_viruses))
        print 'sql_undefined', sql_undefined
        taxon2undefined = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_undefined))

        try:
            non_phylum = taxon2n_non_chlamydiae[str(taxon_id)]
        except KeyError:
            non_phylum = 0
        try:
            euk = taxon2n_euk[str(taxon_id)]
        except KeyError:
            euk= 0
        try:
            archae = taxon2n_arch[str(taxon_id)]
        except KeyError:
            archae = 0
        try:
            virus = taxon2n_virus[str(taxon_id)]
        except KeyError:
            virus = 0
        try:
            undef = taxon2undefined[str(taxon_id)]
        except KeyError:
            undef = 0


        sql = 'insert into blastnr.BBH_taxo_hit_non_identical_%s values (%s,%s,%s,%s,%s,%s,%s)' % (                                                                                            biodb,
                                                                                               taxon_id,
                                                                                               taxon2n_chlamydiae[str(taxon_id)],
                                                                                               non_phylum,
                                                                                               euk,
                                                                                               archae,
                                                                                               virus,
                                                                                               undef)
        server.adaptor.execute(sql,)
        server.commit()

def count_majority_phylum_consensus(biodb):
    # ATTENTION PAS FINI!
    # claculate majority rule top n hits
    # then
    # n Chlamydiae (or any top ranking phylum)
    # n non (or any top ranking phylum)
    # n eukaryotes
    # n viruses
    # n archae
    # n undefined

    import biosql_own_sql_tables as bsql
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table blastnr.BBH_taxo_%s (taxon_id INT, n_top_phylum INT, ' \
          ' n_non_top_phylum INT, n_euk INT, n_arch INT, n_virus INT, undefined INT)' % (biodb)
    print sql
    server.adaptor.execute(sql,)
    server.commit()

    sql = 'select taxon_id, count(*) as n, best_hit_phylum from blastnr_blastnr_majority_phylum' \
          ' group by taxon_id,best_hit_phylum order by taxon_id,n DESC;' % biodb
    print sql
    data = server.adaptor.execute_and_fetchall(sql,)
    taxon2most_freq_phylum = {}
    for row in data:
        if row[0] not in taxon2most_freq_phylum:
            taxon2most_freq_phylum[row[0]] = row[2]

    sql = 'select taxon_id from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id' \
          ' where t1.name="%s" group by taxon_id' % biodb

    taxon_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]


    for taxon_id in taxon_list:
        print 'taxon id, most freq!', taxon_id, taxon2most_freq_phylum[taxon_id]

        # for best hit and for majority rule
        sql_phylum = 'select t1.seqfeature_id, t2.phylum, count(*) as n from blastnr_blastnr t1 inner join blastnr_blastnr_taxonomy t2 ' \
                                  ' on t1.subject_taxid=t2.taxon_id ' \
                                  ' where t1.query_taxon_id=%s and t2.superkingdom="Bacteria" ' \
                                  ' and t2.phylum="%s" group by t1.seqfeature_id,t2.phylum order by t1.seqfeature_id,n DESC;' % (biodb,
                                                                                       taxon_id,
                                                                                       taxon2most_freq_phylum[taxon_id]
                                                                                       )
        print sql_phylum
        sql_superkingdom = 'select t1.seqfeature_id, t2.superkingdom, count(*) as n from blastnr_blastnr t1 inner join blastnr_blastnr_taxonomy t2 ' \
                                  ' on t1.subject_taxid=t2.taxon_id ' \
                                  ' where t1.hit_number=1 and t1.query_taxon_id=%s and t2.superkingdom!="Bacteria" ' \
                                  ' and t2.phylum!="%s" group by t1.query_taxon_id;' % (biodb,
                                                                                        taxon_id,
                                                                                        taxon2most_freq_phylum[taxon_id])

        n_chlam = 0
        n_non_chlamydiae = 0
        n_euk = 0
        n_arch = 0
        n_virus = 0
        n_undefined = 0


        sql = 'insert into blastnr.blast_hits_taxonomy_overview_%s values (%s,%s,%s,%s,%s,%s,%s)' % (biodb,
                                                                                             taxon_id,
                                                                                             n_chlam,
                                                                                             n_non_chlamydiae,
                                                                                             n_euk,
                                                                                             n_arch,
                                                                                             n_virus,
                                                                                             n_undefined)
        server.adaptor.execute(sql,)
        server.commit()

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--biodb', type=str, help="biodatabase name")

    args = parser.parse_args()

    #majority_phylum(args.biodb, 100)
    #count_less_than_n_hits(args.biodb)
    #count_majority_phylum(args.biodb, 1)
    #count_majority_phylum(args.biodb, 2)
    #get_best_non_top_phylum_hit(args.biodb)

    #count_majority_phylum_non_identical(args.biodb)
    #best_blast_hit_majority_species(args.biodb)


    accession2excluded_taxon = {
 "CWGJ01000001":    	483423,
 "CCSC01000000":     	1444712,
 "LJUH01000000":     	1704235,
 "MKSG01000000":     	1895742,
 "MKSK01000000 ":    	1895743,
 "NC_002620":        	83560,
 "NC_003361":        	83557,
 "NC_000117":        	813,
 "NC_000922":        	83558,
 "NC_007899":        	83556,
 "NC_005861":        	362787,
 "NC_010655":        	239935,
 "NC_004552":        	83555,
 "NC_014225":        	71667,
 "NZ_ACZE00000000":  	83552,
 "NC_015470":        	83554,
 "NC_015408":        	85991,
 "NC_015713":        	83561,
 "NC_015702":        	83552,
 "NZ_APJW00000000":  	1405396,
 "NZ_CP015840 ":     	1457153,
 "NZ_AYKJ01000000":  	83559,
 "NZ_CP006571 ":     	1457141,
 "NZ_BASK00000000":  	1353976,
 "NZ_BASL00000000":  	1353977,
 "NZ_CCEJ000000000": 	340071,
 "NZ_CCJF00000000":  	1444711,
 "NZ_JSAM00000000":  	83552,
 "NZ_JSAN00000000":  	362787,
 "NZ_JRXI00000000":  	1478174,
 "NZ_JSDQ00000000":  	1478175,
 "NZ_BBPT00000000":  	1314958,
 "NZ_BAWW00000000":  	83552,
 "NZ_LN879502":      	389348,
 "NZ_FCNU00000000":  	1414722,
 "NZ_CP014639":      	1806891,
 "NZ_BCPZ00000000":  	1785087,
 "NZ_FLYO00000000":  	1871324,
 "NZ_FLYF00000000":  	1871322,
 "NZ_FLYP00000000":  	1871323,
 "LNES01000000":     	1752208,
 "MGLO01000000":     	1797595,
 "MGLP01000000":     	1797596,
 "MGLQ01000000":     	1797597,
 "MGLR01000000":     	1797598,
 "MGLS01000000":     	1797599,
 "MGLT01000000":     	1797600,
 "MGLU01000000":     	1797601,
 "MGLV01000000":     	1797602,
 "MGLW01000000":     	1797603,
 "MGLX01000000":     	1797604,
 "MGLY01000000":     	1797605,
 "MGLZ01000000":     	1797606,
 "MGMA01000000":     	1797607,
 "MGMB01000000":     	1797608,
 "MGMC01000000":     	1797609,
 "MGMD01000000":     	1797610,
 "MGME01000000":     	1797611}

    count_majority_phylum_non_identical(args.biodb, accession2excluded_taxon)