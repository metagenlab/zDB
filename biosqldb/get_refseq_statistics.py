#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


def majority_phylum(biodb, n_hits):
    import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select taxon_id from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id' \
                   ' where t1.name="%s" group by taxon_id' % biodb

    taxon_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    print taxon_list
    sql = 'CREATE table blastnr.blastnr_majority_phylum_%s (taxon_id INT, seqfeature_id INT, best_hit_phylum varchar(400), majority_phylum varchar(400))' % biodb

    server.adaptor.execute(sql,)
    server.commit()

    for one_taxon in taxon_list:
        print one_taxon
        sql = 'select seqfeature_id,phylum, count(*) as n from blastnr.blastnr_%s A ' \
              ' inner join blastnr.blastnr_taxonomy B on A.subject_taxid=B.taxon_id ' \
              ' where A.query_taxon_id=%s and hit_number<=%s group by seqfeature_id,phylum order by seqfeature_id, n DESC;' % (biodb,
                                                                                                                       one_taxon,
                                                                                                                       n_hits)
        sql2 = 'select seqfeature_id,phylum from blastnr.blastnr_%s A ' \
               ' inner join blastnr.blastnr_taxonomy B on A.subject_taxid=B.taxon_id  ' \
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
                        sql = 'insert into blastnr.blastnr_majority_phylum_%s values(%s,%s,"%s", "%s")' % (biodb,
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

    import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    # n 0 hits
    # n 1-99 hits
    # n 100 hits

    sql = 'CREATE TABLE blastnr.count_n_blast_%s (taxon_id INT, total INT, n_no_hits INT, n_100_hits INT, n_less_100_hits INT)' % biodb

    server.adaptor.execute(sql)

    sql = 'select query_taxon_id, seqfeature_id, count(*) as n from blastnr.blastnr_%s group by query_taxon_id, seqfeature_id;' % biodb
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

    sql = 'select taxon_id, count(*) from orthology_detail_%s group by taxon_id' % biodb

    taxon2n_proteins = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    for taxon in taxon2n_proteins:
        total = taxon2n_proteins[taxon]
        n_100 = taxon2n_100[int(taxon)]
        n_less_than_100 = taxon2n_less_than_100[int(taxon)]
        sql = 'insert into blastnr.count_n_blast_%s values (%s,%s,%s,%s, %s)' % (biodb,
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
    import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table blastnr.BBH_taxo_hit_number_%s_%s (taxon_id INT, n_top_phylum INT, ' \
          ' n_non_top_phylum INT, n_euk INT, n_arch INT, n_virus INT, undefined INT)' % (hit_number,biodb)
    print sql
    server.adaptor.execute(sql,)
    server.commit()

    sql = 'select taxon_id, count(*) as n, best_hit_phylum from blastnr.blastnr_majority_phylum_%s' \
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
        sql_bacteria_chlamydiae = 'select t1.query_taxon_id, count(*) as n from blastnr.blastnr_%s t1 inner join blastnr.blastnr_taxonomy t2 ' \
                                  ' on t1.subject_taxid=t2.taxon_id ' \
                                  ' where t1.hit_number=%s and t1.query_taxon_id=%s and t2.superkingdom="Bacteria" ' \
                                  ' and t2.phylum="%s" group by t1.query_taxon_id;' % (biodb,
                                                                                       hit_number,
                                                                                       taxon_id,
                                                                                       taxon2most_freq_phylum[taxon_id]
                                                                                       )

        sql_bacteria_non_chlamydia = 'select t1.query_taxon_id, count(*) as n from blastnr.blastnr_%s t1 inner join blastnr.blastnr_taxonomy t2 ' \
                                  ' on t1.subject_taxid=t2.taxon_id ' \
                                  ' where t1.hit_number=%s and t1.query_taxon_id=%s and t2.superkingdom="Bacteria" ' \
                                  ' and t2.phylum!="%s" group by t1.query_taxon_id;' % (biodb,
                                                                                        hit_number,
                                                                                        taxon_id,
                                                                                        taxon2most_freq_phylum[taxon_id])

        sql_euk = 'select t1.query_taxon_id, count(*) as n from blastnr.blastnr_%s t1 inner join blastnr.blastnr_taxonomy t2 ' \
                  ' on t1.subject_taxid=t2.taxon_id where t1.hit_number=%s and t1.query_taxon_id=%s' \
                  ' and t2.superkingdom="Eukaryota" group by t1.query_taxon_id;' % (biodb,
                                                                                    hit_number,
                                                                                    taxon_id)

        sql_archae = 'select t1.query_taxon_id, count(*) as n from blastnr.blastnr_%s t1 inner join blastnr.blastnr_taxonomy t2 ' \
                  ' on t1.subject_taxid=t2.taxon_id where t1.hit_number=%s and t1.query_taxon_id=%s' \
                  ' and t2.superkingdom="Archaea" group by t1.query_taxon_id;' % (biodb,
                                                                                  hit_number,
                                                                                  taxon_id)

        sql_viruses = 'select t1.query_taxon_id, count(*) as n from blastnr.blastnr_%s t1 inner join blastnr.blastnr_taxonomy t2 ' \
                  ' on t1.subject_taxid=t2.taxon_id where t1.hit_number=%s and t1.query_taxon_id=%s' \
                  ' and t2.superkingdom="Viruses" group by t1.query_taxon_id;' % (biodb,
                                                                                  hit_number,
                                                                                  taxon_id)

        sql_undefined = 'select t1.query_taxon_id, count(*) as n from blastnr.blastnr_%s t1 inner join blastnr.blastnr_taxonomy t2 ' \
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
    import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table blastnr.BBH_taxo_%s (taxon_id INT, n_top_phylum INT, ' \
          ' n_non_top_phylum INT, n_euk INT, n_arch INT, n_virus INT, undefined INT)' % (biodb)
    print sql
    server.adaptor.execute(sql,)
    server.commit()

    sql = 'select taxon_id, count(*) as n, best_hit_phylum from blastnr.blastnr_majority_phylum_%s' \
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
        sql_phylum = 'select t1.seqfeature_id, t2.phylum, count(*) as n from blastnr.blastnr_%s t1 inner join blastnr.blastnr_taxonomy t2 ' \
                                  ' on t1.subject_taxid=t2.taxon_id ' \
                                  ' where t1.query_taxon_id=%s and t2.superkingdom="Bacteria" ' \
                                  ' and t2.phylum="%s" group by t1.seqfeature_id,t2.phylum order by t1.seqfeature_id,n DESC;' % (biodb,
                                                                                       taxon_id,
                                                                                       taxon2most_freq_phylum[taxon_id]
                                                                                       )
        print sql_phylum
        sql_superkingdom = 'select t1.seqfeature_id, t2.superkingdom, count(*) as n from blastnr.blastnr_%s t1 inner join blastnr.blastnr_taxonomy t2 ' \
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
    count_majority_phylum(args.biodb, 2)
