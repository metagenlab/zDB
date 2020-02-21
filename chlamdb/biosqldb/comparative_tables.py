#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


def create_comparative_tables(db_name, table_name):

    # create id column + one_column per taxon_id
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(db_name)

    taxon_id_list = manipulate_biosqldb.get_taxon_id_list(server, db_name)

    sql = "CREATE TABLE comparative_tables_%s(id VARCHAR(100) NOT NULL" % (table_name, db_name)
    for i in taxon_id_list:
        sql+=" ,`%s` INT" % i
    sql+=")"
    server.adaptor.execute(sql)

def get_all_accessions(db_name):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'select accession from bioentry t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
          ' where t2.name="%s"' % db_name
    accession_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
    return accession_list

def create_comparative_tables_accession(db_name, table_name):

    '''
    create a presence/absence matrix based on accession (and not taxon_ids)
    '''

    # create id column + one_column per taxon_id
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(db_name)

    accession_list = get_all_accessions(db_name)

    sql = "CREATE TABLE comparative_tables_%s_accessions(id VARCHAR(100) NOT NULL" % (table_name, db_name)

    for i in accession_list:
        sql+=" ,%s INT" % i

    #for i in accession_list:
    #    sql+=" ,index %s(%s)" % (i, i)
    sql+=")"
    server.adaptor.execute(sql)


def collect_pfam(db_name):

    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(db_name)
    taxon_id_list = manipulate_biosqldb.get_taxon_id_list(server, db_name)
    sql_head = 'INSERT INTO comparative_tables_Pfam (id,' % db_name

    for taxon in taxon_id_list:
        sql_head += '`%s`,' % taxon
    sql_head = sql_head[0:-1] + ') values ('

    all_pfam_ids_sql = 'select signature_accession from interpro where analysis="Pfam" ' \
                       'group by signature_accession;' % db_name

    all_pfam_ids = [i[0] for i in server.adaptor.execute_and_fetchall(all_pfam_ids_sql,)]

    i = 0
    for accession in all_pfam_ids:
        print i,'/', len(all_pfam_ids), accession
        i+=1
        sql= 'select taxon_id, count(*) from interpro ' \
             ' where analysis="Pfam" and signature_accession="%s" group by taxon_id;' % (db_name, accession)

        data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

        sql_temp = sql_head + '"%s",' % accession



        for taxon in taxon_id_list:
            try:
                sql_temp += '%s,' % data[str(taxon)]
            except:
                sql_temp += '0,'
        sql_temp = sql_temp[0:-1] + ');'

        print sql_head
        server.adaptor.execute(sql_temp,)
        server.adaptor.commit()

def collect_Pfam_accession(db_name):

    '''

    collect presence/absence of Pfam domains for each genome/plasmid accession

    :param db_name:
    :return:
    '''

    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(db_name)


    sql_head = 'INSERT INTO comparative_tables_Pfam_accessions (id,'

    accession_list = get_all_accessions(db_name)
    sql_head += ','.join(accession_list) + ') values ('

    all_pfam_ids_sql = 'select signature_accession from interpro where analysis="Pfam" ' \
                       'group by signature_accession;' % db_name

    all_pfam_ids = [i[0] for i in server.adaptor.execute_and_fetchall(all_pfam_ids_sql,)]

    sql = 'select t1.accession, t1.description from bioentry t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id' \
          ' where t2.name="%s"' % db_name
    accession2organism = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    i = 0
    for accession in all_pfam_ids:
        print i,'/', len(all_pfam_ids), accession
        i+=1
        sql= 'select organism, count(*) from interpro ' \
             ' where analysis="Pfam" and signature_accession="%s" group by organism;' % (db_name, accession)

        data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

        sql_temp = sql_head + '"%s",' % accession

        for accession in accession_list:
            try:
                sql_temp += '%s,' % data[accession2organism[str(accession)]]
            except:
                sql_temp += '0,'
        sql_temp = sql_temp[0:-1] + ');'
        print sql_temp
        server.adaptor.execute(sql_temp,)
        server.adaptor.commit()

def collect_COGs(db_name):

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(db_name)

    taxon_id_list = manipulate_biosqldb.get_taxon_id_list(server, db_name)

    sql_head = 'INSERT INTO comparative_tables_COG (id,' % db_name

    for taxon in taxon_id_list:
        sql_head += '`%s`,' % taxon
    sql_head = sql_head[0:-1] + ') values ('

    all_cogs_ids_sql = 'select COG_id from COG_locus_tag2gi_hit group by COG_id;' % db_name

    all_cogs_ids = [i[0] for i in server.adaptor.execute_and_fetchall(all_cogs_ids_sql,)]

    database_id = server.adaptor.execute_and_fetchall('select biodatabase_id from biosqldb.biodatabase where name="%s"' % db_name)[0][0]

    i = 0
    for accession in all_cogs_ids:
        print i,'/', len(all_cogs_ids), accession
        i+=1
        sql = 'select B.taxon_id, A.count from (select accession,count(*) as count from ' \
             'COG_locus_tag2gi_hit ' \
             'where COG_id = "%s" group by accession) A ' \
             'left join biosqldb.bioentry as B on A.accession=B.accession and biodatabase_id = %s;' % (db_name, accession, database_id)
        print sql
        data = server.adaptor.execute_and_fetchall(sql,)

        taxon2count = {}
        for n in data:
            if n[0] not in taxon2count:
                taxon2count[n[0]] = int(n[1])
            else:
                taxon2count[n[0]] += int(n[1])
        
        sql_temp = sql_head + '"%s",' % accession

        for taxon in taxon_id_list:
            try:
                sql_temp += '%s,' % taxon2count[int(taxon)]
            except:
                sql_temp += '0,'
        sql_temp = sql_temp[0:-1] + ');'

        print sql_temp
        server.adaptor.execute(sql_temp,)
        server.adaptor.commit()
        #sys.exit()

def collect_COGs_accession(db_name):

    from chlamdb.biosqldb import manipulate_biosqldb
    import sys

    server, db = manipulate_biosqldb.load_db(db_name)

    accession_list = get_all_accessions(db_name)

    sql_head = 'INSERT INTO comparative_tables_COG_accessions (id,' % db_name

    for accession in accession_list:
        sql_head += '%s,' % accession
    sql_head = sql_head[0:-1] + ') values ('

    all_cogs_ids_sql = 'select COG_id from COG_locus_tag2gi_hit group by COG_id;' % db_name

    all_cogs_ids = [i[0] for i in server.adaptor.execute_and_fetchall(all_cogs_ids_sql,)]

    database_id = server.adaptor.execute_and_fetchall('select biodatabase_id from biosqldb.biodatabase where name="%s"' % db_name)[0][0]

    i = 0
    for accession in all_cogs_ids:
        print i,'/', len(all_cogs_ids), accession
        i+=1
        sql= 'select accession,count(*) as c from ' \
             'COG_locus_tag2gi_hit ' \
             'where COG_id="%s" group by accession' % (db_name, accession)

        print sql
        data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

        sql_temp = sql_head + '"%s",' % accession

        for accession in accession_list:
            try:
                sql_temp += '%s,' % data[str(accession)]
            except:
                sql_temp += '0,'
        sql_temp = sql_temp[0:-1] + ');'

        print sql_temp
        server.adaptor.execute(sql_temp,)
        server.adaptor.commit()
        #sys.exit()

def collect_interpro(db_name):

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(db_name)

    taxon_id_list = manipulate_biosqldb.get_taxon_id_list(server, db_name)

    sql_head = 'INSERT INTO comparative_tables_interpro (id,' % db_name

    for taxon in taxon_id_list:
        sql_head += '`%s`,' % taxon
    sql_head = sql_head[0:-1] + ') values ('

    print sql_head
    import sys
    #sys.exit()

    all_interpro_ids_sql = 'select interpro_accession from interpro ' \
                       ' where interpro_accession != "0" group by interpro_accession;' % db_name

    all_interpro_ids = [i[0] for i in server.adaptor.execute_and_fetchall(all_interpro_ids_sql,)]

    i = 0
    for accession in all_interpro_ids:
        print i,'/', len(all_interpro_ids), accession
        i+=1
        sql= 'select A.taxon_id, count(*) as n from (select taxon_id, locus_tag, interpro_accession ' \
             ' from interpro where interpro_accession="%s" group by locus_tag) A group by taxon_id;' % (db_name, accession)

        data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

        sql_temp = sql_head + '"%s",' % accession



        for taxon in taxon_id_list:
            try:
                sql_temp += '%s,' % data[str(taxon)]
            except:
                sql_temp += '0,'
        sql_temp = sql_temp[0:-1] + ');'

        print sql_head
        server.adaptor.execute(sql_temp,)
        server.adaptor.commit()
        #sys.exit()

def collect_interpro_accession(db_name):

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(db_name)

    accession_list = get_all_accessions(db_name)

    sql_head = 'INSERT INTO comparative_tables_interpro_accessions (id,'
    
    for taxon in accession_list:
        sql_head += '%s,' % taxon
    sql_head = sql_head[0:-1] + ') values ('


    all_interpro_ids_sql = 'select interpro_accession from interpro ' \
                       ' where interpro_accession != "0" group by interpro_accession;' % db_name

    all_interpro_ids = [i[0] for i in server.adaptor.execute_and_fetchall(all_interpro_ids_sql,)]

    i = 0
    for accession in all_interpro_ids:
        print i,'/', len(all_interpro_ids), accession
        i+=1
        # group first by locus then by organism (mu.tiple occurances in a single locus_tag count as 1)
        sql= 'select A.organism, count(*) as n from (select organism, locus_tag, count(*) as n ' \
             ' from interpro where interpro_accession="%s" group by organism, locus_tag) A group by organism;' % (db_name, accession)

        data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

        sql_temp = sql_head + '"%s",' % accession

        for accession in accession_list:
            try:
                sql_temp += '%s,' % data[str(accession)]
            except:
                sql_temp += '0,'
        sql_temp = sql_temp[0:-1] + ');'

        print sql_temp
        server.adaptor.execute(sql_temp,)
        server.adaptor.commit()
        #sys.exit()

def collect_ko(db_name):

    from chlamdb.biosqldb import manipulate_biosqldb



    server, db = manipulate_biosqldb.load_db(db_name)

    sql_db_id = 'select biodatabase_id from biodatabase where name="%s"' % db_name

    biodb_id = server.adaptor.execute_and_fetchall(sql_db_id,)[0][0]

    taxon_id_list = manipulate_biosqldb.get_taxon_id_list(server, db_name)

    sql_head = 'INSERT INTO comparative_tables_ko (id,' % db_name

    for taxon in taxon_id_list:
        sql_head += '`%s`,' % taxon
    sql_head = sql_head[0:-1] + ') values ('

    all_ko_ids_sql = 'select distinct ko_id from enzyme_locus2ko;' % (db_name)

    all_ko_ids = [i[0] for i in server.adaptor.execute_and_fetchall(all_ko_ids_sql,)]

    i = 0
    for accession in all_ko_ids:
        print i,'/', len(all_ko_ids), accession
        i+=1

        sql = 'select taxon_id, count(*) from enzyme_locus2ko where ko_id="%s" group by taxon_id;' % (db_name, accession)

        print sql

        data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
        sql_temp = sql_head + '"%s",' % accession

        for taxon in taxon_id_list:
            try:
                sql_temp += '%s,' % data[str(taxon)]
            except:
                sql_temp += '0,'
        sql_temp = sql_temp[0:-1] + ');'

        print sql_temp
        server.adaptor.execute(sql_temp,)
        server.adaptor.commit()

def collect_ko_accession(db_name):

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(db_name)

    sql_db_id = 'select biodatabase_id from biodatabase where name="%s"' % db_name

    biodb_id = server.adaptor.execute_and_fetchall(sql_db_id,)[0][0]

    accession_list = get_all_accessions(db_name)

    sql_head = 'INSERT INTO comparative_tables_ko_accessions (id,'

    for accession in accession_list:
        sql_head += '%s,' % accession
    sql_head = sql_head[0:-1] + ') values ('

    all_ko_ids_sql = 'select distinct ko_id from enzyme_locus2ko;' % (db_name)

    all_ko_ids = [i[0] for i in server.adaptor.execute_and_fetchall(all_ko_ids_sql,)]

    i = 0
    for accession in all_ko_ids:
        print i,'/', len(all_ko_ids), accession
        i+=1

        sql = 'select B.accession, count(*) as n from (select locus_tag, ko_id from enzyme_locus2ko ' \
              ' where ko_id="%s") A inner join orthology_detail as B on A.locus_tag=B.locus_tag ' \
              ' group by accession;' % (db_name,
                                        accession,
                                        db_name)

        print sql

        data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
        sql_temp = sql_head + '"%s",' % accession

        for accession in accession_list:
            try:
                sql_temp += '%s,' % data[str(accession)]
            except:
                sql_temp += '0,'
        sql_temp = sql_temp[0:-1] + ');'

        print sql_temp
        server.adaptor.execute(sql_temp,)
        server.adaptor.commit()


def collect_EC(db_name):

    from chlamdb.biosqldb import manipulate_biosqldb



    server, db = manipulate_biosqldb.load_db(db_name)

    sql_db_id = 'select biodatabase_id from biodatabase where name="%s"' % db_name

    biodb_id = server.adaptor.execute_and_fetchall(sql_db_id,)[0][0]

    print "biodb_id", biodb_id

    taxon_id_list = manipulate_biosqldb.get_taxon_id_list(server, db_name)

    sql_head = 'INSERT INTO comparative_tables_EC (id,' % db_name

    for taxon in taxon_id_list:
        sql_head += '`%s`,' % taxon
    sql_head = sql_head[0:-1] + ') values ('

    all_EC_ids_sql = 'select distinct ec from enzyme_locus2ec as t1 inner join enzyme_enzymes as t2 on t1.ec_id=t2.enzyme_id;' % (db_name)

    all_ec_ids = [i[0] for i in server.adaptor.execute_and_fetchall(all_EC_ids_sql,)]

    i = 0
    for accession in all_ec_ids:
        print i,'/', len(all_ec_ids), accession
        i+=1

        sql = 'select taxon_id, count(*) from (select taxon_id, t1.accession, ec_id from enzyme_locus2ec as t1 ' \
                     ' inner join biosqldb.bioentry as t2 on t1.accession=t2.accession where biodatabase_id=%s) A ' \
                     ' inner join enzyme_enzymes as B on A.ec_id=B.enzyme_id where ec="%s" group by taxon_id;' % (db_name, biodb_id, accession)

        print sql
        #sql= 'select taxon_id, count(*) from enzyme_locus2ec ' \
        #     ' where ec = "%s" group by taxon_id;' % (db_name, accession)

        data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
        sql_temp = sql_head + '"%s",' % accession

        for taxon in taxon_id_list:
            try:
                sql_temp += '%s,' % data[str(taxon)]
            except:
                sql_temp += '0,'
        sql_temp = sql_temp[0:-1] + ');'

        print sql_temp
        server.adaptor.execute(sql_temp,)
        server.adaptor.commit()

def collect_EC_accession(db_name):

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(db_name)

    accession_list = get_all_accessions(db_name)

    sql_head = 'INSERT INTO comparative_tables_EC_accessions (id,'

    for accession in accession_list:
        sql_head += '%s,' % accession
    sql_head = sql_head[0:-1] + ') values ('

    all_EC_ids_sql = 'select distinct ec from enzyme_locus2ec as t1 inner join enzyme_enzymes as t2 on t1.ec_id=t2.enzyme_id;' % (db_name)

    all_ec_ids = [i[0] for i in server.adaptor.execute_and_fetchall(all_EC_ids_sql,)]

    i = 0
    for accession in all_ec_ids:
        print i,'/', len(all_ec_ids), accession
        i+=1

        sql = 'select t1.accession, count(*) as n from enzyme_locus2ec as t1 inner join enzyme_enzymes as t2 ' \
              ' on t1.ec_id=t2.enzyme_id where ec="%s" group by accession;' % (db_name, accession)

        print sql

        data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
        sql_temp = sql_head + '"%s",' % accession

        for accession in accession_list:
            try:
                sql_temp += '%s,' % data[str(accession)]
            except:
                sql_temp += '0,'
        sql_temp = sql_temp[0:-1] + ');'

        print sql_temp
        server.adaptor.execute(sql_temp,)
        server.adaptor.commit()

def collect_orthogroup_accession(db_name):

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(db_name)

    accession_list = get_all_accessions(db_name)

    sql_head = 'INSERT INTO comparative_tables_orthology_accessions (id,'

    for accession in accession_list:
        sql_head += '%s,' % accession
    sql_head = sql_head[0:-1] + ') values ('

    all_orthogroup_ids_sql = 'select distinct orthogroup from orthology_detail;'
    
    all_orthogroup_ids = [i[0] for i in server.adaptor.execute_and_fetchall(all_orthogroup_ids_sql,)]

    i = 0
    for accession in all_orthogroup_ids:
        print i,'/', len(all_orthogroup_ids), accession
        i+=1

        sql = 'select accession, count(*) as n from orthology_detail ' \
              ' where orthogroup="%s" group by accession;' % (db_name, accession)

        print sql

        data = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
        sql_temp = sql_head + '"%s",' % accession

        for accession in accession_list:
            try:
                sql_temp += '%s,' % data[str(accession)]
            except:
                sql_temp += '0,'
        sql_temp = sql_temp[0:-1] + ');'

        print sql_temp
        server.adaptor.execute(sql_temp,)
        server.adaptor.commit()




def get_mysql_table(db_name, table_name):
    import numpy as np

    server, db = manipulate_biosqldb.load_db(db_name)

    all_taxons_id = manipulate_biosqldb.get_taxon_id_list(server, db_name)

    sql_taxons = "id, "
    for i in range(0, len(all_taxons_id)-1):
        sql_taxons += ' `%s`,' % all_taxons_id[i]
    sql_taxons += ' `%s`' % all_taxons_id[-1]

    sql = "select %s from comparative_tables_%s" % (sql_taxons, table_name)
    #print sql
    mat = np.array(server.adaptor.execute_and_fetchall(sql,))
    #np.chararray
    f = open("%s_matrix.tab" % table_name, "w")

    taxonid2genome = manipulate_biosqldb.taxon_id2genome_description(server, db_name, True)

    taxons_ids = [taxonid2genome[int(i)] for i in all_taxons_id]

    f.write('"id"\t"' + '"\t"'.join(taxons_ids) + '"\n')

    for row in mat:
        row = [str(i) for i in row]
        f.write("\t".join(row) + "\n")


def n_shared_orthogroup_table(db_name):

    from chlamdb.biosqldb import manipulate_biosqldb
    import sys

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = "CREATE TABLE comparative_tables_shared_orthogroups(taxon_1 INT NOT NULL," \
          " taxon_2 INT NOT NULL," \
          " n_shared_orthogroups INT)"

    server.adaptor.execute(sql)

    orthodico = manipulate_biosqldb.get_orthogroup_count_dico(server, db_name)

    for taxon_1 in orthodico.keys():
        for taxon_2 in orthodico.keys():
            sys.stdout.write("%s\t%s\n" % (taxon_1, taxon_2))
            sql = 'insert into comparative_tables_shared_orthogroups(taxon_1, taxon_2, n_shared_orthogroups) ' \
                  ' VALUES ("%s", "%s", %s)' % (taxon_1, taxon_2, orthodico[taxon_1][taxon_2])
            server.adaptor.execute(sql)
            server.adaptor.commit()


def identity_closest_homolog(db_name):

    from chlamdb.biosqldb import manipulate_biosqldb
    import biosql_own_sql_tables
    import sys

    server, db = manipulate_biosqldb.load_db(db_name)

    sql1 = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % db_name
    locus2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql1,))


    sql2 = "CREATE TABLE comparative_tables_identity_closest_homolog2(taxon_1 INT NOT NULL," \
          " taxon_2 INT NOT NULL," \
          " locus_1 INT NOT NULL," \
          " locus_2 INT NOT NULL," \
          " identity FLOAT, index locus_1(locus_1)," \
          " index locus_2(locus_2), index taxon_1(taxon_1), index taxon_2(taxon_2))"

    server.adaptor.execute(sql2)

    #identitydico = biosql_own_sql_tables.calculate_average_protein_identity_new_tables(db_name)
    taxon2description = manipulate_biosqldb.taxon_id2genome_description(server, biodatabase_name=db_name)

    all_taxons = taxon2description.keys()
    for i, taxon_1 in enumerate(all_taxons):

        locus2identity = biosql_own_sql_tables.circos_locus2taxon_highest_identity(db_name, taxon_1)
        for taxon_2 in all_taxons:
            if taxon_1 == taxon_2:
                continue

            for locus in locus2identity:
                try:
                    #print taxon_1, taxon_2, locus, locus2identity[locus][long(taxon_2)][1], locus2identity[locus][long(taxon_2)][0]
                    sys.stdout.write("%s\t%s\n" % (taxon_1, taxon_2))
                    sql = 'insert into comparative_tables_identity_closest_homolog2(taxon_1, taxon_2, locus_1, locus_2, identity) ' \
                          ' VALUES ("%s", "%s", "%s", "%s", %s)' % (db_name,
                                                                    taxon_1,
                                                                    taxon_2,
                                                                    locus2seqfeature_id[locus],
                                                                    locus2seqfeature_id[locus2identity[locus][long(taxon_2)][1]],
                                                                    locus2identity[locus][long(taxon_2)][0])
                    server.adaptor.execute(sql)


                except:
                    #print 'no homologs'
                    continue
            server.adaptor.commit()

def shared_orthogroups_average_identity(db_name):

    from chlamdb.biosqldb import manipulate_biosqldb
    import biosql_own_sql_tables
    import sys
    import numpy

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = "CREATE TABLE comparative_tables_shared_og_av_id(taxon_1 INT NOT NULL," \
          " taxon_2 INT NOT NULL," \
          " average_identity FLOAT," \
          " median_identity FLOAT," \
          " n_pairs INT)" % (db_name)
    server.adaptor.execute(sql)


    taxon2description = manipulate_biosqldb.taxon_id2genome_description(server, biodatabase_name=db_name)
    print len(taxon2description)
    all_taxons = taxon2description.keys()
    for i, taxon_1 in enumerate(all_taxons):
        for taxon_2 in all_taxons[i+1:]:
            data_sql = 'select identity from comparative_tables_identity_closest_homolog2 where taxon_1=%s and taxon_2=%s' % (taxon_1,
                                                                                                                              taxon_2)
            data = list([i[0] for i in server.adaptor.execute_and_fetchall(data_sql,)])

            sql = 'insert into comparative_tables_shared_og_av_id(taxon_1, taxon_2, average_identity,' \
                  ' median_identity, n_pairs) values (%s, %s, %s, %s, %s)' % (taxon_1,
                                                                              taxon_2,
                                                                              numpy.average(data),
                                                                              numpy.median(data),
                                                                              len(data))
            server.adaptor.execute_and_fetchall(sql,)
        server.adaptor.commit()

if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--database_name', type=str, help="Database name")
    parser.add_argument("-o", '--orthology', help="orthology tables (n shared ortho, identity closest, average ID)", action="store_true")
    parser.add_argument("-c", '--cog', help="cog table", action="store_true")
    parser.add_argument("-p", '--pfam', help="pfam table", action="store_true")
    parser.add_argument("-e", '--ec', help="priam EC table", action="store_true")
    parser.add_argument("-i", '--interpro', help="interpro table", action="store_true")
    parser.add_argument("-k", '--ko', help="KEGG ko table", action="store_true")
    parser.add_argument("-a", '--accessions', help="accession tables", action="store_true")
    
    args = parser.parse_args()

    if args.orthology:
        identity_closest_homolog(args.database_name)
        n_shared_orthogroup_table(args.database_name)
        
        shared_orthogroups_average_identity(args.database_name)

        if args.accessions:
            create_comparative_tables_accession(args.database_name, "orthology")
            collect_orthogroup_accession(args.database_name)

    if args.cog:
        create_comparative_tables(args.database_name, "COG")
        collect_COGs(args.database_name)
        
        if args.accessions:
            create_comparative_tables_accession(args.database_name, "COG")
            collect_COGs_accession(args.database_name)
        
        
    if args.pfam:
        create_comparative_tables(args.database_name, "Pfam")
        collect_pfam(args.database_name)
        if args.accessions:
            create_comparative_tables_accession(args.database_name, 'Pfam')
            collect_Pfam_accession(args.database_name)

    if args.ec:
        
        create_comparative_tables(args.database_name, "EC")
        collect_EC(args.database_name)
        if args.accessions:
            create_comparative_tables_accession(args.database_name, "EC")
            collect_EC_accession(args.database_name)
        
    if args.interpro:
        create_comparative_tables(args.database_name, "interpro")
        collect_interpro(args.database_name)
        if args.accessions:
            create_comparative_tables_accession(args.database_name, "interpro")
            collect_interpro_accession(args.database_name)
        
        
    if args.ko:
        create_comparative_tables(args.database_name, "ko")
        collect_ko(args.database_name)
        if args.accessions:
            create_comparative_tables_accession(args.database_name, "ko")
            collect_ko_accession(args.database_name)

    #get_mysql_table("chlamydia_03_15", "Pfam")
