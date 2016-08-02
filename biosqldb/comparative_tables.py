#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


def create_comparative_tables(db_name, table_name):

    # create id column + one_column per taxon_id
    import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(db_name)

    taxon_id_list = manipulate_biosqldb.get_taxon_id_list(server, db_name)

    sql = "CREATE TABLE comparative_tables.%s_%s(id VARCHAR(100) NOT NULL" % (table_name, db_name)
    for i in taxon_id_list:
        sql+=" ,`%s` INT" % i
    sql+=")"
    server.adaptor.execute(sql)

def collect_pfam(db_name):

    import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(db_name)
    taxon_id_list = manipulate_biosqldb.get_taxon_id_list(server, db_name)
    sql_head = 'INSERT INTO comparative_tables.Pfam_%s (id,' % db_name

    for taxon in taxon_id_list:
        sql_head += '`%s`,' % taxon
    sql_head = sql_head[0:-1] + ') values ('

    all_pfam_ids_sql = 'select signature_accession from interpro_%s where analysis="Pfam" ' \
                       'group by signature_accession;' % db_name

    all_pfam_ids = [i[0] for i in server.adaptor.execute_and_fetchall(all_pfam_ids_sql,)]

    i = 0
    for accession in all_pfam_ids:
        print i,'/', len(all_pfam_ids), accession
        i+=1
        sql= 'select taxon_id, count(*) from biosqldb.interpro_%s ' \
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



def collect_COGs(db_name):

    import manipulate_biosqldb
    import sys

    server, db = manipulate_biosqldb.load_db(db_name)

    taxon_id_list = manipulate_biosqldb.get_taxon_id_list(server, db_name)

    sql_head = 'INSERT INTO comparative_tables.COG_%s (id,' % db_name

    for taxon in taxon_id_list:
        sql_head += '`%s`,' % taxon
    sql_head = sql_head[0:-1] + ') values ('

    all_cogs_ids_sql = 'select COG_id from COG.locus_tag2gi_hit_%s group by COG_id;' % db_name

    all_cogs_ids = [i[0] for i in server.adaptor.execute_and_fetchall(all_cogs_ids_sql,)]

    database_id = server.adaptor.execute_and_fetchall('select biodatabase_id from biosqldb.biodatabase where name="%s"' % db_name)[0][0]

    i = 0
    for accession in all_cogs_ids:
        print i,'/', len(all_cogs_ids), accession
        i+=1
        sql= 'select B.taxon_id, A.count from (select accession,count(*) as count from ' \
             'COG.locus_tag2gi_hit_%s ' \
             'where COG_id = "%s" group by accession) A ' \
             'left join biosqldb.bioentry as B on A.accession=B.accession and biodatabase_id = %s;' % (db_name, accession, database_id)

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
        #sys.exit()


def collect_interpro(db_name):

    import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(db_name)

    taxon_id_list = manipulate_biosqldb.get_taxon_id_list(server, db_name)

    sql_head = 'INSERT INTO comparative_tables.interpro_%s (id,' % db_name

    for taxon in taxon_id_list:
        sql_head += '`%s`,' % taxon
    sql_head = sql_head[0:-1] + ') values ('

    print sql_head
    import sys
    #sys.exit()

    all_interpro_ids_sql = 'select interpro_accession from biosqldb.interpro_%s ' \
                       ' where interpro_accession != "0" group by interpro_accession;' % db_name

    all_interpro_ids = [i[0] for i in server.adaptor.execute_and_fetchall(all_interpro_ids_sql,)]

    i = 0
    for accession in all_interpro_ids:
        print i,'/', len(all_interpro_ids), accession
        i+=1
        sql= 'select taxon_id, count(*) from biosqldb.interpro_%s ' \
             ' where interpro_accession = "%s" group by taxon_id;' % (db_name, accession)

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



def collect_ko(db_name):

    import manipulate_biosqldb



    server, db = manipulate_biosqldb.load_db(db_name)

    sql_db_id = 'select biodatabase_id from biodatabase where name="%s"' % db_name

    biodb_id = server.adaptor.execute_and_fetchall(sql_db_id,)[0][0]

    taxon_id_list = manipulate_biosqldb.get_taxon_id_list(server, db_name)

    sql_head = 'INSERT INTO comparative_tables.ko_%s (id,' % db_name

    for taxon in taxon_id_list:
        sql_head += '`%s`,' % taxon
    sql_head = sql_head[0:-1] + ') values ('

    all_ko_ids_sql = 'select distinct ko_id from enzyme.locus2ko_%s;' % (db_name)

    all_ko_ids = [i[0] for i in server.adaptor.execute_and_fetchall(all_ko_ids_sql,)]

    i = 0
    for accession in all_ko_ids:
        print i,'/', len(all_ko_ids), accession
        i+=1

        sql = 'select taxon_id, count(*) from enzyme.locus2ko_%s where ko_id="%s" group by taxon_id;' % (db_name, accession)

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

def collect_EC(db_name):

    import manipulate_biosqldb



    server, db = manipulate_biosqldb.load_db(db_name)

    sql_db_id = 'select biodatabase_id from biodatabase where name="%s"' % db_name

    biodb_id = server.adaptor.execute_and_fetchall(sql_db_id,)[0][0]

    print "biodb_id", biodb_id

    taxon_id_list = manipulate_biosqldb.get_taxon_id_list(server, db_name)

    sql_head = 'INSERT INTO comparative_tables.EC_%s (id,' % db_name

    for taxon in taxon_id_list:
        sql_head += '`%s`,' % taxon
    sql_head = sql_head[0:-1] + ') values ('

    all_EC_ids_sql = 'select distinct ec from enzyme.locus2ec_%s as t1 inner join enzyme.enzymes as t2 on t1.ec_id=t2.enzyme_id;' % (db_name)

    all_ec_ids = [i[0] for i in server.adaptor.execute_and_fetchall(all_EC_ids_sql,)]

    i = 0
    for accession in all_ec_ids:
        print i,'/', len(all_ec_ids), accession
        i+=1

        sql = 'select taxon_id, count(*) from (select taxon_id, t1.accession, ec_id from enzyme.locus2ec_%s as t1 ' \
                     ' inner join biosqldb.bioentry as t2 on t1.accession=t2.accession where biodatabase_id=%s) A ' \
                     ' inner join enzyme.enzymes as B on A.ec_id=B.enzyme_id where ec="%s" group by taxon_id;' % (db_name, biodb_id, accession)

        print sql
        #sql= 'select taxon_id, count(*) from enzyme.locus2ec_%s ' \
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

def get_mysql_table(db_name, table_name):
    import numpy as np

    server, db = manipulate_biosqldb.load_db(db_name)

    all_taxons_id = manipulate_biosqldb.get_taxon_id_list(server, db_name)

    sql_taxons = "id, "
    for i in range(0, len(all_taxons_id)-1):
        sql_taxons += ' `%s`,' % all_taxons_id[i]
    sql_taxons += ' `%s`' % all_taxons_id[-1]

    sql = "select %s from comparative_tables.%s_%s" % (sql_taxons, table_name, db_name)
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

    import manipulate_biosqldb
    import sys

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = "CREATE TABLE comparative_tables.shared_orthogroups_%s(taxon_1 INT NOT NULL," \
          " taxon_2 INT NOT NULL," \
          " n_shared_orthogroups INT)" % (db_name)

    server.adaptor.execute(sql)

    orthodico = manipulate_biosqldb.get_orthogroup_count_dico(server, db_name)

    for taxon_1 in orthodico.keys():
        for taxon_2 in orthodico.keys():
            sys.stdout.write("%s\t%s\n" % (taxon_1, taxon_2))
            sql = 'insert into comparative_tables.shared_orthogroups_%s(taxon_1, taxon_2, n_shared_orthogroups) ' \
                  ' VALUES ("%s", "%s", %s)' % (db_name, taxon_1, taxon_2, orthodico[taxon_1][taxon_2])
            server.adaptor.execute(sql)
            server.adaptor.commit()


def identity_closest_homolog(db_name):

    import manipulate_biosqldb
    import biosql_own_sql_tables
    import sys

    server, db = manipulate_biosqldb.load_db(db_name)


    sql2 = "CREATE TABLE comparative_tables.identity_closest_homolog_%s(taxon_1 INT NOT NULL," \
          " taxon_2 INT NOT NULL," \
          " locus_1 VARCHAR(100) NOT NULL," \
          " locus_2 VARCHAR(100) NOT NULL," \
          " identity FLOAT)" % (db_name)



    server.adaptor.execute(sql2)

    #identitydico = biosql_own_sql_tables.calculate_average_protein_identity_new_tables(db_name)
    taxon2description = manipulate_biosqldb.taxon_id2genome_description(server, biodatabase_name=db_name)
    print len(taxon2description)
    all_taxons = taxon2description.keys()
    for i, taxon_1 in enumerate(all_taxons):

        locus2identity = biosql_own_sql_tables.circos_locus2taxon_highest_identity(db_name, taxon_1)
        for taxon_2 in all_taxons[i+1:]:
            print taxon_1, taxon_2
            for locus in locus2identity:
                try:
                    #print taxon_1, taxon_2, locus, locus2identity[locus][long(taxon_2)][1], locus2identity[locus][long(taxon_2)][0]
                    sys.stdout.write("%s\t%s\n" % (taxon_1, taxon_2))
                    sql = 'insert into comparative_tables.identity_closest_homolog_%s(taxon_1, taxon_2, locus_1, locus_2, identity) ' \
                          ' VALUES ("%s", "%s", "%s", "%s", %s)' % (db_name,
                                                                    taxon_1,
                                                                    taxon_2,
                                                                    locus,
                                                                    locus2identity[locus][long(taxon_2)][1],
                                                                    locus2identity[locus][long(taxon_2)][0])
                    server.adaptor.execute(sql)


                except:
                    #print 'no homologs'
                    continue
            server.adaptor.commit()

def shared_orthogroups_average_identity(db_name):

    import manipulate_biosqldb
    import biosql_own_sql_tables
    import sys
    import numpy

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = "CREATE TABLE comparative_tables.shared_orthogroups_average_identity_%s(taxon_1 INT NOT NULL," \
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
            data_sql = 'select identity from comparative_tables.identity_closest_homolog_%s where taxon_1=%s and taxon_2=%s' % (db_name,
                                                                                                             taxon_1,
                                                                                                             taxon_2)
            data = list([i[0] for i in server.adaptor.execute_and_fetchall(data_sql,)])

            sql = 'insert into comparative_tables.shared_orthogroups_average_identity_%s(taxon_1, taxon_2, average_identity,' \
                  ' median_identity, n_pairs) values (%s, %s, %s, %s, %s)' % (db_name,
                                                                              taxon_1,
                                                                              taxon_2,
                                                                              numpy.average(data),
                                                                              numpy.median(data),
                                                                              len(data))
            server.adaptor.execute_and_fetchall(sql,)
        server.adaptor.commit()

if __name__ == '__main__':
    import argparse
    import manipulate_biosqldb

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--database_name', type=str, help="Database name")

    args = parser.parse_args()

    #create_comparative_tables(args.database_name, "EC")
    

    '''
    create_comparative_tables(args.database_name, "Pfam")
    create_comparative_tables(args.database_name, "COG")
    create_comparative_tables(args.database_name, "EC")
    create_comparative_tables(args.database_name, "interpro")
    collect_interpro(args.database_name)
    collect_pfam(args.database_name)
    collect_interpro(args.database_name)
    collect_COGs(args.database_name)
    collect_EC(args.database_name)
    '''

    create_comparative_tables(args.database_name, "ko")
    collect_ko(args.database_name)

    #n_shared_orthogroup_table(args.database_name)
    #identity_closest_homolog(args.database_name)
    #shared_orthogroups_average_identity(args.database_name)


    #get_mysql_table("chlamydia_03_15", "Pfam")
    
