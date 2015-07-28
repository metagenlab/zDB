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

    print sql_head
    import sys
    #sys.exit()

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
        #sys.exit()



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
             ' where interpro_accession = "%s" group by taxon_id;;' % (db_name, accession)

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

    taxonid2genome = manipulate_biosqldb.taxon_id2genome_description(server, "chlamydia_03_15", True)

    taxons_ids = [taxonid2genome[int(i)] for i in all_taxons_id]

    f.write('"id"\t"' + '"\t"'.join(taxons_ids) + '"\n')

    for row in mat:
        row = [str(i) for i in row]
        f.write("\t".join(row) + "\n")


if __name__ == '__main__':
    import argparse
    import manipulate_biosqldb

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--database_name', type=str, help="Database name")

    args = parser.parse_args()

    #create_comparative_tables(args.database_name, "Pfam")
    #create_comparative_tables(args.database_name, "COG")
    #create_comparative_tables(args.database_name, "EC")
    #create_comparative_tables(args.database_name, "interpro")
    collect_interpro(args.database_name)
    collect_pfam(args.database_name)
    collect_interpro(args.database_name)

    #get_mysql_table("chlamydia_03_15", "Pfam")