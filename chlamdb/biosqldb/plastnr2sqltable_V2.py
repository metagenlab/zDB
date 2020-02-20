#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# Load blast results into biosql databas
# Tablulated blast files: qgi qacc sgi sacc sscinames sskingdoms staxids evalue nident pident positive gaps length qstart qend qcovs sstart send sstrand stitle
# Genbank files of query proteins must be leaded in the biosql database already (used to fetch locus tag, taxon_id, organism name, ...)
# Create 3 new tables (if they do not exit): - blastnr_taxonomy with txaon_id and corresponding full taxonomic path
#                                            todo: use taxonomy implemented in biosqldb shema instead: currently fetch taxonomic path from ncbi and add it to blast_nr_taxonomy
#                                           - blast_nr_`biodb_name`: contain data foreach hit
#                                           - blast_nr_taxonomy_`biodb_name`: contain taxonomical data oreach hit (as one hit can now include multiple taxons/MULTISPECIES hits)
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: april 2015
# ---------------------------------------------------------------------------

'''

    #### general data ####

    hit_number
    taxon_id
    locus_tag
    organism
    orthogroup

    #### NCBI taxonomical ranks ####

    query_gi
    subject_gi
    subject_scientific_name
    subject_taxid
    evalue
    n_identical
    percent_identity
    positive
    gaps
    length
    query_start
    query_end
    query_cov
    subject_start
    subject_end]
    subject_strand
    subject_title

    no_rank
    superkingdom
    kingdom
    subkingdom
    superphylum
    phylum
    subphylum
    superclass
    class
    subclass
    superorder
    order
    suborder
    superfamily
    family
    subfamily
    genus
    subgenus
    species
    species_subgroup
    species_group
    subspecies

    tribe
    infraorder
    subtribe
    forma
    infraclass
    varietas
    parvorder
'''

def _chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

def load_blastnr_file_into_db(locus_tag2taxon_id,
                                locus_tag2seqfeature_id,
                                protein_id2seqfeature_id,
                                locus_tag2bioentry_id,
                                mysql_host,
                                mysql_user,
                                mysql_pwd,
                                mysql_db,
                                input_blast_files,
                                biodb,
                                accession2descriptions):

    import accession2taxon_id

    '''
    Load tabulated blast results into sql table blastnr_`db_name`
    Ab unique identifier (primary sql key) is attributed to each blast hsp

    :param seqfeature_id2locus_tag: dictionnary with whole data for `biodb` biodatabase
    :param locus_tag2seqfeature_id: dictionnary with whole data for `biodb` biodatabase
    :param protein_id2seqfeature_id: dictionnary with whole data for `biodb` biodatabase
    :param db_name: name of the biodatabase
    :param input_blast_files: all input tabulated blast files
    :return: None
    '''

    import time
    import MySQLdb

    n_file = 0
    for one_blast_file in input_blast_files:
        n_file +=1
        conn = MySQLdb.connect(host=mysql_host, # your host, usually localhost
                                    user=mysql_user, # your username
                                    passwd=mysql_pwd, # your password
                                    db=mysql_db) # name of the data base
        cursor = conn.cursor()


        with open(one_blast_file, 'r') as f:
            print 'Loading', n_file, one_blast_file, '...'
            input_file = [i.rstrip().split('\t') for i in f]



            print 'loading blast results into database...'
            for n, line in enumerate(input_file):
                #print n, line[0]
                # qgi qacc sgi sacc sscinames sskingdoms staxids evalue nident pident positive gaps length qstart qend qcovs sstart send sstrand stitle
                if n%1000 == 0:
                    print time.ctime() + ': ' + str(round((float(n)/len(input_file))*100,2)) + '%...'

                query_accession = line[0]

                try:
                    seqfeature_id = protein_id2seqfeature_id[query_accession]
                except KeyError:
                    seqfeature_id = locus_tag2seqfeature_id[query_accession]

                evalue = line[10]
                #n_identical = int(line[8]) # remove
                percent_identity = float(line[2])
                #positive = int(line[10]) # remove
                gaps = int(line[5])
                length = int(line[3])
                query_start = int(line[6])
                query_end = int(line[7])
                #query_cov = int(line[]) # calculate
                subject_start = int(line[8])
                subject_end = int(line[9])
                bit_score = float(line[11])
                #subject_strand = line[18] # remove

                subject_accession = line[1].split("|")[1]
                subject_gi = accession2descriptions[subject_accession]['id']#int(line[1].split("|")[1]) # remove
                subject_taxon_id = accession2descriptions[subject_accession]['taxon_id']

                #subject_gi = int(line[1].split("|")[1]) # parse it from plastp result

                try:
                    subject_kingdom = accession2descriptions[subject_accession]['taxonomy']  # get from accession2annotations
                except:
                    subject_kingdom = '-'
                subject_scientific_names = accession2descriptions[subject_accession]['source'] # get from accession2annotations
                subject_title = accession2descriptions[subject_accession]['description'] # get from accession2annotations

                # first hit of the file
                if n == 0:
                    hit_n = 1
                # check if new query accession)
                elif line[0] != input_file[n-1][0]:
                    # new query => insert its first hit into blastnrdb
                    hit_n = 1
                # same query, new hit
                else:
                    hit_n+=1

                '''
                1 query_taxon_id int
                2 query_bioentry_id INT
                3 seqfeature_id INT
                4 hit_number int
                5 subject_gi int
                6 subject_accession varchar(200)
                7 subject_kingdom varchar(200)
                8 subject_scientific_name TEXT(2000000)
                9 subject_taxid INT
                10 subject_title VARCHAR(2000)
                11 evalue varchar(200)
                12 bit_score float
                12 percent_identity float
                14 gaps int
                15 length int
                16 query_start int
                17 query_end int
                18 subject_start int
                19 subject_end
                '''

                values = '(%s,%s,%s,%s,%s,"%s","%s","%s",%s,"%s","%s",%s,%s,%s,%s,%s,%s,%s,%s);' % (locus_tag2taxon_id[query_accession],
                                                                                                    locus_tag2bioentry_id[query_accession],
                                                                                                    seqfeature_id,
                                                                                                    hit_n,
                                                                                                    subject_gi,
                                                                                                    subject_accession,
                                                                                                    subject_kingdom,
                                                                                                    subject_scientific_names,
                                                                                                    subject_taxon_id,
                                                                                                    subject_title,
                                                                                                    evalue,
                                                                                                    bit_score,
                                                                                                    percent_identity,
                                                                                                    gaps,
                                                                                                    length,
                                                                                                    query_start,
                                                                                                    query_end,
                                                                                                    subject_start,
                                                                                                    subject_end
                                                                                                    )

                sql = 'insert into blastnr_%s values %s' % (biodb,values)
                #print sql
                cursor.execute(sql,)
                conn.commit()


def create_sql_plastnr_tables(db_name, mysql_host, mysql_user, mysql_pwd, mysql_db='blastnr'):
    import MySQLdb

    conn = MySQLdb.connect(host=mysql_host, # your host, usually localhost
                                user=mysql_user, # your username
                                passwd=mysql_pwd, # your password
                                db=mysql_db) # name of the data base
    cursor = conn.cursor()



    sql_plast = 'CREATE TABLE IF NOT EXISTS blastnr_%s (query_taxon_id INT,' \
                ' query_bioentry_id INT,' \
                ' seqfeature_id INT,' \
                ' hit_number int,' \
                ' subject_gi int,' \
                ' subject_accession varchar(200),' \
                ' subject_kingdom varchar(200),' \
                ' subject_scientific_name TEXT(2000000), ' \
                ' subject_taxid INT,' \
                ' subject_title VARCHAR(2000),' \
                ' evalue varchar(200),' \
                ' bit_score float,' \
                ' percent_identity float,' \
                ' gaps int,' \
                ' length int,' \
                ' query_start int,' \
                ' query_end int,' \
                ' subject_start int,' \
                ' subject_end int,' \
                ' INDEX query_taxon_id (query_taxon_id),' \
                ' INDEX query_bioentry_id (query_bioentry_id),' \
                ' INDEX hit_number (hit_number),' \
                ' INDEX seqfeature_id (seqfeature_id),' \
                ' INDEX subject_taxid(subject_taxid))' % (db_name)

    try:

        cursor.execute(sql_plast)
        print 'sql hits ok'
        conn.commit()
    except:
        print sql_plast
        print 'not created'


    sql_taxonomy1 = 'CREATE TABLE IF NOT EXISTS blastnr_taxonomy (taxon_id int unique PRIMARY KEY, ' \
                    ' no_rank VARCHAR(200) default "-", ' \
                    ' superkingdom VARCHAR(200) default "-", ' \
                    ' kingdom VARCHAR(200) default "-", ' \
                    ' subkingdom VARCHAR(200) default "-", ' \
                    ' superphylum VARCHAR(200) default "-", ' \
                    ' phylum VARCHAR(200) default "-", ' \
                    ' subphylum VARCHAR(200) default "-", ' \
                    ' superclass VARCHAR(200) default "-", ' \
                    ' class VARCHAR(200) default "-", ' \
                    ' subclass VARCHAR(200) default "-", ' \
                    ' superorder VARCHAR(200) default "-", ' \
                    ' `order` VARCHAR(200) default "-", ' \
                    ' suborder VARCHAR(200) default "-", ' \
                    ' superfamily VARCHAR(200) default "-", ' \
                    ' family VARCHAR(200) default "-", ' \
                    ' subfamily VARCHAR(200) default "-", ' \
                    ' genus VARCHAR(200) default "-", ' \
                    ' subgenus VARCHAR(200) default "-", ' \
                    ' species VARCHAR(200) default "-", ' \
                    ' species_subgroup VARCHAR(200) default "-", ' \
                    ' species_group VARCHAR(200) default "-", ' \
                    ' subspecies VARCHAR(200) default "-", ' \
                    ' tribe VARCHAR(200) default "-", ' \
                    ' infraorder VARCHAR(200) default "-", ' \
                    ' subtribe VARCHAR(200) default "-", ' \
                    ' forma VARCHAR(200) default "-", ' \
                    ' infraclass VARCHAR(200) default "-", ' \
                    ' varietas VARCHAR(200) default "-", ' \
                    ' parvorder VARCHAR(200) default "-")'
    try:

        cursor.execute(sql_taxonomy1)
        print 'sql hits ok'
        conn.commit()
    except:
        print sql_taxonomy1
        print 'not created'


def update_blastnr_taxonomy_table(blast_table_name,
                                  mysql_host="localhost",
                                  mysql_user="root",
                                  mysql_pwd="baba",
                                  mysql_db="blastnr"):
    import MySQLdb
    from plastnr2sqltable import insert_taxons_into_sqldb
    conn = MySQLdb.connect(host=mysql_host, # your host, usually localhost
                                user=mysql_user, # your username
                                passwd=mysql_pwd, # your password
                                db=mysql_db) # name of the data base
    cursor = conn.cursor()

    sql_taxonomy1 = 'CREATE TABLE IF NOT EXISTS blastnr_taxonomy (taxon_id int unique PRIMARY KEY, ' \
                    ' no_rank VARCHAR(200) default "-", ' \
                    ' superkingdom VARCHAR(200) default "-", ' \
                    ' kingdom VARCHAR(200) default "-", ' \
                    ' subkingdom VARCHAR(200) default "-", ' \
                    ' superphylum VARCHAR(200) default "-", ' \
                    ' phylum VARCHAR(200) default "-", ' \
                    ' subphylum VARCHAR(200) default "-", ' \
                    ' superclass VARCHAR(200) default "-", ' \
                    ' class VARCHAR(200) default "-", ' \
                    ' subclass VARCHAR(200) default "-", ' \
                    ' superorder VARCHAR(200) default "-", ' \
                    ' `order` VARCHAR(200) default "-", ' \
                    ' suborder VARCHAR(200) default "-", ' \
                    ' superfamily VARCHAR(200) default "-", ' \
                    ' family VARCHAR(200) default "-", ' \
                    ' subfamily VARCHAR(200) default "-", ' \
                    ' genus VARCHAR(200) default "-", ' \
                    ' subgenus VARCHAR(200) default "-", ' \
                    ' species VARCHAR(200) default "-", ' \
                    ' species_subgroup VARCHAR(200) default "-", ' \
                    ' species_group VARCHAR(200) default "-", ' \
                    ' subspecies VARCHAR(200) default "-", ' \
                    ' tribe VARCHAR(200) default "-", ' \
                    ' infraorder VARCHAR(200) default "-", ' \
                    ' subtribe VARCHAR(200) default "-", ' \
                    ' forma VARCHAR(200) default "-", ' \
                    ' infraclass VARCHAR(200) default "-", ' \
                    ' varietas VARCHAR(200) default "-", ' \
                    ' parvorder VARCHAR(200) default "-")'
    try:

        cursor.execute(sql_taxonomy1)
        print 'sql hits ok'
        conn.commit()
    except:
        print sql_taxonomy1
        print 'not created'


    sql = 'select A.subject_taxid from (select distinct subject_taxid from %s) A ' \
          ' left join blastnr_taxonomy as B on A.subject_taxid=B.taxon_id where B.taxon_id is NULL;' % blast_table_name

    cursor.execute(sql,)
    taxon_list = [str(i[0]) for i in cursor.fetchall()]
    if len(taxon_list) > 0:
        print 'Fetching ncbi taxonomy... for %s taxons' % str(len(taxon_list))
        print taxon_list[0:10]

        # subdivide the taxon list in smaller lists, otherwise NCBI will limit the results to? 10000 (observed once only)
        insert_taxons_into_sqldb(taxon_list, 300, mysql_pwd=mysql_pwd)

def blastnr2biosql( locus_tag2seqfeature_id,
                    protein_id2seqfeature_id,
                    db_name,
                    n_procs,
                    mysql_host,
                    mysql_user,
                    mysql_pwd,
                    mysql_db,
                    input_blast_files,
                    load_blastnr_tables=True):

    import numpy
    import re
    import accession2taxon_id
    import time
    import re
    from datetime import datetime
    from multiprocessing import Process
    from plastnr2sqltable import update2biosql_blastnr_table

    print 'host1', mysql_host
    print 'user1', mysql_user

    create_sql_plastnr_tables(db_name,mysql_host,mysql_user,mysql_pwd,mysql_db)

    print "add eventual new taxons to the main blastnr taxonomy table"
    #update2biosql_blastnr_table.update2biosql_blastnr_table(mysql_host, mysql_user, mysql_pwd, mysql_db, *input_blast_files)

    # load blast data into blastnr_ db_name table
    n_cpu = n_procs
    n_poc_per_list = int(numpy.ceil(len(input_blast_files)/float(n_cpu)))
    query_lists = _chunks(input_blast_files, n_poc_per_list)
    print 'n lists: %s' % len(query_lists)



    sql = 'create table if not exists blastnr.accession2taxon_and_description (gi INT AUTO_INCREMENT PRIMARY KEY, accession varchar(600),taxon_id INT, description TEXT,' \
          ' taxonomy TEXT, source TEXT);'
    server.adaptor.execute(sql,)

    print 'get locus2taxon_id'
    sql = 'select locus_tag, taxon_id from orthology_detail' % db_name
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    print 'get locus2bioentry'
    sql2 = 'select locus_tag,bioentry_id from biodatabase t1 ' \
           ' inner join bioentry as t2 on t1.biodatabase_id=t2.biodatabase_id' \
           ' inner join orthology_detail t3 on t2.accession=t3.accession where t1.name="%s"' % (db_name,db_name)


    locus_tag2bioentry_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

    sql = 'select accession from blastnr.accession2taxon_and_description'
    accession_in_db = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    protein_accession_list = []
    for one_blast_file in input_blast_files:
                print 'blast file %s' % one_blast_file
                with open(one_blast_file, 'r') as f:
                    all_accession_ids = [i.rstrip().split("\t")[1].split('|')[1] for i in f]
                    protein_accession_list+=all_accession_ids

    print 'filtering gi from db'
    nr_protein_accession_list = list(set(protein_accession_list) - set(accession_in_db))
    #nr_protein_gi_list = list(set(protein_gi_list))
    '''
    sql = 'select accession, taxon_id from blastnr.accession2taxon_and_description'

    accession2taxon_id_dico = {}
    data = server.adaptor.execute_and_fetchall(sql,)
    for row in data:
        accession2taxon_id_dico[int(row[0])] = str(row[1])
    '''

    print 'getting protein 2 taxon id for %s proteins (out of %s)' % (len(nr_protein_accession_list), len(set(protein_accession_list)))

    id_lists = _chunks(nr_protein_accession_list, 50)

    for i, one_list in enumerate(id_lists):
        print i, " lists /", len(id_lists), str(datetime.now())
        #print one_list[0:10]
        if i % 100 == 0 and i != 0:
            time.sleep(60)
        try:
            #print 'getting taxon dico'
            accession2taxon_id_dico = accession2taxon_id.gi2taxon_id(one_list, "protein")
            #print len(accession2taxon_id_dico), accession2taxon_id_dico
            #print 'getting accession dico'
            accession2descriptions = accession2taxon_id.gi2description(one_list, "protein")
            #print len(accession2descriptions), accession2descriptions
        except:
            print 'problem getting accession2data...'
        for accession in one_list:
            try:
                descr = accession2descriptions[str(accession)]
                try:
                    descr['description'] = re.sub("'", "''", descr['description'])
                except:
                    print descr
                    import sys
                    sys.exit()
                try:
                    taxonomy = descr['taxonomy'][0]
                except:
                    taxonomy = '-'

                taxon_id = accession2taxon_id_dico[accession]
                description = re.sub('"','',descr['description'])
                source = re.sub('"','',descr['source'])

                sql = 'insert into blastnr.accession2taxon_and_description (accession, taxon_id, description, taxonomy, source) ' \
                      ' values("%s",%s,"%s", "%s", "%s")' % (accession,
                                                           taxon_id,
                                                           description,
                                                           taxonomy,
                                                           source)
                #try:
                server.adaptor.execute(sql, )
                #except:
                #    print sql
                #    import sys
                #    sys.exit()
            except:
                print sql
                print 'problem with:', accession
        server.adaptor.commit()

    if load_blastnr_tables:
        # collect annotation for all proteins in the database
        print 'collecting annotation for blastnr_table'
        sql = 'select * from blastnr.accession2taxon_and_description'
        data = server.adaptor.execute_and_fetchall(sql,)
        accession2descriptions = {}
        for row in data:
            descr = {}
            descr['description'] = row[3]
            descr['source'] = row[5]
            descr['taxonomy'] = row[4]
            descr['taxon_id'] = row[2]
            descr['id'] = row[0]
            accession2descriptions[str(row[1])] = descr
            #gi2taxon_id[row[0]] = row[2]



        if len(input_blast_files)>n_poc_per_list:

            procs = []
            for one_list in query_lists:
                proc = Process(target=load_blastnr_file_into_db, args=(locus_tag2taxon_id,
                                                                       locus_tag2seqfeature_id,
                                                                       protein_id2seqfeature_id,
                                                                       locus_tag2bioentry_id,
                                                                       mysql_host,
                                                                       mysql_user,
                                                                       mysql_pwd,
                                                                       mysql_db,
                                                                       one_list,
                                                                       db_name,
                                                                       accession2descriptions))
                procs.append(proc)
                proc.start()

            # Wait for all worker processes to finish
            for proc in procs:
                proc.join()
        else:
            print 'Only %s input blast file(s), not working in paralell mode' % (len(input_blast_files))
            load_blastnr_file_into_db(locus_tag2taxon_id,
                                       locus_tag2seqfeature_id,
                                       protein_id2seqfeature_id,
                                       locus_tag2bioentry_id,
                                       mysql_host,
                                       mysql_user,
                                       mysql_pwd,
                                       mysql_db,
                                       input_blast_files,
                                       db_name,
                                       gi2descriptions,
                                       gi2taxon_id)

    #sys.stdout.write("done!")





if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    import sys
    import os
    import json

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_blast', type=str, help="input blast tab files", nargs='+')
    parser.add_argument("-d", '--mysql_database', type=str, help="Biosql biodatabase name")
    parser.add_argument("-t", '--create_tables', action='store_true', help="Create SQL tables")
    parser.add_argument("-p", '--n_procs', type=int, help="Number of threads to use (default=8)", default=4)
    parser.add_argument("-c", '--clean_tables', action='store_true', help="delete all sql tables")
    parser.add_argument("-l", '--load_tables', action='store_true', help="load tab files into biodatabase")
    parser.add_argument("-f", '--filter_n_hits', action='store_true', help="filter_n_hits (max 100 hits/locus)")
    parser.add_argument("-u", '--update_taxo_table', action='store_true', help="update blastnr taxonomy table based on all hits taxons")
    #create_protein_accession2taxon_id_table()
    #import sys
    #sys.exit()


    args = parser.parse_args()


    mysql_host = 'localhost'
    mysql_user = 'root'

    mysql_pwd = os.environ['SQLPSW']
    mysql_db = 'blastnr'

    biodb = args.mysql_database

    if args.clean_tables:
        del_blastnr_table_content(args.mysql_database)

    '''

    Structure of the data in the biosql database: 3 tables

    - one main taxonomical table /blastnr_taxonomy/ with
            taxon id: taxonomical data

    - for each biodatabase:
        - one table storing the blast data (blastnr_ db_name)
        - one table storing taxon ids of each hsps (blastnr_taxonomy_ db_name)
        => TODO this last table is redundant as each hsp of a single hit should have the same taxon ids??

    '''
    if args.load_tables:

        server, db = manipulate_biosqldb.load_db(biodb)

        sys.stdout.write("creating locus_tag2seqfeature_id")
        locus_tag2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)

        sys.stdout.write("creating protein_id2seqfeature_id")
        protein_id2seqfeature_id = manipulate_biosqldb.protein_id2seqfeature_id_dict(server, biodb)

        blastnr2biosql(locus_tag2seqfeature_id,
                        protein_id2seqfeature_id,
                        biodb,
                        args.n_procs,
                        mysql_host,
                        mysql_user,
                        mysql_pwd,
                        mysql_db,
                        args.input_blast)

    if args.update_taxo_table:
        import os
        sqlpsw = os.environ['SQLPSW']
        update_blastnr_taxonomy_table('blastnr_%s' % biodb,  mysql_pwd=sqlpsw)
        update_blastnr_taxonomy_table('blast_swissprot_%s' % biodb,  mysql_pwd=sqlpsw)
        #update_blastnr_taxonomy_table('blast_swissprot_%s' % biodb)
