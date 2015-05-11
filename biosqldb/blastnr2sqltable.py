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


def insert_hit(conn,
               cursor,
               db_name,
               locus_tag2accession,
               hit_n,
               query_accession,
               line,
               locus_tag):

        import re
        subject_accession = line[3]
        query_gi = int(line[0])
        subject_gi = int(line[2])

        subject_kingdom = line[5]
        subject_scientific_names = line[4]
        subject_taxids = line[6]
        subject_title = re.sub('\"', '', line[19]).split('[')[0]






        sql_blast_hit = 'INSERT INTO blastnr_hits_%s_%s(hit_number, ' \
                    'query_accession, ' \
                    'locus_tag, ' \
                    'query_gi, ' \
                    'subject_gi, ' \
                    'subject_accession, ' \
                    'subject_kingdom, ' \
                    'subject_scientific_name, ' \
                    'subject_taxid, ' \
                    'subject_title)' \
                    ' values (%s,"%s","%s",%s,%s,"%s","%s","%s","%s","%s");' % (db_name,
                                                    locus_tag2accession[locus_tag],
                                                    hit_n,
                                                    query_accession,
                                                    locus_tag,
                                                    query_gi,
                                                    subject_gi,
                                                    subject_accession,
                                                    subject_kingdom,
                                                    subject_scientific_names,
                                                    subject_taxids,
                                                    subject_title
                                                    )
        '''
        sql_id = 'SELECT `AUTO_INCREMENT`' \
                 ' FROM  INFORMATION_SCHEMA.TABLES' \
                 ' WHERE TABLE_SCHEMA = "blastnr"' \
                 ' AND   TABLE_NAME   = "blastnr_hits_%s_%s";' % (db_name, locus_tag2accession[locus_tag])
        '''

        cursor.execute(sql_blast_hit)

        conn.commit()
        return conn.insert_id()

        #except:
        #print 'problem with:'
        #print sql_blast_hit
        #import sys
        #sys.exit()


def _load_blastnr_file_into_db(seqfeature_id2locus_tag,
                                locus_tag2seqfeature_id,
                                protein_id2seqfeature_id,
                                locus_tag2accession,
                                db_name,
                                mysql_host,
                                mysql_user,
                                mysql_pwd,
                                mysql_db,
                                input_blast_files):

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
    import re

    import MySQLdb


    sql_blast_hsp_head = 'INSERT INTO blastnr_hsps_%s_%s(nr_hit_id, ' \
            'evalue, ' \
            'n_identical, ' \
            'percent_identity, ' \
            'positive, ' \
            'gaps, ' \
            'length, ' \
            'query_start, ' \
            'query_end, ' \
            'query_cov, ' \
            'subject_start, ' \
            'subject_end, ' \
            'subject_strand)' \
            ' values %s'


    n_file = 0
    for one_blast_file in input_blast_files:
        n_file +=1
        conn = MySQLdb.connect(host=mysql_host, # your host, usually localhost
                                    user=mysql_user, # your username
                                    passwd=mysql_pwd, # your password
                                    db=mysql_db) # name of the data base
        cursor = conn.cursor()

        values = ''

        with open(one_blast_file, 'r') as f:
            print 'Loading', n_file, one_blast_file, '...'
            input_file = [i.rstrip().split('\t') for i in f]

            print 'loading blast results into database...'
            for n, line in enumerate(input_file):
                # qgi qacc sgi sacc sscinames sskingdoms staxids evalue nident pident positive gaps length qstart qend qcovs sstart send sstrand stitle
                if n%1000 == 0:
                    print time.ctime() + ': ' + str(round((float(n)/len(input_file))*100,2)) + '%...'

                query_accession = line[1].split("|")[3]
                try:
                    seqfeature_id = protein_id2seqfeature_id[query_accession]
                except KeyError:
                    seqfeature_id = locus_tag2seqfeature_id[query_accession]
                locus_tag = seqfeature_id2locus_tag[str(seqfeature_id)]

                if n != 0:
                    query_accession_previous_line = input_file[n-1][1].split("|")[3]
                    try:
                        seqfeature_id_previous_line = protein_id2seqfeature_id[query_accession_previous_line]
                    except KeyError:
                        seqfeature_id_previous_line = locus_tag2seqfeature_id[query_accession_previous_line]
                    locus_tag_previous_line = seqfeature_id2locus_tag[str(seqfeature_id_previous_line)]

                    if locus_tag2accession[locus_tag_previous_line] != locus_tag2accession[locus_tag]:
                            # insert previous hsps data
                            sql_hsps = sql_blast_hsp_head % (db_name, locus_tag2accession[locus_tag_previous_line], values[0:-1])

                            #print sql_hsps
                            cursor.execute(sql_hsps)
                            conn.commit()
                            values = ''

                evalue = line[7]
                n_identical = int(line[8])
                percent_identity = float(line[9])
                positive = int(line[10])
                gaps = int(line[11])
                length = int(line[12])
                query_start = int(line[13])
                query_end = int(line[14])
                query_cov = int(line[15])
                subject_start = int(line[16])
                subject_end = int(line[17])
                subject_strand = line[18]


                # first hit of the file
                if n == 0:
                    blast_n = 1
                    hit_n = 1
                    hsp_n = 1
                    # new hit, insert new entry into blast_hit table


                    hit_id = insert_hit(conn,
                                       cursor,
                                       db_name,
                                       locus_tag2accession,
                                       hit_n,
                                       query_accession,
                                       line,
                                       locus_tag)

                    values += '(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,"%s"),' % (hit_id,
                                                                    evalue,
                                                                    n_identical,
                                                                    percent_identity,
                                                                    positive,
                                                                    gaps,
                                                                    length,
                                                                    query_start,
                                                                    query_end,
                                                                    query_cov,
                                                                    subject_start,
                                                                    subject_end,
                                                                    subject_strand,
                                                                    )

                # check if new query accession)
                elif line[1] != input_file[n-1][1]:
                    # new query => insert its first hit into blastnrdb
                    blast_n += 1
                    hit_n = 1
                    hsp_n = 1

                    hit_id = insert_hit(conn,
                                       cursor,
                                       db_name,
                                       locus_tag2accession,
                                       hit_n,
                                       query_accession,
                                       line,
                                       locus_tag)

                    values += '(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,"%s"),' % (hit_id,
                                                                    evalue,
                                                                    n_identical,
                                                                    percent_identity,
                                                                    positive,
                                                                    gaps,
                                                                    length,
                                                                    query_start,
                                                                    query_end,
                                                                    query_cov,
                                                                    subject_start,
                                                                    subject_end,
                                                                    subject_strand,
                                                                    )

                # same query
                else:
                    # same hit ==> multiple hsps ==> hit id remain the same
                    if line[3] == input_file[n-1][3]:
                        # same hit
                        hsp_n+=1
                        values += '(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,"%s"),' % (hit_id,
                                                                        evalue,
                                                                        n_identical,
                                                                        percent_identity,
                                                                        positive,
                                                                        gaps,
                                                                        length,
                                                                        query_start,
                                                                        query_end,
                                                                        query_cov,
                                                                        subject_start,
                                                                        subject_end,
                                                                        subject_strand,
                                                                        )
                    # different hit of the same query
                    else:
                        hit_n+=1
                        # new hit, insert new entry into blast_hit table

                        hit_id = insert_hit(conn,
                                           cursor,
                                           db_name,
                                           locus_tag2accession,
                                           hit_n,
                                           query_accession,
                                           line,
                                           locus_tag)


                        # initiate new value string
                        values += '(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,"%s"),' % (hit_id,
                                                                        evalue,
                                                                        n_identical,
                                                                        percent_identity,
                                                                        positive,
                                                                        gaps,
                                                                        length,
                                                                        query_start,
                                                                        query_end,
                                                                        query_cov,
                                                                        subject_start,
                                                                        subject_end,
                                                                        subject_strand,
                                                                        )
            # end fo the file, inserting hsps
            sql_hsps = sql_blast_hsp_head % (db_name, locus_tag2accession[locus_tag], values[0:-1])

            #print sql_hsps
            cursor.execute(sql_hsps)
            conn.commit()
            values = ''








                # todo check if blast hit already loaded
                #sql = 'select * from blastnr_%s where (query_accession="%s" and subject_accession="%s" and subject_start=%s and subject_end=%s)' % (biodb, query_accession, subject_accession,subject_start, subject_end)

                #check = server.adaptor.execute_and_fetchall(sql,)

                #if len(check) != 0:
                #    print 'Hit %s VS %s already into database' % (query_accession, subject_accession)
                #    continue







                #print sql_blast_hsp




def _load_taxonomic_data(biodb, mysql_host, mysql_user, mysql_pwd, mysql_db, accession_list):
    '''

    Creating new entries for a chunck of nr_id_list

    :param biodb:
    :param nrhit2taxon_id_list: dictionnary of unique identifier for each blast hsp and the corresponding hits taxon ids (semi colon separeted)
    taxon1;taxon2;taxon3;...
    :param nr_id_list: list of unique hsps hits to add into the sql table
    :return:
    '''

    import MySQLdb
    conn = MySQLdb.connect(host=mysql_host, # your host, usually localhost
                                user=mysql_user, # your username
                                passwd=mysql_pwd, # your password
                                db=mysql_db) # name of the data base
    cursor = conn.cursor()

    for accession in accession_list:

        sql1 = 'select nr_hit_id, subject_taxid from blastnr_hits_%s_%s' % (biodb, accession)
        cursor.execute(sql1)
        nr_hit_id2taxon_ids = manipulate_biosqldb.to_dict(cursor.fetchall())
        for nr_hit in nr_hit_id2taxon_ids.keys():
            # new table for the multiples taxons ids/hits
            for taxon in nr_hit_id2taxon_ids[nr_hit].split(';'):
                if taxon != 'N/A':
                    sql = 'INSERT INTO blastnr_hits_taxonomy_%s_%s (' \
                        'nr_hit_id, subject_taxon_id) values (%s, %s)' % (biodb,
                                                                  accession,
                                                                  nr_hit,
                                                                  taxon)
                #try:
                    cursor.execute(sql)
                    conn.commit()
                #except:
                #    print "problem with"
                #    print sql
                else:
                    print 'N/A taxon for hit %s' % nr_hit


def blastnr2biodb_taxonomic_table(db_name,
                                  locus_tag2accession,
                                  mysql_host,
                                  mysql_user,
                                  mysql_pwd,
                                  mysql_db,
                                  n_procs):
    '''
    Add data to table blastnr_taxonomy_`db_name`
    Foreach hit, list the taxonomic ids (each hit can have multiple taxon id => MULTISPECIES)
    :param db_name: biodatabase name
    :return:
    '''

    import MySQLdb
    conn = MySQLdb.connect(host=mysql_host, # your host, usually localhost
                                user=mysql_user, # your username
                                passwd=mysql_pwd, # your password
                                db=mysql_db) # name of the data base
    cursor = conn.cursor()


    import numpy
    from multiprocessing import Process

    sql = 'select accession from bioentry' \
          ' inner join biodatabase on bioentry.biodatabase_id=biodatabase.biodatabase_id where biodatabase.name="%s"' % db_name

    all_accessions = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    n_cpu = n_procs
    n_poc_per_list = int(numpy.ceil(len(all_accessions)/float(n_cpu)))
    query_lists = _chunks(all_accessions, n_poc_per_list)

    procs = []
    for one_list in query_lists:
        proc = Process(target=_load_taxonomic_data, args=(db_name, mysql_host, mysql_user, mysql_pwd, mysql_db, one_list))
        procs.append(proc)
        proc.start()

    # Wait for all worker processes to finish
    for proc in procs:
        proc.join()



def create_sql_blastnr_tables(db_name, mysql_host, mysql_user, mysql_pwd, mysql_db='blastnr', main_blastnr_table=False, alternate_tables=True):
    import MySQLdb

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'select accession from bioentry' \
          ' inner join biodatabase on bioentry.biodatabase_id=biodatabase.biodatabase_id where biodatabase.name="%s"' % db_name

    all_accessions = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    conn = MySQLdb.connect(host=mysql_host, # your host, usually localhost
                                user=mysql_user, # your username
                                passwd=mysql_pwd, # your password
                                db=mysql_db) # name of the data base
    cursor = conn.cursor()


    if alternate_tables:

        for accession in all_accessions:

            sql_blast_hsps = 'CREATE TABLE IF NOT EXISTS blastnr_hsps_%s_%s (hsp_id INT AUTO_INCREMENT PRIMARY KEY, ' \
                        ' nr_hit_id INT,' \
                        ' evalue varchar(200),' \
                        ' n_identical int,' \
                        ' percent_identity float,' \
                        ' positive int,' \
                        ' gaps int,' \
                        ' length int,' \
                        ' query_start int,' \
                        ' query_end int,' \
                        ' query_cov float,' \
                        ' subject_start int,' \
                        ' subject_end int,' \
                        ' subject_strand VARCHAR(5))' % (db_name, accession)

            sql_blast_hsps2 = 'ALTER TABLE blastnr.blastnr_hsps_%s_%s ADD CONSTRAINT fk_blast_hsp_hit_id_%s ' \
                              'FOREIGN KEY (nr_hit_id) REFERENCES blastnr_hits_%s_%s(nr_hit_id);' % (db_name, accession, accession, db_name, accession)


            sql_blast_hits = 'CREATE TABLE IF NOT EXISTS blastnr_hits_%s_%s (nr_hit_id INT AUTO_INCREMENT PRIMARY KEY, ' \
                            ' query_accession varchar(200),' \
                            ' hit_number int,' \
                            ' locus_tag VARCHAR(100), ' \
                            ' query_gi int,' \
                            ' subject_gi int,' \
                            ' subject_accession varchar(200),' \
                            ' subject_kingdom varchar(200),' \
                            ' subject_scientific_name TEXT(2000000), ' \
                            ' subject_taxid TEXT(2000000),' \
                            ' subject_title VARCHAR(2000))' % (db_name, accession)


            sql_blast_taxonomy = 'CREATE TABLE blastnr.blastnr_hits_taxonomy_%s_%s (nr_hit_id INT, ' \
                                 'subject_taxon_id int)' % (db_name, accession)

            sql_blast_taxonomy2 = 'ALTER TABLE blastnr.blastnr_hits_taxonomy_%s_%s ADD CONSTRAINT fk_blast_hit_id_%s ' \
                                  'FOREIGN KEY (nr_hit_id) REFERENCES blastnr_hits_%s_%s(nr_hit_id);' % (db_name, accession, accession, db_name, accession)

            try:

                cursor.execute(sql_blast_hits)
                print 'sql hits ok'
                conn.commit()
            except:
                print sql_blast_hits
                print 'not created'

            try:

                cursor.execute(sql_blast_hsps)
                conn.commit()
                print 'sql hsps1 ok'
                cursor.execute(sql_blast_hsps2)
                print 'sql hsps2 ok'
                conn.commit()
            except:
                print sql_blast_hsps
                print 'not created'

            try:

                cursor.execute(sql_blast_taxonomy)
                conn.commit()
                cursor.execute(sql_blast_taxonomy2)
                print "sql_taxonomy ok "
                conn.commit()
            except:
                print sql_blast_taxonomy
                print 'not created'

    sql_taxonomy1 = 'CREATE TABLE IF NOT EXISTS blastnr_taxonomy (taxon_id int unique, ' \
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


    if main_blastnr_table:
        cursor.execute(sql_taxonomy1)
        conn.commit()



def update2biosql_blastnr_table(mysql_host, mysql_user, mysql_pwd, mysql_db, *input_blast_files):
    '''
    Update the main blastnr taxonomical table containing taxon ids and their full taxonomical path:
    kindom, class, order, family,...  => include all current (april 2015) NCBI taxonomical ranks

    Full taxonomical path obtained from NCBI using the taxon_id2scientific_classification function using entrez utilities

    If path is not complete, missing rank are indicated with minus symbol '-',
    ie. taxon  51291 (chlamydiales order) only contain:
        cellular organisms; Bacteria; Chlamydiae/Verrucomicrobia group; Chlamydiae; Chlamydiia

    Structure of the table:
    taxon_id, rank1, rank2, rank3,..., rankn

    ==> for each taxon added in the various biodatabases, a new entry is created (if it does not already exist)

    :param server:
    :param input_blast_files:
    :return:
    '''

    import sequence_id2scientific_classification
    import MySQLdb

    print 'host', mysql_host

    print 'user', mysql_user

    conn = MySQLdb.connect(host=mysql_host, # your host, usually localhost
                                user=mysql_user, # your username
                                passwd=mysql_pwd, # your password
                                db=mysql_db) # name of the data base
    cursor = conn.cursor()


    sql = 'SELECT taxon_id from blastnr_taxonomy'
    all_taxon_ids = []
    cursor.execute(sql,)
    database_taxon_ids = [str(i[0]) for i in cursor.fetchall()]
    taxid2classification = {}
    print "Number of taxons into database: ", len(database_taxon_ids)
    for one_blast_file in input_blast_files:
                print 'blast file %s' % one_blast_file
                with open(one_blast_file, 'r') as f:
                    input_file = [i.rstrip().split('\t') for i in f]
                    for n, line in enumerate(input_file):
                        temp_subject_taxids = line[6].split(';')
                        # store taxon_ids
                        for i in temp_subject_taxids:
                            if i not in database_taxon_ids and i not in all_taxon_ids:
                                all_taxon_ids.append(i)
                print 'Number of noew taxons:', len(all_taxon_ids)

    print 'Fetching ncbi taxonomy... for %s taxons' % str(len(all_taxon_ids))

    # subdivide the taxon list in smaller lists, otherwise NCBI will limit the results to? 10000 (observed once only)
    id_lists = _chunks(all_taxon_ids, 5000)
    for one_list in id_lists:
        taxid2classification.update(sequence_id2scientific_classification.taxon_id2scientific_classification(one_list))


    print 'Number of taxon id retrieved:', len(taxid2classification.keys())

    print 'Updating blastnr_taxonomy table with %s new taxons' % str(len(all_taxon_ids))
    for taxon_id in all_taxon_ids:
            if taxon_id == 'N/A':
                continue
            import re
            import MySQLdb
            sql = 'INSERT INTO blastnr_taxonomy(taxon_id) values (%s)' % (taxon_id)
            try:

                cursor.execute(sql)
                conn.commit()
            except MySQLdb.IntegrityError:
                print 'Taxon %s already in database' % str(taxon_id)
                continue
            for rank in taxid2classification[taxon_id].keys():

                sql_id = re.sub(' ', '_', rank)

                value = taxid2classification[taxon_id][rank]
                sql = 'UPDATE blastnr_taxonomy SET `%s`="%s" where taxon_id=%s' % (sql_id, value, taxon_id)
                #print sql

                cursor.execute(sql)
                conn.commit()


def blastnr2biosql(seqfeature_id2locus_tag,
                    locus_tag2seqfeature_id,
                    protein_id2seqfeature_id,
                    locus_tag2accession,
                    db_name,
                    n_procs,
                    mysql_host,
                    mysql_user,
                    mysql_pwd,
                    mysql_db,
                    *input_blast_files):

    import numpy
    from multiprocessing import Process
    print 'host1', mysql_host
    print 'user1', mysql_user

    # add eventual new taxons to the main blastnr taxonomy table
    #update2biosql_blastnr_table(mysql_host, mysql_user, mysql_pwd, mysql_db, *input_blast_files)

    # load blast data into blastnr_ db_name table
    n_cpu = n_procs
    n_poc_per_list = int(numpy.ceil(len(input_blast_files)/float(n_cpu)))
    query_lists = _chunks(input_blast_files, n_poc_per_list)

    procs = []
    for one_list in query_lists:
        proc = Process(target=_load_blastnr_file_into_db, args=(seqfeature_id2locus_tag,
                                                                locus_tag2seqfeature_id,
                                                                protein_id2seqfeature_id,
                                                                locus_tag2accession,
                                                                db_name,
                                                                mysql_host,
                                                                mysql_user,
                                                                mysql_pwd,
                                                                mysql_db,
                                                                one_list))
        procs.append(proc)
        proc.start()

    # Wait for all worker processes to finish
    for proc in procs:
        proc.join()

    blastnr2biodb_taxonomic_table(db_name, locus_tag2accession, mysql_host, mysql_user, mysql_pwd, mysql_db, n_procs)


def del_blastnr_table_content(db_name):
    server, db = manipulate_biosqldb.load_db(db_name)
    sql = 'select accession from bioentry' \
      ' inner join biodatabase on bioentry.biodatabase_id=biodatabase.biodatabase_id where biodatabase.name="%s"' % db_name

    all_accessions = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    for accession in all_accessions:
        sql1 = 'delete from blastnr.blastnr_hsps_%s_%s' % (db_name, accession)
        sql2 = 'delete from blastnr.blastnr_hits_%s_%s' % (db_name, accession)
        sql3 = 'delete from blastnr.blastnr_hits_taxonomy_%s_%s' % (db_name, accession)
        server.adaptor.execute(sql1)
        server.adaptor.execute(sql2)
        server.adaptor.execute(sql3)



if __name__ == '__main__':
    import argparse
    import manipulate_biosqldb
    import json

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_blast', type=str, help="input blast tab files", nargs='+')
    parser.add_argument("-d", '--mysql_database', type=str, help="Biosql biodatabase name")
    parser.add_argument("-t", '--create_tables', action='store_true', help="Create SQL tables")
    parser.add_argument("-p", '--n_procs', type=int, help="Number of threads to use (default=8)", default=8)
    parser.add_argument("-c", '--clean_tables', action='store_true', help="delete all sql tables")



    args = parser.parse_args()

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

    mysql_host = 'localhost'
    mysql_user = 'root'
    mysql_pwd = 'agnathe3'
    mysql_db = 'blastnr'

    biodb = args.mysql_database

    server, db = manipulate_biosqldb.load_db(biodb)

    print "creating locus_tag2seqfeature_id"
    locus_tag2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)

    print "creating protein_id2seqfeature_id"
    protein_id2seqfeature_id = manipulate_biosqldb.protein_id2seqfeature_id_dict(server, biodb)

    print "getting seqfeature_id2locus_tag"
    seqfeature_id2locus_tag = manipulate_biosqldb.seqfeature_id2locus_tag_dico(server, biodb)

    print "getting locus_tag2accession"
    locus_tag2accession = manipulate_biosqldb.locus_tag2accession(server, args.mysql_database)

    if args.create_tables:
        create_sql_blastnr_tables(args.mysql_database, mysql_host, mysql_user, mysql_pwd, mysql_db, main_blastnr_table=False, alternate_tables=True)

    blastnr2biosql(seqfeature_id2locus_tag,
                    locus_tag2seqfeature_id,
                    protein_id2seqfeature_id,
                    locus_tag2accession,
                    biodb,
                    args.n_procs,
                    mysql_host,
                    mysql_user,
                    mysql_pwd,
                    mysql_db,
                    *args.input_blast)


