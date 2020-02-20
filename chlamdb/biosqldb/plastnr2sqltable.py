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

def create_protein_accession2taxon_id_table():
    sql = 'CREATE TABLE ncbi_nr.protein_id2taxon_id (' \
    ' protein_id VARCHAR(60) UNIQUE,' \
    ' taxon_id INT);'
    print (sql)

def insert_hit(conn,
               cursor,
               db_name,
               locus_tag2accession,
               hit_n,
               query_accession,
               locus_tag,
               subject_accession,
               subject_gi,
               subject_kingdom,
               subject_scientific_names,
               subject_taxids,
               subject_title):

        sql_blast_hit = 'INSERT INTO blastnr_hits_%s_%s(hit_number, ' \
                    'query_accession, ' \
                    'locus_tag, ' \
                    'subject_gi, ' \
                    'subject_accession, ' \
                    'subject_kingdom, ' \
                    'subject_scientific_name, ' \
                    'subject_taxid, ' \
                    'subject_title)' \
                    ' values (%s,"%s","%s",%s,"%s","%s","%s","%s","%s");' % (db_name,
                                                    locus_tag2accession[locus_tag],
                                                    hit_n,
                                                    query_accession,
                                                    locus_tag,
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
        index = conn.insert_id()
        conn.commit()
        return index

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

    import accession2taxon_id
    from datetime import datetime

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
            'evalue,' \
            'bit_score,' \
            'percent_identity, ' \
            'gaps, ' \
            'length, ' \
            'query_start, ' \
            'query_end, ' \
            'subject_start, ' \
            'subject_end)' \
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
            print ('Loading', n_file, one_blast_file, '...')
            input_file = [i.rstrip().split('\t') for i in f]

            input_file_hit_gi = [i[1].split('|')[1] for i in input_file]

            

            print ('getting protein 2 taxon id for %s proteins' % len(input_file_hit_gi))
            gi2taxon_id = {}
            gi2descriptions = {}

            id_lists = _chunks(input_file_hit_gi, 300)
            for i, one_list in enumerate(id_lists):
                #print i, "/", len(id_lists)
                if i % 100 == 0:
                    print (i, " lists /", len(id_lists))
                    time.sleep(60)
                gi2taxon_id.update(accession2taxon_id.gi2taxon_id(one_list, "protein"))
                gi2descriptions.update(accession2taxon_id.gi2description(one_list, "protein"))
            #print gi2descriptions.keys()
            #print 'getting protein 2 description dir %s proteins' % len(input_file_hit_accessions)

            print ('loading blast results into database...')
            for n, line in enumerate(input_file):
                #print n, line[0]
                # qgi qacc sgi sacc sscinames sskingdoms staxids evalue nident pident positive gaps length qstart qend qcovs sstart send sstrand stitle
                if n%1000 == 0:
                    print (time.ctime() + ': ' + str(round((float(n)/len(input_file))*100,2)) + '%...')

                query_accession = line[0]

                try:
                    seqfeature_id = protein_id2seqfeature_id[query_accession]
                except KeyError:
                    seqfeature_id = locus_tag2seqfeature_id[query_accession]
                locus_tag = seqfeature_id2locus_tag[str(seqfeature_id)]

                if n != 0:
                    query_accession_previous_line = input_file[n-1][0]
                    try:
                        seqfeature_id_previous_line = protein_id2seqfeature_id[query_accession_previous_line]
                    except KeyError:
                        seqfeature_id_previous_line = locus_tag2seqfeature_id[query_accession_previous_line]
                    locus_tag_previous_line = seqfeature_id2locus_tag[str(seqfeature_id_previous_line)]

                    # check if new from the same accession or new accession
                    if locus_tag2accession[locus_tag_previous_line] != locus_tag2accession[locus_tag]:
                            # insert previous hsps data
                            sql_hsps = sql_blast_hsp_head % (db_name, locus_tag2accession[locus_tag_previous_line], values[0:-1])

                            #print sql_hsps
                            try:
                                cursor.execute(sql_hsps)
                                conn.commit()
                                values = ''
                            except:
                                print (sql_hsps)
                    #else:
                    #    
                    #    print 'previous loocus accession: %s --- new locus: %s' % (locus_tag2accession[locus_tag_previous_line], locus_tag2accession[locus_tag])

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

                subject_accession = line[1].split("|")[3]
                #query_gi = int(line[1].split("|")[1]) # remove
                subject_gi = int(line[1].split("|")[1]) # parse it from plastp result
                try:
                    subject_kingdom = gi2descriptions[str(subject_gi)]['taxonomy'][0]  # get from accession2annotations
                except:
                    subject_kingdom = '-'
                subject_scientific_names = gi2descriptions[str(subject_gi)]['source'] # get from accession2annotations
                subject_taxids = gi2taxon_id[subject_gi] # get from accession2taxon_id
                subject_title = gi2descriptions[str(subject_gi)]['description'] # get from accession2annotations

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
                                       locus_tag,
                                       subject_accession,
                                       subject_gi,
                                       subject_kingdom,
                                       subject_scientific_names,
                                       subject_taxids,
                                       subject_title)

                    values += '(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s),' % (hit_id,
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

                # check if new query accession)
                elif line[0] != input_file[n-1][0]:
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
                                       locus_tag,
                                       subject_accession,
                                       subject_gi,
                                       subject_kingdom,
                                       subject_scientific_names,
                                       subject_taxids,
                                       subject_title)

                    values += '(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s),' % (hit_id,
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

                # same query
                else:
                    # same hit ==> multiple hsps ==> hit id remain the same
                    if line[1] == input_file[n-1][1]:
                        # same hit
                        hsp_n+=1
                        values += '(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s),' % (hit_id,
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
                                       locus_tag,
                                       subject_accession,
                                       subject_gi,
                                       subject_kingdom,
                                       subject_scientific_names,
                                       subject_taxids,
                                       subject_title)


                        # initiate new value string
                        values += '(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s),' % (hit_id,
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
    from MySQLdb import ProgrammingError
    conn = MySQLdb.connect(host=mysql_host, # your host, usually localhost
                                user=mysql_user, # your username
                                passwd=mysql_pwd, # your password
                                db=mysql_db) # name of the data base
    cursor = conn.cursor()



    
    for accession in accession_list:
        print ('Loading...', accession)

        sql1 = 'select nr_hit_id, subject_taxid from blastnr_hits_%s_%s' % (biodb, accession)
        cursor.execute(sql1)
        nr_hit_id2taxon_ids = manipulate_biosqldb.to_dict(cursor.fetchall())
        for nr_hit in nr_hit_id2taxon_ids.keys():
            # new table for the multiples taxons ids/hits
            for taxon in nr_hit_id2taxon_ids[nr_hit].split(';'):
                if taxon != 'N/A':
                    try:
                        sql = 'INSERT INTO blastnr_hits_taxonomy_%s_%s (' \
                            'nr_hit_id, subject_taxon_id) values (%s, %s)' % (biodb,
                                                                      accession,
                                                                      nr_hit,
                                                                      taxon)
                    except ProgrammingError:
                        print ('PROBLEM WITH:')
                        print (sql)
                #try:
                    cursor.execute(sql)
                    conn.commit()
                #except:
                #    print "problem with"
                #    print sql
                else:
                    print ('N/A taxon for hit %s' % nr_hit)


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

    server, db = manipulate_biosqldb.load_db(db_name)

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
                        ' bit_score float,' \
                        ' percent_identity float,' \
                        ' gaps int,' \
                        ' length int,' \
                        ' query_start int,' \
                        ' query_end int,' \
                        ' subject_start int,' \
                        ' subject_end int,' \
                        ' INDEX nr_hit_id (nr_hit_id))' % (db_name, accession)

            sql_blast_hsps2 = 'ALTER TABLE blastnr_blastnr_hsps_%s ADD CONSTRAINT fk_blast_hsp_hit_id_%s ' \
                              'FOREIGN KEY (nr_hit_id) REFERENCES blastnr_hits_%s_%s(nr_hit_id);' % (db_name, accession, accession, db_name, accession)


            sql_blast_hits = 'CREATE TABLE IF NOT EXISTS blastnr_hits_%s_%s (nr_hit_id INT AUTO_INCREMENT PRIMARY KEY, ' \
                            ' query_accession varchar(200),' \
                            ' hit_number int,' \
                            ' locus_tag VARCHAR(100), ' \
                            ' subject_gi int,' \
                            ' subject_accession varchar(200),' \
                            ' subject_kingdom varchar(200),' \
                            ' subject_scientific_name TEXT(2000000), ' \
                            ' subject_taxid TEXT(2000000),' \
                            ' subject_title VARCHAR(2000),' \
                            ' FOREIGN KEY (locus_tag) REFERENCES orthology_detail(locus_tag),' \
                            ' INDEX locus_tag (locus_tag))' % (db_name, accession, db_name)


            sql_blast_taxonomy = 'CREATE TABLE blastnr_blastnr_hits_taxonomy_%s (nr_hit_id INT, ' \
                                 ' subject_taxon_id int,' \
                                 ' INDEX nr_hit_id (nr_hit_id),' \
                                 ' INDEX subject_taxon_id (subject_taxon_id))' % (db_name, accession)

            sql_blast_taxonomy2 = 'ALTER TABLE blastnr_blastnr_hits_taxonomy_%s ADD CONSTRAINT fk_blast_hit_id_%s ' \
                                  'FOREIGN KEY (nr_hit_id) REFERENCES blastnr_hits_%s_%s(nr_hit_id);' % (db_name, accession, accession, db_name, accession)

            try:

                cursor.execute(sql_blast_hits)
                print ('sql hits ok')
                conn.commit()
            except:
                print (sql_blast_hits)
                print ('not created')

            try:

                cursor.execute(sql_blast_hsps)
                conn.commit()
                print ('sql hsps1 ok')
                cursor.execute(sql_blast_hsps2)
                print ('sql hsps2 ok')
                conn.commit()
            except:
                print (sql_blast_hsps)
                print ('not created')

            try:

                cursor.execute(sql_blast_taxonomy)
                conn.commit()
                cursor.execute(sql_blast_taxonomy2)
                print ("sql_taxonomy ok ")
                conn.commit()
            except:
                print (sql_blast_taxonomy)
                print ('not created')

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


    if main_blastnr_table:
        cursor.execute(sql_taxonomy1)
        conn.commit()


def filter_blast_number(one_blast_file, out_name,max_hits_per_locus=100):
    locus2count = {}
    g = open(out_name, 'w')

    with open(one_blast_file, 'r') as f:
        for row in f:
            locus = row.rstrip().split("\t")[0]
            #print 'locus', locus
            if locus not in locus2count:
                locus2count[locus] = 1
            else:
                if locus2count[locus] <=max_hits_per_locus:
                    g.write(row)
                    locus2count[locus] += 1


def insert_taxons_into_sqldb(taxon_id_list,
                             chunk_size=300,
                             mysql_host = 'localhost',
                             mysql_user = 'root',
                             mysql_pwd = 'wnkonwn',
                             mysql_db = 'blastnr'):
    import time
    import MySQLdb
    import sequence_id2scientific_classification
    import re

    conn = MySQLdb.connect(host=mysql_host, # your host, usually localhost
                                user=mysql_user, # your username
                                passwd=mysql_pwd, # your password
                                db=mysql_db) # name of the data base
    cursor = conn.cursor()

    taxid2classification = {}

    id_lists = _chunks(taxon_id_list, chunk_size)
    for i, one_list in enumerate(id_lists):
        print (i, "/", len(id_lists))
        if i % 100 == 0 and i != 0:
            time.sleep(60)
        taxid2classification.update(sequence_id2scientific_classification.taxon_id2scientific_classification(one_list))

    print ('Number of taxon id retrieved:', len(taxid2classification.keys()))

    print ('Updating blastnr_taxonomy table with %s new taxons' % str(len(taxon_id_list)))
    for taxon_id in taxon_id_list:
            if taxon_id == 'N/A':
                continue

            sql = 'INSERT INTO blastnr_taxonomy(taxon_id) values (%s)' % (taxon_id)
            try:

                cursor.execute(sql)
                conn.commit()
            except MySQLdb.IntegrityError:
                print ('Taxon %s already in database' % str(taxon_id))
                continue
            try:
                for rank in taxid2classification[taxon_id].keys():

                    sql_id = re.sub(' ', '_', rank)

                    value = taxid2classification[taxon_id][rank]
                    sql = 'UPDATE blastnr_taxonomy SET `%s`="%s" where taxon_id=%s' % (sql_id, value, taxon_id)
                    #print sql
                    try:
                        cursor.execute(sql)
                        conn.commit()
                    except:
                        print ('could not insert rank', sql)
            except KeyError:
                print ('Could not add the following taxon: %s, trying again to get the data...' % taxon_id)
                temp_dico = sequence_id2scientific_classification.taxon_id2scientific_classification([taxon_id])

                try:
                    for rank in temp_dico[taxon_id].keys():

                        sql_id = re.sub(' ', '_', rank)

                        value = temp_dico[taxon_id][rank]
                        sql = 'UPDATE blastnr_taxonomy SET `%s`="%s" where taxon_id=%s' % (sql_id, value, taxon_id)

                        cursor.execute(sql)
                        conn.commit()
                    print ('sucess!')
                except KeyError:
                    print ('Could not add %s' % taxon_id)
                    print ("temp_dico", temp_dico)


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


    from time import gmtime, strftime
    import re
    import MySQLdb
    import time
    import accession2taxon_id
    from chlamdb.biosqldb import manipulate_biosqldb

    conn = MySQLdb.connect(host=mysql_host, # your host, usually localhost
                                user=mysql_user, # your username
                                passwd=mysql_pwd, # your password
                                db=mysql_db) # name of the data base
    cursor = conn.cursor()

    sql1 = 'select gi, taxon_id from gi2taxon_id'
    cursor.execute(sql1,)
    gi2taxon_id_database = manipulate_biosqldb.to_dict(cursor.fetchall())
    sql = 'SELECT taxon_id from blastnr_taxonomy'
    all_taxon_ids = []
    cursor.execute(sql,)
    database_taxon_ids = [str(i[0]) for i in cursor.fetchall()]

    print ("Number of taxons into database: ", len(database_taxon_ids))

    protein_gi_list = []
    
    for one_blast_file in input_blast_files:

                print ('blast file %s' % one_blast_file)
                with open(one_blast_file, 'r') as f:

                    #hit_column =
                    all_gi_ids = [i.rstrip().split("\t")[1].split('|')[1] for i in f]
                    protein_gi_list+=all_gi_ids

    nr_protein_gi_list = list(set(protein_gi_list))
    gi_def_list = []
    for gi in nr_protein_gi_list:
        # if gi already into database with corresponding taxon_id, skip
        try:
            if gi2taxon_id_database[gi] in database_taxon_ids:
                continue
            else:
                gi_def_list.append(gi)
        except KeyError:
            gi_def_list.append(gi)
    print ('getting gi 2 taxon id for %s proteins' % len(gi_def_list))
    gi2taxon_id = {}

    id_lists = _chunks(gi_def_list, 300)
    for i, one_list in enumerate(id_lists):
        #print i, "/", len(id_lists)
        if i % 100 == 0:
            print (i, "lists ok / %s lists" % len(id_lists), strftime("%Y-%m-%d %H:%M:%S", gmtime()))
            time.sleep(60)
            gi2taxon_id.update(accession2taxon_id.gi2taxon_id(one_list, "protein"))
            #print "protein_accession2taxon_id", protein_accession2taxon_id
            temp_subject_taxids = [str(i) for i in gi2taxon_id.values()]

            for i in temp_subject_taxids:
                if i not in database_taxon_ids and i not in all_taxon_ids:
                    all_taxon_ids.append(i)

    # add new values to gi2taxon_id                
    for gi in gi2taxon_id:
        sql = 'insert into gi2taxon_id values(%s,%s)' % (gi, gi2taxon_id[gi])
        cursor.execute(sql,)
    conn.commit()
                    
    print ('Number of new taxons:', len(all_taxon_ids))
    print ('i.e: ', all_taxon_ids[1:10])
    time.sleep(60)

    print ('Fetching ncbi taxonomy... for %s taxons' % str(len(all_taxon_ids)))

    # subdivide the taxon list in smaller lists, otherwise NCBI will limit the results to? 10000 (observed once only)
    insert_taxons_into_sqldb(all_taxon_ids, 300, mysql_pwd=mysql_pwd)




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
    import biosql_own_sql_tables
    print ('host1', mysql_host)
    print ('user1', mysql_user)

    print ("add eventual new taxons to the main blastnr taxonomy table")
    #update2biosql_blastnr_table(mysql_host, mysql_user, mysql_pwd, mysql_db, *input_blast_files)

    # load blast data into blastnr_ db_name table
    n_cpu = n_procs
    n_poc_per_list = int(numpy.ceil(len(input_blast_files)/float(n_cpu)))
    query_lists = _chunks(input_blast_files, n_poc_per_list)
    print ('n lists: %s' % len(query_lists))
    

    if len(input_blast_files)>n_poc_per_list:

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
    else:
        print ('Only %s input blast file(s), not working in paralell mode' % (len(input_blast_files)))
        _load_blastnr_file_into_db(seqfeature_id2locus_tag,
                                locus_tag2seqfeature_id,
                                protein_id2seqfeature_id,
                                locus_tag2accession,
                                db_name,
                                mysql_host,
                                mysql_user,
                                mysql_pwd,
                                mysql_db,
                                input_blast_files)

    blastnr2biodb_taxonomic_table(db_name, locus_tag2accession, mysql_host, mysql_user, mysql_pwd, mysql_db, n_procs)
    sys.stdout.write("cleaning multispecies blastnr record...")
    biosql_own_sql_tables.clean_multispecies_blastnr_record(db_name, create_new_sql_tables=True)

def del_blastnr_table_content(db_name):
    server, db = manipulate_biosqldb.load_db(db_name)
    sql = 'select accession from bioentry' \
      ' inner join biodatabase on bioentry.biodatabase_id=biodatabase.biodatabase_id where biodatabase.name="%s"' % db_name

    all_accessions = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    for accession in all_accessions:
        sql1 = 'DROP TABLE IF EXISTS blastnr_blastnr_hsps_%s' % (db_name, accession)
        sql2 = 'DROP TABLE IF EXISTS blastnr_blastnr_hits_%s' % (db_name, accession)
        sql3 = 'DROP TABLE IF EXISTS blastnr_blastnr_hits_taxonomy_%s' % (db_name, accession)
        sql4 = 'DROP TABLE IF EXISTS blastnr_blastnr_hits_taxonomy_filtered_%s' % (db_name, accession)


        print (sql1)
        server.adaptor.execute(sql1)
        server.adaptor.commit()
        print (sql3)
        server.adaptor.execute(sql3)
        server.adaptor.commit()
        print (sql4)
        server.adaptor.execute(sql4)
        server.adaptor.commit()
        print (sql2)
        server.adaptor.execute(sql2)
        server.adaptor.commit()




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
    parser.add_argument("-p", '--n_procs', type=int, help="Number of threads to use (default=8)", default=8)
    parser.add_argument("-c", '--clean_tables', action='store_true', help="delete all sql tables")
    parser.add_argument("-l", '--load_tables', action='store_true', help="load tab files into biodatabase")
    parser.add_argument("-f", '--filter_n_hits', action='store_true', help="filter_n_hits (max 100 hits/locus)")
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



    if args.filter_n_hits:
        for i in args.input_blast:
            outname = i.split('.')[0] + '_filtered_100.tab'
            print (outname)
            filter_blast_number(i,outname,100)


    if args.create_tables:



        create_sql_blastnr_tables(args.mysql_database, mysql_host, mysql_user, mysql_pwd, mysql_db, main_blastnr_table=True, alternate_tables=True)


    if args.load_tables:


        server, db = manipulate_biosqldb.load_db(biodb)

        sys.stdout.write("creating locus_tag2seqfeature_id")
        locus_tag2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)

        sys.stdout.write("creating protein_id2seqfeature_id")
        protein_id2seqfeature_id = manipulate_biosqldb.protein_id2seqfeature_id_dict(server, biodb)

        sys.stdout.write("getting seqfeature_id2locus_tag")
        seqfeature_id2locus_tag = manipulate_biosqldb.seqfeature_id2locus_tag_dico(server, biodb)

        sys.stdout.write("getting locus_tag2accession")
        locus_tag2accession = manipulate_biosqldb.locus_tag2accession(server, args.mysql_database)


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


