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


def get_swissprot_annotation(accession_list):

    import urllib
    from urllib.error import URLError

    link = "http://www.uniprot.org/uniprot/?query=id:%s&columns=id,taxon,annotation score,protein names,genes,organism&format=tab" % ('+OR+'.join(accession_list))
    link = link.replace(' ', '%20')

    req = urllib.request.Request(link)
    try:
        page = urllib.request.urlopen(req)
        data = page.read().decode('utf-8').split('\n')
        rows = [i.rstrip().split('\t') for i in data]
        accession2data = {}
        for n,row in enumerate(rows):
            if n==0:
                continue
            elif len(row)<5:
                continue
            else:
                accession2data[row[0]] = row[1:]
        return accession2data
    except URLError:
        print ('echec')
        print (link)
        return (False)




def load_blastswissprot_file_into_db(locus_tag2taxon_id,
                                locus_tag2seqfeature_id,
                                locus_tag2bioentry_id,
                                mysql_host,
                                mysql_user,
                                mysql_pwd,
                                mysql_db,
                                input_blast_files,
                                biodb):


    '''

    '''

    import plastnr2sqltable
    import MySQLdb
    from chlamdb.biosqldb import manipulate_biosqldb
    import time

    n_file = 0
    for one_blast_file in input_blast_files:
        n_file +=1

        server, db = manipulate_biosqldb.load_db(biodb)
        conn = server.adaptor.conn
        cursor = server.adaptor.cursor

        with open(one_blast_file, 'r') as f:
            print ('Loading', n_file, one_blast_file, '...')
            input_file = [i.rstrip().split('\t') for i in f]
            # sp|P0ADG4|SUHB_ECOLI
            sp_accessions = [i[1].split('|')[1] for i in input_file]

            hit_lists = _chunks(sp_accessions, 300)

            print ('getting annotation')
            accession2annotation = {}
            for n, one_list in enumerate(hit_lists):
                if n % 100 == 0:
                    print ('%s lists/ %s' % (n, len(hit_lists)))
                temp_dico = get_swissprot_annotation(one_list)
                if temp_dico:
                    accession2annotation.update(temp_dico)
                else:
                    print ('waiting 20 sec...')
                    time.sleep(20)
                    temp_dico = get_swissprot_annotation(one_list)
                    if temp_dico:
                        accession2annotation.update(temp_dico)
                    else:
                        print ('1-waiting another 60 sec...')
                        time.sleep(60)
                        temp_dico = get_swissprot_annotation(one_list)
                        if temp_dico:
                            accession2annotation.update(temp_dico)
                        else:
                            print ('2-waiting another 60 sec...')
                            print (one_list)
                            time.sleep(60)
                            temp_dico = get_swissprot_annotation(one_list)                            
                        

                        

            all_taxid = [str(accession2annotation[i][0]) for i in accession2annotation]

            sql = 'select taxon_id from blastnr_blastnr_taxonomy;'
            taxid_in_db = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql,)]
            missing_taxons = list(set(all_taxid)-set(taxid_in_db))

            print ('downloading taxonomy data for %s taxons...' % len(missing_taxons))
            plastnr2sqltable.insert_taxons_into_sqldb(missing_taxons, mysql_pwd=mysql_pwd)

            print ("getting taxid2superkingdom")
            sql = 'select taxon_id,superkingdom from blastnr_blastnr_taxonomy'
            cursor.execute(sql,)
            taxon_id2superkingdom = manipulate_biosqldb.to_dict(cursor.fetchall())

            print ("getting locus2protein_length")
            sql = 'select locus_tag,char_length(translation) from orthology_detail' % biodb
            cursor.execute(sql,)
            locus_tag2protein_length = manipulate_biosqldb.to_dict(cursor.fetchall())

            print ('loading blast results into database...')
            for n, line in enumerate(input_file):

                sp_accession = line[1].split('|')[1]
                # taxon
                # annotation score
                # protein names
                # genes
                # organism

                query_accession = line[0]
                annot = accession2annotation[sp_accession]
                subject_taxon_id = annot[0]
                subject_annot_score = annot[1].split(' ')[0]
                subject_description = annot[2]
                subject_genes = annot[3]
                subject_organism = annot[4]
                subject_kingdom = taxon_id2superkingdom[subject_taxon_id]

                seqfeature_id = locus_tag2seqfeature_id[query_accession]

                evalue = line[10]
                percent_identity = float(line[2])
                gaps = int(line[5])
                length = int(line[3])
                query_start = int(line[6])
                query_end = int(line[7])
                query_cov = round(((query_end-query_start)/float(locus_tag2protein_length[query_accession]))*100,2)
                subject_start = int(line[8])
                subject_end = int(line[9])
                bit_score = float(line[11])

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
                5 subject_accession varchar(200)
                6 subject_kingdom varchar(200)
                7 subject_scientific_name TEXT(2000000)
                8 subject_taxid INT
                9 subject_title VARCHAR(2000)
                10 evalue varchar(200)
                11 bit_score float
                12 percent_identity float
                13 gaps int
                14 length int
                15 query_start int
                16 query_end int
                17 query_cov float
                18 subject_start int
                19 subject_end
                20 genes
                21 annot score
                '''

                values = '(%s,%s,%s,%s,"%s","%s","%s",%s,"%s","%s",%s,%s,%s,%s,%s,%s,%s,%s,%s,"%s",%s);' % (locus_tag2taxon_id[query_accession],
                                                                                                    locus_tag2bioentry_id[query_accession],
                                                                                                    seqfeature_id,
                                                                                                    hit_n,
                                                                                                    sp_accession,
                                                                                                    subject_kingdom,
                                                                                                    subject_organism,
                                                                                                    subject_taxon_id,
                                                                                                    subject_description,
                                                                                                    evalue,
                                                                                                    bit_score,
                                                                                                    percent_identity,
                                                                                                    gaps,
                                                                                                    length,
                                                                                                    query_start,
                                                                                                    query_end,
                                                                                                    query_cov,
                                                                                                    subject_start,
                                                                                                    subject_end,
                                                                                                    subject_genes,
                                                                                                    subject_annot_score
                                                                                                    )
                try:
                    sql = 'insert into blast_swissprot_%s values %s' % (biodb,values)
                    cursor.execute(sql,)
                    conn.commit()
                except:
                    print ('problem with sql query')
                    print (sql)

def create_sql_blast_swissprot_tables(db_name, mysql_host, mysql_user, mysql_pwd, mysql_db='blastnr'):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(db_name)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    '''
    1 query_taxon_id int
    2 query_bioentry_id INT
    3 seqfeature_id INT
    4 hit_number int
    5 subject_accession varchar(200)
    6 subject_kingdom varchar(200)
    7 subject_scientific_name TEXT(2000000)
    8 subject_taxid INT
    9 subject_title VARCHAR(2000)
    10 evalue varchar(200)
    11 bit_score float
    12 percent_identity float
    13 gaps int
    14 length int
    15 query_start int
    16 query_end int
    17 query_cov float
    18 subject_start int
    19 subject_end
    20 genes
    21 annot score
    '''

    sql_plast = 'CREATE TABLE IF NOT EXISTS blast_swissprot_%s (query_taxon_id INT,' \
                ' query_bioentry_id INT,' \
                ' seqfeature_id INT,' \
                ' hit_number int,' \
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
                ' query_cov float,' \
                ' subject_start int,' \
                ' subject_end int,' \
                ' genes TEXT,' \
                ' annot_score int,' \
                ' INDEX query_taxon_id (query_taxon_id),' \
                ' INDEX query_bioentry_id (query_bioentry_id),' \
                ' INDEX hit_number (hit_number),' \
                ' INDEX seqfeature_id (seqfeature_id),' \
                ' INDEX subject_taxid(subject_taxid))' % (db_name)

    try:

        cursor.execute(sql_plast)
        print ('sql hits ok')
        conn.commit()
    except:
        print (sql_plast)
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





def blastswiss2biosql( locus_tag2seqfeature_id,
                    db_name,
                    n_procs,
                    mysql_host,
                    mysql_user,
                    mysql_pwd,
                    mysql_db,
                    *input_blast_files):

    import numpy
    from multiprocessing import Process

    create_sql_blast_swissprot_tables(db_name,mysql_host,mysql_user,mysql_pwd,mysql_db)

    # load blast data into blastnr_ db_name table
    n_cpu = n_procs
    n_poc_per_list = int(numpy.ceil(len(input_blast_files)/float(n_cpu)))
    query_lists = _chunks(input_blast_files, n_poc_per_list)
    print ('n lists: %s' % len(query_lists))



    sql = 'create table if not exists blastnr.gi2taxon_and_description (gi INT, taxon_id INT, description TEXT,' \
          ' taxonomy TEXT, source TEXT);'
    server.adaptor.execute(sql,)

    print ('get locus2taxon_id')
    sql = 'select locus_tag, taxon_id from orthology_detail' % db_name
    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    print ('get locus2bioentry')
    sql2 = 'select locus_tag,bioentry_id from biodatabase t1 ' \
           ' inner join bioentry as t2 on t1.biodatabase_id=t2.biodatabase_id' \
           ' inner join orthology_detail t3 on t2.accession=t3.accession where t1.name="%s"' % (db_name,db_name)

    locus_tag2bioentry_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

    
    if len(input_blast_files)>n_poc_per_list:

        procs = []
        for one_list in query_lists:
            proc = Process(target=load_blastswissprot_file_into_db, args=(locus_tag2taxon_id,
                                                                        locus_tag2seqfeature_id,
                                                                        locus_tag2bioentry_id,
                                                                        mysql_host,
                                                                        mysql_user,
                                                                        mysql_pwd,
                                                                        mysql_db,
                                                                        one_list,
                                                                        biodb))
            procs.append(proc)
            proc.start()

        # Wait for all worker processes to finish
        for proc in procs:
            proc.join()
    else:
        print ('Only %s input blast file(s) for %s cpus, not working in paralell mode' % (len(input_blast_files),n_poc_per_list))
        load_blastswissprot_file_into_db(locus_tag2taxon_id,
                                locus_tag2seqfeature_id,
                                locus_tag2bioentry_id,
                                mysql_host,
                                mysql_user,
                                mysql_pwd,
                                mysql_db,
                                input_blast_files,
                                biodb)

    sys.stdout.write("done!")





if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    import sys
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_blast', type=str, help="input blast tab files", nargs='+')
    parser.add_argument("-d", '--mysql_database', type=str, help="Biosql biodatabase name")
    parser.add_argument("-t", '--create_tables', action='store_true', help="Create SQL tables")
    parser.add_argument("-p", '--n_procs', type=int, help="Number of threads to use (default=8)", default=8)
    parser.add_argument("-c", '--clean_tables', action='store_true', help="delete all sql tables")
    parser.add_argument("-l", '--load_tables', action='store_true', help="load tab files into biodatabase")
    parser.add_argument("-f", '--filter_n_hits', action='store_true', help="filter_n_hits (max 100 hits/locus)")

    args = parser.parse_args()


    mysql_host = 'localhost'
    mysql_user = 'root'

    mysql_pwd = os.environ['SQLPSW']
    mysql_db = 'blastnr'

    biodb = args.mysql_database

    if args.load_tables:

        server, db = manipulate_biosqldb.load_db(biodb)

        sys.stdout.write("creating locus_tag2seqfeature_id")
        locus_tag2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)

        sys.stdout.write("creating protein_id2seqfeature_id")
        protein_id2seqfeature_id = manipulate_biosqldb.protein_id2seqfeature_id_dict(server, biodb)

        blastswiss2biosql(locus_tag2seqfeature_id,
                        biodb,
                        args.n_procs,
                        mysql_host,
                        mysql_user,
                        mysql_pwd,
                        mysql_db,
                        *args.input_blast)


