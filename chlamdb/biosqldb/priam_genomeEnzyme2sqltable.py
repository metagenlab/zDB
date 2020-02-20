#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# create mysql 2 tables from dat file
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2016
# ---------------------------------------------------------------------------





def locus2ko_table(locus_tag2ko_dico,
                   biodatabase,
                   ko_accession2ko_id):

    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodatabase)

    sql2 = 'select locus_tag, seqfeature_id from annotation_seqfeature_id2locus' % biodatabase

    locus2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2))

    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme_seqfeature_id2ko (seqfeature_id INT,' \
           ' ko_id INT, ' \
           ' index ko_id (ko_id),' \
           ' index seqid (seqfeature_id));' % (biodatabase)

    server.adaptor.execute_and_fetchall(sql2,)

    for locus in locus_tag2ko_dico:
        ko = locus_tag2ko_dico[locus]
        ko_id = ko_accession2ko_id[ko]
        seqfeature_id = locus2seqfeature_id[locus]

        sql = 'insert into enzyme_seqfeature_id2ko (seqfeature_id, ko_id) values (%s, %s)' % (biodatabase,
                                                                                                 seqfeature_id,
                                                                                                 ko_id)


        server.adaptor.execute(sql,)
    server.commit()


def parse_blast_koala_output(result_file):
    locus2ko = {}
    with open(result_file, 'r') as f:
        for line in f:
            data = line.rstrip().split('\t')
            if len(data)<2:
                continue
            else:
                locus2ko[data[0]] = data[1]
    return locus2ko






def locus2ec_table(locus_tag2ec_dico, biodatabase):

    from chlamdb.biosqldb import manipulate_biosqldb
    import biosql_own_sql_tables

    server, db = manipulate_biosqldb.load_db(biodatabase)

    #sql = 'select locus_tag, accession from orthology_detail' % biodatabase
    #sql2 = 'select locus_tag, orthogroup from orthology_detail' % biodatabase
    #locus2bioentry_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql))
    #locus2orthogroup = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2))

    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme.seqfeature_id2ec_%s (enzyme_id INT AUTO_INCREMENT PRIMARY KEY,' \
           ' seqfeature_id INT,' \
           ' ec_id INT,' \
           ' FOREIGN KEY(ec_id) REFERENCES enzymes(enzyme_id));' % (biodatabase)

    server.adaptor.execute_and_fetchall(sql2,)

    sql = 'select locus_tag, seqfeature_id from annotation_seqfeature_id2locus' % biodatabase

    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    for locus in locus_tag2ec_dico:
        for ec_data in locus_tag2ec_dico[locus]:

            sql = 'select enzyme_id from enzyme_enzymes where ec="%s"' % ec_data[0]
            ec_id = server.adaptor.execute_and_fetchall(sql,)[0][0]
            seqfeature_id = locus_tag2seqfeature_id[locus]
            sql = 'insert into enzyme.seqfeature_id2ec_%s (seqfeature_id, ec_id) values (%s, %s)' % (biodatabase,
                                                                                             seqfeature_id,
                                                                                             ec_id)
            server.adaptor.execute(sql,)
    server.commit()

def get_microbial_metabolism_in_diverse_environments_kegg01120():

    import urllib2
    import re
    import MySQLdb
    from bs4 import BeautifulSoup
    import os
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="enzyme") # name of the data base
    cursor = conn.cursor()

    conn.set_character_set('utf8')
    cursor.execute('SET NAMES utf8;')
    cursor.execute('SET CHARACTER SET utf8;')
    cursor.execute('SET character_set_connection=utf8;')

    sql = 'CREATE TABLE IF NOT EXISTS microbial_metabolism_map01120 (module_name varchar(400));'

    cursor.execute(sql,)
    conn.commit()

    adress = "http://www.genome.jp/kegg-bin/show_pathway?ko01120"

    html = urllib2.urlopen(adress).read()

    html = re.sub("\&\#","-" ,html)
    soup = BeautifulSoup(html, "html")
    #html = soup.encode('utf-8')#.encode('latin-1') #encode('utf-8') # prettify()

    div = soup.findAll("div", { "class" : "control" })[0]
    input_list = soup.findAll("input")

    list_of_modules = []
    begin = False
    temp_list = []
    for one_input in input_list:
        if 'c_level' in str(one_input):
            if begin == False:
                begin = True
            else:
                list_of_modules.append(temp_list)
                temp_list = []
            continue
        if begin == True:
            module = str(one_input).split('value="')[1][0:-3].split('_')[1]
            sql = 'insert into microbial_metabolism_map01120 values("%s")' % module
            cursor.execute(sql,)

            print (sql)
        else:
            print ('---', one_input)
    conn.commit()


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    import parse_priam_EC

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_priam_files', type=str, help="input interpro csv file", nargs='+', default=False)
    parser.add_argument("-k", '--ko_table', type=str, help="input blastGhost file", default=False)
    parser.add_argument("-d", '--database_name', type=str, help="database name")

    args = parser.parse_args()


    if args.input_priam_files:
        locus2ec={}

        for priam_file in args.input_priam_files:
            locus2ec.update(parse_priam_EC.locus2EC(priam_file))

        locus2ec_table(locus2ec, args.database_name)


    if args.ko_table:
        server, db = manipulate_biosqldb.load_db(args.database_name)
        sql = 'select ko_accession, ko_id from enzyme_ko_annotation'
        ko_accession2ko_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
        locus2ko = parse_blast_koala_output(args.ko_table)
        locus2ko_table(locus2ko,
                       args.database_name,
                       ko_accession2ko_id)
