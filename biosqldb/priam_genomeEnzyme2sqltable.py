#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# create mysql 2 tables from dat file
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2016
# ---------------------------------------------------------------------------


def get_ec_data_from_IUBMB(ec):
    import BeautifulSoup
    import urllib2
    import re
    import MySQLdb

    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd="estrella3", # your password
                                db="enzyme") # name of the data base
    cursor = conn.cursor()

    name = re.compile(".*Accepted name.*")
    alname = re.compile(".*Other name.*")
    reaction = re.compile(".*Reaction:\<\/b\>.*")
    comments = re.compile(".*Comments.*")


    ec_sep = ec.split('.')
    adress = "http://www.chem.qmul.ac.uk/iubmb/enzyme/EC%s/%s/%s/%s.html" % (ec_sep[0], ec_sep[1], ec_sep[2], ec_sep[3])
    #print adress
    html = urllib2.urlopen(adress).read()
    for i, data in enumerate(list(html.split('<p>'))):
        if re.match(name, data):
            name = data.split("</b>")[-1]
        elif re.match(alname, data):
            altname = data.split("):</b>")[1].split(';')
        elif re.match(reaction, data):
            rr = data.split("<br>")
            reaction_list = []
            for i in rr:
                i = re.sub("\r\r<b>Reaction:</b> ","",i)
                i = re.sub("<[a-z\/]+>","",i)
                reaction_list.append(i)

        elif re.match(comments, data):
            cc = data.split("</b>")[1].split(';')
        else:
            continue




    sql_new = 'INSERT into enzymes (ec) values ("%s");' % ec
    cursor.execute(sql_new, )
    conn.commit()

    sql_id = 'SELECT LAST_INSERT_ID()'
    cursor.execute(sql_id, )
    id = int(cursor.fetchall()[0][0])
    print id

    sql = 'INSERT into enzymes_dat (enzyme_dat_id, line, value) values (%s, "description", "%s");' % (id, name)
    cursor.execute(sql,)
    # alternative names
    for i in altname:
        sql = 'INSERT into enzymes_dat (enzyme_dat_id,line, value) values(%s, "alternative name", "%s");' % (id, re.sub("<[a-z\/]+>", "", i))
        cursor.execute(sql,)

    # Catalytic activity
    for i in reaction_list:
        sql = 'INSERT into enzymes_dat (enzyme_dat_id,line, value) values(%s,"catalytic activity", "%s");' % (id, re.sub("<[a-z\/]+>", "", i))
        cursor.execute(sql,)

    # comments
    for i in cc:
        sql = 'INSERT into enzymes_dat (enzyme_dat_id,line, value) values(%s, "comment", "%s");' % (id, re.sub("<[a-z\/]+>", "", i))
        cursor.execute(sql,)

    conn.commit()
    return id




def get_ec2get_pathway_table(map2category):
    '''
    1. get all kegg pathways from API (http://rest.kegg.jp/) => create enzyme.kegg_pathway table
    2. get all ec associated for each pathway => create enzyme.kegg2ec table
    todo: remove existing tables for uptade if rerun

    :return: nothing
    '''


    import MySQLdb
    import urllib2
    import sys

    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd="estrella3", # your password
                                db="enzyme") # name of the data base
    cursor = conn.cursor()

    sql1 = 'CREATE TABLE IF NOT EXISTS kegg_pathway (pathway_id INT AUTO_INCREMENT PRIMARY KEY,' \
           ' pathway_name VARCHAR(200),' \
           ' pathway_category_short VARCHAR(200),' \
           ' pathway_category VARCHAR(200),' \
           ' description LONG);'

    print sql1
    cursor.execute(sql1,)

    sql2 = 'CREATE TABLE IF NOT EXISTS kegg2ec (pathway_id INT,' \
           ' ec_id INT,' \
           ' FOREIGN KEY(pathway_id) REFERENCES kegg_pathway(pathway_id),' \
           ' FOREIGN KEY(ec_id) REFERENCES enzymes(enzyme_id))' \


    sql3 = ' CONSTRAINT fk_pathway_id' \
           ' FOREIGN KEY(pathway_id) REFERENCES kegg_pathway(id)' \
           ' ON DELETE CASCADE,' \
           ' CONSTRAINT fk_ec_id' \
           ' FOREIGN KEY(ec_id) REFERENCES enzymes(id)' \
           ' ON DELETE CASCADE);'

    print sql2
    cursor.execute(sql2,)



    pathway_file_file = 'http://rest.kegg.jp/list/pathway'
    data = urllib2.urlopen(pathway_file_file)
    for line in data:
        pathway = line.rstrip().split("\t")
        map = pathway[0][5:]
        description = pathway[1]

        try:
            cat = map2category[map][1]
            cat_short = map2category[map][0]
        except:
            cat = 'uncategorized'
            cat_short = 'uncategorized'
        sql = 'INSERT into kegg_pathway (pathway_name, description,pathway_category_short, pathway_category) values ("%s", "%s", "%s", "%s");' % (map,
                                                                                                                                                  description,
                                                                                                                                                  cat_short,
                                                                                                                                                  cat)

        print sql
        cursor.execute(sql,)
        conn.commit()

        sql = 'SELECT LAST_INSERT_ID();'

        cursor.execute(sql, )
        try:
            id = cursor.fetchall()[0][0]
        except:
            pass

        #print 'id', id


        ec_numbers_link = "http://rest.kegg.jp/link/ec/%s" % map
        ec_data = urllib2.urlopen(ec_numbers_link)

        for line in ec_data:
            try:
                ec = line.rstrip().split("\t")[1]
            except:
                print 'problem:', line
                continue
            if '-' in ec:
                continue
            else:
                ec = ec[3:]
                sql_ec = 'select enzyme_id from enzymes where ec = "%s"' % ec
                #print sql_ec
                cursor.execute(sql_ec, )
                try:
                    ec_id = cursor.fetchall()[0][0]
                except:
                    sys.stdout.write("trying to get enzyme data from IUBMB...")
                    try:
                        ec_id = int(get_ec_data_from_IUBMB(ec))
                    except urllib2.HTTPError:
                        sys.stdout.write("NO DATA FOR %s" % ec)
                        continue

                sql = 'INSERT into kegg2ec (pathway_id, ec_id) values (%s,"%s");' % (id, ec_id)
                #print sql
                cursor.execute(sql,)
        conn.commit()

def get_kegg_pathway_classification():
    """
    get kegg pathway bread classification from we page
    :return: table
    """

    import urllib2
    import re

    url = "http://www.genome.jp/kegg/pathway.html"

    data = urllib2.urlopen(url)

    main_category = re.compile("<a name=.*")
    sub_category = re.compile("<b>.*")
    map = re.compile(".*<a href.*")
    map2category = {}
    for line in data:
         if re.match(main_category, line.strip()):
             main_cat = line.rstrip().split("\"")[1]
         elif re.match(sub_category, line):
             try:
                sub_cat = line.rstrip().split("<b>")[1][0:-4]
             except:
                continue
         elif re.match(map, line):
             #print main_cat, sub_cat, line.rstrip()
             try:
                 #print main_cat, sub_cat, line.rstrip()
                 map_number = re.findall("map[0-9]+", line)[0]
                 map_description = re.findall(">(.*)<\/a", line)[0]
                 #print map_description
                 map2category[map_number] = [main_cat, sub_cat, map_description]
             except:
                 pass
         else:
             pass
    return map2category


def load_enzyme_nomenclature_table():


    '''

    download all SIB enzyme nomenclature from FTP (ftp://ftp.expasy.org/databases/enzyme/)
    create the enzyme.enzymes table with the list of all EC with associated description
    create the enzyme.enzymes_dat with detailed information about each EC
    todo: remove existing tables for uptade if rerun

    :return: nothing
    '''

    from Bio.ExPASy import Enzyme
    import MySQLdb
    import urllib2

    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd="estrella3", # your password
                                db="enzyme") # name of the data base
    cursor = conn.cursor()

    '''
    ID  Identification                         (Begins each entry; 1 per entry)
    DE  Description (official name)            (>=1 per entry)
    AN  Alternate name(s)                      (>=0 per entry)
    CA  Catalytic activity                     (>=1 per entry)
    CF  Cofactor(s)                            (>=0 per entry)
    CC  Comments                               (>=0 per entry)
    PR  Cross-references to PROSITE            (>=0 per entry)
    DR  Cross-references to Swiss-Prot         (>=0 per entry)
    '''



    enzyme_file = 'ftp://ftp.expasy.org/databases/enzyme/enzyme.dat'

    data = urllib2.urlopen(enzyme_file)


    sql1 = 'CREATE TABLE IF NOT EXISTS enzymes (enzyme_id INT AUTO_INCREMENT PRIMARY KEY,' \
          ' ec VARCHAR(200));'

    sql2 = 'CREATE TABLE IF NOT EXISTS enzymes_dat (enzyme_dat_id INT,' \
          ' line VARCHAR(20),' \
          ' value LONG,' \
           ' CONSTRAINT fk_enzyme_id' \
           ' FOREIGN KEY(enzyme_dat_id) REFERENCES enzymes(enzyme_id)' \
           ' ON DELETE CASCADE);'

    print 'create enzyme table'
    print sql1
    cursor.execute(sql1,)
    print 'create dat table'
    print sql2
    cursor.execute(sql2)

    for data in Enzyme.parse(data):

        enzyme = data['ID']


        print  data

        # insert enzyme id into primary TABLE
        sql = 'INSERT into enzymes (ec) values ("%s");' % enzyme

        print sql
        cursor.execute(sql,)
        conn.commit()

        sql = 'SELECT LAST_INSERT_ID();'

        cursor.execute(sql, )
        id = cursor.fetchall()[0][0]

        #print id

        # description
        sql = 'INSERT into enzymes_dat (enzyme_dat_id, line, value) values (%s, "description", "%s");' % (id, data['DE'])
        cursor.execute(sql,)
        # alternative names
        for i in data['AN']:
            sql = 'INSERT into enzymes_dat (enzyme_dat_id,line, value) values(%s, "alternative name", "%s");' % (id, i)
            cursor.execute(sql,)

        # Catalytic activity
        sql = 'INSERT into enzymes_dat (enzyme_dat_id,line, value) values(%s,"catalytic activity", "%s");' % (id, data['CA'])
        cursor.execute(sql,)

        # Cofactors
        sql = 'INSERT into enzymes_dat (enzyme_dat_id,line, value) values(%s, "cofactors", "%s");' % (id, data['CF'])
        cursor.execute(sql,)

        # prosite crossref
        for i in data['PR']:
            sql = 'INSERT into enzymes_dat (enzyme_dat_id,line, value) values(%s, "prosite", "%s");' % (id, i)
            cursor.execute(sql,)
        # comments
        for i in data['CC']:
            sql = 'INSERT into enzymes_dat (enzyme_dat_id,line, value) values(%s, "comment", "%s");' % (id, i)
            cursor.execute(sql,)

        conn.commit()


def locus2ec_table(locus_tag2ec_dico, biodatabase):

    import manipulate_biosqldb
    import biosql_own_sql_tables

    server, db = manipulate_biosqldb.load_db(biodatabase)

    sql = 'select locus_tag, accession from orthology_detail_%s' % biodatabase
    sql2 = 'select locus_tag, orthogroup from orthology_detail_%s' % biodatabase

    locus2accession = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql))
    locus2orthogroup = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2))

    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme.locus2ec_%s (enzyme_id INT AUTO_INCREMENT PRIMARY KEY,' \
           ' accession varchar(60),'\
           ' locus_tag VARCHAR(20),' \
           ' orthogroup varchar(20),' \
           ' ec_id INT,' \
           ' FOREIGN KEY(ec_id) REFERENCES enzymes(enzyme_id));' % (biodatabase)

    print sql2
    server.adaptor.execute_and_fetchall(sql2,)

    for locus in locus_tag2ec_dico:
        for ec_data in locus_tag2ec_dico[locus]:
            print locus, ec_data[0], locus2accession[locus], locus2orthogroup[locus]

            sql = 'select enzyme_id from enzyme.enzymes where ec="%s"' % ec_data[0]
            ec_id = server.adaptor.execute_and_fetchall(sql,)[0][0]

            sql = 'insert into enzyme.locus2ec_%s (accession, locus_tag, orthogroup, ec_id) values ("%s", "%s", "%s", %s)' % (biodatabase,
                                                                                                                       locus2accession[locus],
                                                                                                                       locus,
                                                                                                                       locus2orthogroup[locus],
                                                                                                                       ec_id)
            print sql
            server.adaptor.execute(sql,)
            server.commit()

if __name__ == '__main__':
    import argparse
    import manipulate_biosqldb
    import parse_priam_EC

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_priam_files', type=str, help="input interpro csv file", nargs='+', default=False)
    parser.add_argument("-d", '--database_name', type=str, help="database name")
    parser.add_argument("-u", '--update_database', action="store_true", help="update enzyme and KEGG data from online Databases")


    args = parser.parse_args()

    if args.update_database:
        load_enzyme_nomenclature_table()
        get_ec2get_pathway_table(get_kegg_pathway_classification())
        #get_ec_data_from_IUBMB("1.14.13.217")
        #get_ec_data_from_IUBMB("1.1.1.1")

    if args.input_priam_files:
        locus2ec={}

        for priam_file in args.input_priam_files:
            locus2ec.update(parse_priam_EC.locus2EC(priam_file))

        locus2ec_table(locus2ec, args.database_name)

    # get_ec2get_pathway_table(get_kegg_pathway_classification())
