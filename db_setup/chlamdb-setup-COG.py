#!/usr/bin/env python


def create_COG_tables():
    import MySQLdb
    import os
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="COG") # name of the data base
    cursor = conn.cursor()

    sql1 = 'CREATE TABLE cog_2014 (domain_id varchar(100),' \
                            ' genome_name varchar(200),' \
                            ' protein_id INT,' \
                            ' protein_length INT,' \
                            ' domain_start INT,' \
                            ' domain_end INT,' \
                            ' COG_id INT,' \
                            ' membership_class INT,' \
                            ' index COG_id(COG_id),' \
                            ' index protein_id(protein_id))'

    sql2 = 'CREATE table cog_names_2014 (COG_id INT,' \
                                 ' COG_name varchar(100),' \
                                 ' description varchar(200),' \
                                 ' index COG_id(COG_id))'

    sql3 = 'CREATE table cog_id2cog_category (COG_id INT,' \
                                 ' category_id int,' \
                                 ' index COG_id(COG_id),' \
                                 ' index category_id(category_id))'


    sql4 = 'alter table COG.code2category add column `category_id` INT unsigned primary KEY AUTO_INCREMENT;'

    cursor.execute(sql1,)
    cursor.execute(sql2,)
    cursor.execute(sql3,)
    try:
        cursor.execute(sql4,)
    except:
        pass
    conn.commit()


def COG_id2cog_data_from_website(GOG_id):
    import urllib2
    from bs4 import BeautifulSoup
    module_file_file = 'https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid=%s' % GOG_id
    data = urllib2.urlopen(module_file_file)
    # COG0001	H	Glutamate-1-semialdehyde aminotransferase
    soup = BeautifulSoup(data, "lxml")
    html = soup.encode('utf-8')
    title = soup.findAll("div", { "class" : "desctit" })[0]
    description = title.text.split('[')[0]
    category = title.text.split('[')[1].split(']')[0]
    return [description, category]

def load_cog_tables(cognames_2014, cog_2014):
    '''

    COG names
    COG0001	H	Glutamate-1-semialdehyde aminotransferase

    COG
    158333741,Acaryochloris_marina_MBIC11017_uid58167,158333741,432,1,432,COG0001,0,
    <domain-id>, <genome-name>, <protein-id>,<protein-length>,
    <domain-start>, <domain-end>, <COG-id>, <membership-class>,

    :param cognames_2014:
    :param cog_2014:
    :return:
    '''
    import MySQLdb
    import os
    import re
    from biosqldb import manipulate_biosqldb
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="COG") # name of the data base
    cursor = conn.cursor()

    sql2 = 'select code, category_id from COG.code2category'
    cursor.execute(sql2,)
    cog_category2cog_category_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    sql3 = 'select description, category_id from COG.code2category'
    cursor.execute(sql3,)
    cog_description2cog_category_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    n=0
    with open(cognames_2014, 'r') as f:
        for row in f:
            data = row.rstrip().split('\t')
            if data[0][0] == '#':
                continue

            # split multiple function and insert one row/function
            fonctions = list(data[1])
            for fonction in fonctions:
                sql = 'insert into COG.cog_id2cog_category values (%s, %s)' % (n,
                                                                           cog_category2cog_category_id[fonction])
                cursor.execute(sql,)

            sql = 'insert into COG.cog_names_2014(COG_id, COG_name, description) values(%s,"%s","%s")' % (n,
                                                                                                          data[0],
                                                                                                          re.sub('"','',data[2]))
            cursor.execute(sql,)
            n+=1
    conn.commit()

    sql = 'select COG_name, COG_id from COG.cog_names_2014'

    cursor.execute(sql,)
    cog_name2cog_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    problem = []
    with open(cog_2014, 'r') as f:
        f = [i for i in f]
        for n, row in enumerate(f):
            if n%100 ==0:
                print ('%s / %s' % (n, len(f)))
            data = row.rstrip().split(',')
            if data[0][0] == '#':
                continue
            try:
                sql = 'insert into COG.cog_2014(COG_id, ' \
                      ' domain_id, ' \
                      ' genome_name, ' \
                      ' protein_id,' \
                      ' protein_length,' \
                      ' domain_start,' \
                      ' domain_end,' \
                      ' membership_class' \
                      ' ) values(%s,%s,"%s", %s, %s, %s, %s, %s)' % (cog_name2cog_id[data[6]],
                                                                         data[0],
                                                                         data[1],
                                                                         data[2],
                                                                         data[3],
                                                                         data[4],
                                                                         data[5],
                                                                         data[7])

                cursor.execute(sql,)
            except:
                # cog missing from the cognames table
                # fetch data from the NCBI website and insert new COG name entry
                # mail to Yuri Wolf, should not happen again (they cleaned tables from "old" removed COGs)
                print ('problem with COG id %s' % data[6])
                cog_web_data = COG_id2cog_data_from_website(data[6])
                description = cog_web_data[0]
                categories = cog_web_data[1]

                # get last COG id
                sql = 'select COG_id from cog_names_2014 order by COG_id DESC limit 1'
                cursor.execute(sql,)
                COG_id = int(cursor.fetchall()[0][0]) + 1
                cog_name2cog_id[data[6]] = COG_id
                for cat in cog_description2cog_category_id:
                    if cat in categories:
                        print ('category present!', cat)
                        category_id = cog_description2cog_category_id[cat]
                        sql = 'insert into COG.cog_id2cog_category values (%s, %s)' % (COG_id,
                                                                                       category_id)
                        cursor.execute(sql,)
                        conn.commit()

                # insert COG name
                sql = 'insert into COG.cog_names_2014(COG_id, COG_name, description) values(%s,"%s","%s")' % (cog_name2cog_id[data[6]],
                                                                                                              data[6],
                                                                                                              description)
                cursor.execute(sql,)
                conn.commit()

                # insert entry
                sql = 'insert into COG.cog_2014(COG_id, ' \
                      ' domain_id, ' \
                      ' genome_name, ' \
                      ' protein_id,' \
                      ' protein_length,' \
                      ' domain_start,' \
                      ' domain_end,' \
                      ' membership_class' \
                      ' ) values(%s,%s,"%s", %s, %s, %s, %s, %s)' % (COG_id,
                                                                     data[0],
                                                                     data[1],
                                                                     data[2],
                                                                     data[3],
                                                                     data[4],
                                                                     data[5],
                                                                     data[7])

                cursor.execute(sql,)
    conn.commit()
    print ('number of problems:', len(problem))


if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--cognames_2014', type=str, help="cognames2003-2014.tab")
    parser.add_argument("-c", '--cog_2014', type=str, help="cog2003-2014.csv")

    args = parser.parse_args()
    create_COG_tables()
    load_cog_tables(args.cognames_2014, args.cog_2014)
