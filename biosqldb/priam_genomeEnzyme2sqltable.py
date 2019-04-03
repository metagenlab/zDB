#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# create mysql 2 tables from dat file
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2016
# ---------------------------------------------------------------------------

def get_kegg_module_hierarchy():
    """
    get kegg module broad classification from we page
    :return: table
    """

    import urllib2
    import re

    url = "http://www.genome.jp/kegg-bin/download_htext?htext=ko00002.keg&format=htext&filedir="

    data = urllib2.urlopen(url)
    m = re.compile("A<b>([A-Za-z ]*)</b>.*")
    m2 = re.compile("B  <b>([A-Za-z ]*)</b>.*")
    m3 = re.compile("C    (.*)")
    main_category = re.compile("^A.*")
    sub_category = re.compile("^B.*")
    sub_sub_category = re.compile("^C.*")
    module = re.compile("D      (M.*)")
    module2category = {}
    for line in data:

         if re.match(main_category, line.strip()):

             main_cat = m.match(line).group(1)
         elif re.match(sub_category, line):
             try:
                sub_cat = m2.match(line).group(1)
             except:
                 continue
         elif re.match(sub_sub_category, line):
            sub_sub_cat = m3.match(line).group(1)
         elif re.match(module, line):
            #print main_cat, sub_cat, line.rstrip()

            mod = module.match(line).group(1).split(' ')[0]
            description = module.match(line).group(1).split('  ')[-1]
            print (description)
            #print main_cat, sub_cat, sub_sub_cat, mod
            module2category[mod] = [main_cat, sub_cat, sub_sub_cat, description]

         else:
             pass
    return module2category


def get_ko2ec(biodb):
    import urllib2
    import re
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'CREATE TABLE IF NOT EXISTS enzyme.ko2ec (ko_id VARCHAR(20),' \
           ' enzyme_id INT,' \
           ' index enzyme_id (enzyme_id), ' \
           ' index ko_id (ko_id));'
    server.adaptor.execute_and_fetchall(sql)


    url = 'http://rest.kegg.jp/list/ko'
    data = [line for line in urllib2.urlopen(url)]

    total = len(data)
    #data = ["ko:K00512"]
    for n, ko_data in enumerate(data):
        print ("%s / %s" % (n, total))

        ko_id = ko_data.split('\t')[0][3:]
        print (ko_id)

        url_ko = 'http://rest.kegg.jp/link/ec/%s' % ko_id

        ko_data = [line.rstrip().split('\t') for line in urllib2.urlopen(url_ko)]

        for one_ec in ko_data:
            try:
                ec_name = one_ec[1][3:]
            except:
                print ('problem with:', ko_id, 'no ec?')
                continue

            try:
                sql_ec_id = 'select enzyme_id from enzyme.enzymes where ec="%s";' % ec_name
                ec_id = server.adaptor.execute_and_fetchall(sql_ec_id,)[0][0]

                sql = 'INSERT INTO enzyme.ko2ec(ko_id, enzyme_id) VALUES ("%s",%s);' % (one_ec[0][3:], ec_id)
                #print sql
                server.adaptor.execute_and_fetchall(sql)

            except IndexError:
                print ('problem getting ec ID for:', ec_name, 'not in biosqldb?')
                print ('trying to add ec data from IUBMB')
                try:
                    ec_id = get_ec_data_from_IUBMB(ec_name)
                    print ("new_ec_id:", ec_id)
                    sql = 'INSERT INTO enzyme.ko2ec(ko_id, enzyme_id) VALUES ("%s",%s);' % (one_ec[0][3:], ec_id)
                    #print sql
                    server.adaptor.execute_and_fetchall(sql)
                except:
                    print ('FAIL')
        server.commit()




def get_complete_ko_table(biodb):
    import urllib2
    import re
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'CREATE TABLE IF NOT EXISTS enzyme.ko_annotation (ko_id INT AUTO_INCREMENT PRIMARY KEY, ' \
          ' ko_accession VARCHAR(20),' \
          ' name varchar(60),' \
          ' definition TEXT,' \
          ' EC TEXT,' \
          ' pathways TEXT,' \
          ' modules TEXT, ' \
          ' dbxrefs TEXT, ' \
          ' index ko_id (ko_id),' \
          ' index koa(ko_accession));'
    server.adaptor.execute_and_fetchall(sql)

    url = 'http://rest.kegg.jp/list/ko'
    data = [line for line in urllib2.urlopen(url)]

    sql = 'select ko_accession from enzyme.ko_annotation'
    ko_already_in_db = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    print ('test:', ko_already_in_db[0:10])

    total = len(data)
    for n, line in enumerate(data):
        print ("%s / %s" % (n, total))
        ko = line.rstrip().split('\t')[0][3:]
        print (ko)
        if ko in ko_already_in_db:
            print ('%s already in DB' % ko)
            continue
        url_ko = 'http://rest.kegg.jp/get/%s' % ko

        ko_data = [line for line in urllib2.urlopen(url_ko)]

        refs = ko2dbxrefs(ko_data)
        if refs:
            refs_str = ''
            for database in refs:
                if database != 'RN':
                    refs_str+='%s:%s,' % (database, refs[database])
            refs_str=refs_str[0:-1]
        else:
            refs_str = '-'
        definition = ko2definition(ko_data)

        if definition:
            definition = re.sub('"','', definition)
            if '[EC:' in definition:
                ec = definition.split('[EC:')[1][0:-1]
            else:
                ec = '-'
        else:
            definition = '-'
            ec = '-'
        modules = ko2MODULE(ko_data)
        if modules:
            modules = ','.join(modules)
        else:
            modules = '-'
        name = ko2name(ko_data)
        if name:
            name = ','.join(name)

        pathways = ko2PATHWAY(ko_data)
        if pathways:
            pathways = ','.join(pathways)
        else:
            pathways = '-'

        sql = 'insert into enzyme.ko_annotation (ko_accession, name, definition, EC, pathways, modules, dbxrefs)' \
              ' values ("%s", "%s", "%s","%s", "%s", "%s", "%s")' % (ko,
                                                                name,
                                                                definition,
                                                                ec,
                                                                pathways,
                                                                modules,
                                                                refs_str)

        print (sql)
        server.adaptor.execute(sql,)
        server.commit()



def ko2name(ko_record):
    import re
    m = re.compile("^NAME")


    for line in ko_record:
        if re.match(m, line):
            names=line.rstrip().split('        ')[1].split(',')
            return names

def ko2PATHWAY(ko_record):
    import re
    m = re.compile("^PATHWAY")
    m2 = re.compile("^ ")
    pathways =  []
    pathway=False
    for line in ko_record:
        if re.match(m, line):
            pw=line.rstrip().split('     ')[1].split('  ')[0]
            pathways.append(pw)
            pathway = True
        elif pathway and re.match(m2, line):
            pw=line.rstrip().split('            ')[1].split('  ')[0]
            pathways.append(pw)
        elif pathway and not re.match(m2, line):
            return pathways
        else:
            continue

def ko2dbxrefs(ko_record):
    import re
    m = re.compile("^DBLINKS")
    m2 = re.compile("^ ")
    refs =  {}
    ref=False
    for line in ko_record:
        if re.match(m, line):
            rr = line.rstrip().split('     ')[1].split(':')
            refs[rr[0]] = rr[1]
            ref = True
        elif ref and re.match(m2, line):
            rr = line.rstrip().split('            ')[1].split(':')
            refs[rr[0]] = rr[1]
        elif ref and not re.match(m2, line):
            return refs
        else:
            continue



def ko2MODULE(ko_record):
    import re
    m = re.compile("^MODULE")
    m2 = re.compile("^ ")
    modules =  []
    module=False
    for line in ko_record:
        if re.match(m, line):
            pw=line.rstrip().split('      ')[1].split('  ')[0]
            modules.append(pw)
            module = True
        elif module and re.match(m2, line):
            pw=line.rstrip().split('            ')[1].split('  ')[0]
            modules.append(pw)
        elif module and not re.match(m2, line):
            return modules
        else:
            continue


def ko2definition(ko_record):
    import re

    m = re.compile("^DEFINITION")

    for line in ko_record:
        if re.match(m, line):
            definition=line.rstrip().split('  ')[1]
            return definition


def locus2ko_table(locus_tag2ko_dico,
                   biodatabase,
                   ko_accession2ko_id):

    import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodatabase)

    sql2 = 'select locus_tag, seqfeature_id from annotation.seqfeature_id2locus_%s' % biodatabase

    locus2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2))

    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme.seqfeature_id2ko_%s (seqfeature_id INT,' \
           ' ko_id INT, ' \
           ' index ko_id (ko_id),' \
           ' index seqid (seqfeature_id));' % (biodatabase)

    server.adaptor.execute_and_fetchall(sql2,)

    for locus in locus_tag2ko_dico:
        ko = locus_tag2ko_dico[locus]
        ko_id = ko_accession2ko_id[ko]
        seqfeature_id = locus2seqfeature_id[locus]

        sql = 'insert into enzyme.seqfeature_id2ko_%s (seqfeature_id, ko_id) values (%s, %s)' % (biodatabase,
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


def definition2module_list(definition_string):
    import re
    return re.findall('M[0-9]+', definition_string)


def get_module2ko_data(cursor, ko_data, module_id, ko_accession2ko_id):
        import urllib2
        import re

        for line in ko_data:
            try:
                ko = line.rstrip().split("\t")[1][3:]
            except:
                print ('ko not found in line:', line)
                print ('ko data:', ko_data)
                import sys
                sys.exit()
            try:
                ko_url = "http://rest.kegg.jp/get/%s" % ko
                ko_data = urllib2.urlopen(ko_url)
                definition = ko2definition(ko_data)
                ko_id = ko_accession2ko_id[ko]
                sql = 'INSERT into module2ko (module_id, ko_id) values (%s,%s);' % (module_id,
                                                                                    ko_id)
                print (sql)
                cursor.execute(sql,)
            except:
                print (line)
                import sys
                sys.exit()


def get_module_table(module2category,ko_accession2ko_id):
    '''
    1. get all kegg pathways from API (http://rest.kegg.jp/) => create enzyme.kegg_module table
    2. get all KO associated for each module => create enzyme.module2ko table
    todo: remove existing tables for uptade if rerun

    :return: nothing
    '''


    import MySQLdb
    import urllib2
    import re
    import sys
    import os
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="enzyme") # name of the data base
    cursor = conn.cursor()

    sql1 = 'CREATE TABLE IF NOT EXISTS kegg_module (module_id INT AUTO_INCREMENT PRIMARY KEY,' \
           ' module_name VARCHAR(200),' \
           ' module_cat VARCHAR(200),' \
           ' module_sub_cat VARCHAR(200),' \
           ' module_sub_sub_cat VARCHAR(200),' \
           ' description LONG,' \
           ' index mdn(module_name),' \
           ' index mdc(module_cat),' \
           ' index mdsc(module_sub_cat),' \
           ' index mdssc(module_sub_sub_cat));'

    print (sql1)
    cursor.execute(sql1,)

    sql2 = 'CREATE TABLE IF NOT EXISTS module2ko (module_id INT,' \
           ' ko_id INT,' \
           ' index ko_id(ko_id),' \
           ' index module_id(module_id));'

    sqlm = 'select module_name from kegg_module'
    cursor.execute(sqlm,)
    module_in_db = [i[0] for i in cursor.fetchall()]


    sql3 = ' CONSTRAINT fk_pathway_id' \
           ' FOREIGN KEY(pathway_id) REFERENCES kegg_pathway(id)' \
           ' ON DELETE CASCADE,' \
           ' CONSTRAINT fk_ec_id' \
           ' FOREIGN KEY(ec_id) REFERENCES enzymes(id)' \
           ' ON DELETE CASCADE);'

    print (sql2)
    cursor.execute(sql2,)



    module_file_file = 'http://rest.kegg.jp/list/module'
    data = urllib2.urlopen(module_file_file)
    for line in data:
        pathway = line.rstrip().split("\t")
        module = pathway[0][3:]
        description = pathway[1]
        if module in module_in_db:
            print ('%s already in db' % module)
            continue
        try:
            cat = module2category[module][0]
            sub_cat = module2category[module][1]
            sub_sub_cat = module2category[module][2]
            description = module2category[module][3]
        except:
            cat = 'uncategorized'
            cat_short = 'uncategorized'
            print ('------------------------------------------------')

        sql = 'INSERT into kegg_module (module_name, module_cat,module_sub_cat, module_sub_sub_cat, description) ' \
              'values ("%s", "%s", "%s", "%s", "%s");' % (module,
                                                          cat,
                                                          sub_cat,
                                                          sub_sub_cat,
                                                          description)

        print (sql)
        cursor.execute(sql,)
        conn.commit()

        sql = 'SELECT LAST_INSERT_ID();'

        cursor.execute(sql, )

        try:
            module_id = cursor.fetchall()[0][0]
        except:
            pass

        ko_numbers_link = "http://rest.kegg.jp/link/ko/%s" % module
        ko_data = [l for l in urllib2.urlopen(ko_numbers_link)]

        print( len(ko_data))
        if ko_data[0] == '\n':
            print ('MODULE MADE OF SUBMODULES----------')
            module_link = "http://rest.kegg.jp/get/%s" % module
            module_data = [l for l in urllib2.urlopen(module_link)]
            definition = ko2definition(module_data)
            module_list = definition2module_list(definition)
            print ('modules:', module_list)
            for n,m in enumerate(module_list):
                print (n, m)
                ko_numbers_link2 = "http://rest.kegg.jp/link/ko/%s" % m
                ko_data2 = [l for l in urllib2.urlopen(ko_numbers_link2)]
                get_module2ko_data(cursor, ko_data2, module_id, ko_accession2ko_id)
        else:
            get_module2ko_data(cursor,ko_data,module_id, ko_accession2ko_id)

        conn.commit()


def get_pathway2ko(ko_accession2ko_id):
    '''
    1. get all kegg pathways from API (http://rest.kegg.jp/) => create enzyme.kegg_module table
    2. get all KO associated for each pathway => create enzyme.pathway2ko table
    todo: remove existing tables for uptade if rerun

    :return: nothing
    '''


    import MySQLdb
    import urllib2
    import re
    import sys
    import os
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="enzyme") # name of the data base
    cursor = conn.cursor()



    sql2 = 'CREATE TABLE IF NOT EXISTS pathway2ko (pathway_id INT,' \
           ' ko_id INT,' \
           ' index pathway_id(pathway_id),' \
           ' index ko_id(ko_id));'

    print (sql2)
    cursor.execute(sql2,)
    conn.commit()

    sql = 'select pathway_name, pathway_id from enzyme.kegg_pathway'
    pathway_name2pathway_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    pathway_file = 'http://rest.kegg.jp/list/pathway'
    data = urllib2.urlopen(pathway_file)
    for line in data:
        raw = line.rstrip().split("\t")
        pathway = raw[0][5:]
        description = raw[1]
        #print 'pathway', pathway
        #sql = 'select pathway_id from kegg_pathway where pathway_name="%s";' % pathway
        #print sql

        cursor.execute(sql)

        try:
            pathway_id = pathway_name2pathway_id[pathway] # cursor.fetchall()[0][0]
            print ('path id', pathway_id)
        except:
            print ('no pathway id for: %s, incomplete pathway table?' %  pathway)

        ko_numbers_link = "http://rest.kegg.jp/link/ko/%s" % pathway
        ko_data = urllib2.urlopen(ko_numbers_link)

        for line in ko_data:
            try:
                ko = line.rstrip().split("\t")[1][3:]

            except IndexError:
                print ('No Ko for pathway %s?' % pathway)
                continue
            try:
                ko_id = ko_accession2ko_id[ko]
            except KeyError:
                print ('PROBLEM with KO %s?' % ko)
                continue
            try:
                sql = 'INSERT into pathway2ko (pathway_id, ko_id) values (%s,%s);' % (pathway_id,
                                                                                       ko_id)
                cursor.execute(sql,)
            except:
                print (line)
                import sys
                sys.exit()
        conn.commit()

def get_ec_data_from_IUBMB(ec):
    import BeautifulSoup
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

    name_m = re.compile(u".*Accepted name.*")
    alname = re.compile(u".*Other name.*")
    reaction = re.compile(u".*Reaction:\<\/b\>.*")
    comments = re.compile(u".*Comments.*")

    conn.set_character_set('utf8')
    cursor.execute('SET NAMES utf8;')
    cursor.execute('SET CHARACTER SET utf8;')
    cursor.execute('SET character_set_connection=utf8;')

    ec_sep = ec.split('.')
    import requests

    adress = "http://www.chem.qmul.ac.uk/iubmb/enzyme/EC%s/%s/%s/%s.html" % (ec_sep[0], ec_sep[1], ec_sep[2], ec_sep[3])
    #print adress
    #data = requests.get(adress).text
    #print data
    request = urllib2.Request('http://www.somesite.com')
    html = urllib2.urlopen(adress).read()

    #html = urllib2.urlopen(adress).read()
    html = re.sub("\&\#","-",html)
    soup = BeautifulSoup(html, "lxml")
    html = soup.encode('utf-8')#.encode('latin-1') #encode('utf-8') # prettify()

    #print html
    #all_data = soup.find_all("p")#[i.get_text() for i in soup.find_all("p")]

    for i, data in enumerate(list(html.split('<p>'))):

        if re.match(name_m, data):
            print ('name')
            name = data.split("</b>")[-1]
            #name = re.sub("&alpha;","", name)
        elif re.match(alname, data):
            print ('altname')
            altname = data.split("):</b>")[1].split(';')
        elif re.match(reaction, data):
            print ('reaction')
            rr = data.split("<br>")
            reaction_list = []
            for i in rr:
                i = re.sub("\r\r<b>Reaction:</b> ","",i)
                i = re.sub("<[a-z\/]+>","",i)
                reaction_list.append(i)

        elif re.match(comments, data):
            print ('comment')
            cc = data.split("</b>")[1].split(';')
        else:
            continue




    sql_new = 'INSERT into enzymes (ec) values ("%s");' % ec
    cursor.execute(sql_new, )
    conn.commit()

    sql_id = 'SELECT LAST_INSERT_ID()'
    cursor.execute(sql_id, )
    id = int(cursor.fetchall()[0][0])
    print (id)

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
        #i = re.sub("\r\r<b>Reaction:</b> ","",i)
        sql = 'INSERT into enzymes_dat (enzyme_dat_id,line, value) values(%s, "comment", "%s");' % (id, re.sub("<[ =\"A-Za-z0-9\/\.]+>", "", i))
        cursor.execute(sql,)

    conn.commit()
    return id


def get_pathay_table(map2category):
    import MySQLdb
    import urllib2
    import sys
    import os
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="enzyme") # name of the data base
    cursor = conn.cursor()

    sql1 = 'CREATE TABLE IF NOT EXISTS enzyme.kegg_pathway (pathway_id INT AUTO_INCREMENT PRIMARY KEY,' \
           ' pathway_name VARCHAR(200),' \
           ' pathway_category_short VARCHAR(200),' \
           ' pathway_category VARCHAR(200),' \
           ' description LONG,' \
           ' index pn(pathway_name),' \
           ' index pcs(pathway_category_short),' \
           ' index pc(pathway_category));'

    print (sql1)
    cursor.execute(sql1,)

    conn.commit()

    pathway_file_file = 'http://rest.kegg.jp/list/pathway'
    data = urllib2.urlopen(pathway_file_file)

    print ('iter pathway list...')
    for line in data:
        pathway = line.rstrip().split("\t")
        print (pathway)
        map = pathway[0][5:]
        description = pathway[1]

        try:
            cat = map2category[map][1]
            cat_short = map2category[map][0]
        except:
            cat = 'uncategorized'
            cat_short = 'uncategorized'
        sql = 'INSERT into enzyme.kegg_pathway (pathway_name, description,' \
              ' pathway_category_short, pathway_category) values ("%s", "%s", "%s", "%s");' % (map,
                                                                                               description,
                                                                                               cat_short,
                                                                                               cat)

        cursor.execute(sql,)
        conn.commit()

def get_ec2get_pathway_table():
    '''
    1. get all kegg pathways from API (http://rest.kegg.jp/) => create enzyme.kegg_pathway table
    2. get all ec associated for each pathway => create enzyme.kegg2ec table
    todo: remove existing tables for uptade if rerun

    :return: nothing
    '''


    import MySQLdb
    import urllib2
    import sys
    import manipulate_biosqldb
    import os
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="enzyme") # name of the data base
    cursor = conn.cursor()

    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme.kegg2ec (pathway_id INT,' \
           ' ec_id INT,' \
           ' FOREIGN KEY(pathway_id) REFERENCES kegg_pathway(pathway_id),' \
           ' FOREIGN KEY(ec_id) REFERENCES enzymes(enzyme_id))' \


    sql3 = ' CONSTRAINT fk_pathway_id' \
           ' FOREIGN KEY(pathway_id) REFERENCES kegg_pathway(id)' \
           ' ON DELETE CASCADE,' \
           ' CONSTRAINT fk_ec_id' \
           ' FOREIGN KEY(ec_id) REFERENCES enzymes(id)' \
           ' ON DELETE CASCADE);'

    print (sql2)
    cursor.execute(sql2,)
    conn.commit()


    sql = 'select pathway_name,pathway_id from enzyme.kegg_pathway;'
    cursor.execute(sql,)
    map2id = manipulate_biosqldb.to_dict(cursor.fetchall())
    for map in map2id:
        id = map2id[map]

        ec_numbers_link = "http://rest.kegg.jp/link/ec/%s" % map
        print (ec_numbers_link)
        ec_data = urllib2.urlopen(ec_numbers_link)

        for line in ec_data:
            try:
                ec = line.rstrip().split("\t")[1]
            except:
                print ('problem:', line)
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
                    sys.stdout.write("trying to get enzyme data from IUBMB: %s...\n" % ec)
                    try:
                        ec_id = int(get_ec_data_from_IUBMB(ec))
                    except urllib2.HTTPError:
                        sys.stdout.write("NO DATA FOR %s\n" % ec)
                        continue

                sql = 'INSERT into kegg2ec (pathway_id, ec_id) values (%s,"%s");' % (id, ec_id)
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
                 try:
                    map_number = 'map'+re.findall("ko[0-9]+&", line)[0][2:-1]
                    map_description = re.findall(">(.*)<\/a", line)[0]
                    #print map_description
                    map2category[map_number] = [main_cat, sub_cat, map_description]
                 except:
                     try:
                        map_number = 'map'+re.findall("hsa[0-9]+&", line)[0][3:-1]
                        map_description = re.findall(">(.*)<\/a", line)[0]
                        #print map_description
                        map2category[map_number] = [main_cat, sub_cat, map_description]
                     except:
                         try:
                            map_number = 'map'+re.findall("[a-z]{3}[0-9]+&", line)[0][3:-1]
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

    print ('create enzyme table')
    print (sql1)
    cursor.execute(sql1,)
    print ('create dat table')
    print (sql2)
    cursor.execute(sql2)

    for data in Enzyme.parse(data):

        enzyme = data['ID']


        print (data)

        # insert enzyme id into primary TABLE
        sql = 'INSERT into enzymes (ec) values ("%s");' % enzyme

        print (sql)
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

    #sql = 'select locus_tag, accession from orthology_detail_%s' % biodatabase
    #sql2 = 'select locus_tag, orthogroup from orthology_detail_%s' % biodatabase
    #locus2bioentry_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql))
    #locus2orthogroup = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2))

    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme.seqfeature_id2ec_%s (enzyme_id INT AUTO_INCREMENT PRIMARY KEY,' \
           ' seqfeature_id INT,' \
           ' ec_id INT,' \
           ' FOREIGN KEY(ec_id) REFERENCES enzymes(enzyme_id));' % (biodatabase)

    server.adaptor.execute_and_fetchall(sql2,)

    sql = 'select locus_tag, seqfeature_id from annotation.seqfeature_id2locus_%s' % biodatabase

    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    for locus in locus_tag2ec_dico:
        for ec_data in locus_tag2ec_dico[locus]:

            sql = 'select enzyme_id from enzyme.enzymes where ec="%s"' % ec_data[0]
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
    import manipulate_biosqldb
    import parse_priam_EC

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_priam_files', type=str, help="input interpro csv file", nargs='+', default=False)
    parser.add_argument("-k", '--ko_table', type=str, help="input blastGhost file", default=False)
    parser.add_argument("-d", '--database_name', type=str, help="database name")
    parser.add_argument("-u", '--update_database', action="store_true", help="update enzyme and KEGG data from online Databases")


    args = parser.parse_args()



    if args.update_database:

        server, db = manipulate_biosqldb.load_db(args.database_name)

        #get_complete_ko_table(args.database_name)
        #load_enzyme_nomenclature_table()
        sql = 'select ko_accession, ko_id from enzyme.ko_annotation'
        ko_accession2ko_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
        #print 'getting map2category...'
        #map2category = get_kegg_pathway_classification()
        #print 'getting pathway table...'
        #get_pathay_table(map2category)
        print ('getting ko2pathway...')
        get_pathway2ko(ko_accession2ko_id)
        print ('getting module2ko...')
        get_module_table(get_kegg_module_hierarchy(), ko_accession2ko_id)

        #get_ec2get_pathway_table()
        #get_ec_data_from_IUBMB("1.14.13.217")
        #get_ec_data_from_IUBMB("1.1.1.1")

        #get_microbial_metabolism_in_diverse_environments_kegg01120()
        #get_ko2ec(args.database_name)
        #get_ec_data_from_IUBMB("3.2.1.196")
    if args.input_priam_files:
        locus2ec={}

        for priam_file in args.input_priam_files:
            locus2ec.update(parse_priam_EC.locus2EC(priam_file))

        locus2ec_table(locus2ec, args.database_name)


    if args.ko_table:
        server, db = manipulate_biosqldb.load_db(args.database_name)
        sql = 'select ko_accession, ko_id from enzyme.ko_annotation'
        ko_accession2ko_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
        locus2ko = parse_blast_koala_output(args.ko_table)
        locus2ko_table(locus2ko,
                       args.database_name,
                       ko_accession2ko_id)
