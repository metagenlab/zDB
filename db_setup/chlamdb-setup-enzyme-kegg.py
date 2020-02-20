#!/usr/bin/env python


def connect_db(biodb):
    import MySQLdb
    import os
    import re
    from chlamdb.biosqldb import manipulate_biosqldb
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                           user="root", # your username
                           passwd=sqlpsw,
                           db=biodb) # name of the data base

    cursor = conn.cursor()

    conn.set_character_set('utf8')
    cursor.execute('SET NAMES utf8;')
    cursor.execute('SET CHARACTER SET utf8;')
    cursor.execute('SET character_set_connection=utf8;')

    return conn, cursor


def get_complete_ko_table(biodb):
    import urllib.request
    import re
    conn, cursor = connect_db(biodb)

    sql = 'CREATE TABLE IF NOT EXISTS enzyme_ko_annotation (ko_id INT AUTO_INCREMENT PRIMARY KEY, ' \
          ' ko_accession VARCHAR(20),' \
          ' name varchar(60),' \
          ' definition TEXT,' \
          ' EC TEXT,' \
          ' pathways TEXT,' \
          ' modules TEXT, ' \
          ' dbxrefs TEXT, ' \
          ' index ko_id (ko_id),' \
          ' index koa(ko_accession));'
    cursor.execute(sql)

    url = 'http://rest.kegg.jp/list/ko'
    data = urllib.request.urlopen(url).read().decode('utf-8').split('\n')

    sql = 'select ko_accession from enzyme_ko_annotation;'
    cursor.execute(sql,)

    ko_already_in_db = [i[0] for i in cursor.fetchall()]

    print('Already into DB:', ko_already_in_db[0:10])

    total = len(data)
    for n, line in enumerate(data):
        print("%s / %s" % (n, total))
        # manage empty row(s)
        if len(line) == 0:
            continue
        ko = line.rstrip().split('\t')[0][3:]
        print (ko)
        if ko in ko_already_in_db:
            print ('%s already in DB' % ko)
            continue
        url_ko = 'http://rest.kegg.jp/get/%s' % ko

        ko_data = urllib.request.urlopen(url_ko).read().decode('utf-8').split("\n")

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

        sql = 'insert into enzyme_ko_annotation (ko_accession, name, definition, EC, pathways, modules, dbxrefs)' \
              ' values ("%s", "%s", "%s","%s", "%s", "%s", "%s")' % (ko,
                                                                name,
                                                                definition,
                                                                ec,
                                                                pathways,
                                                                modules,
                                                                refs_str)

        cursor.execute(sql,)
        conn.commit()


def load_enzyme_nomenclature_table(biodb):


    '''

    download all SIB enzyme nomenclature from FTP (ftp://ftp.expasy.org/databases/enzyme/)
    create the enzyme.enzymes table with the list of all EC with associated description
    create the enzyme.enzymes_dat with detailed information about each EC
    todo: remove existing tables for uptade if rerun

    :return: nothing
    '''

    from Bio.ExPASy import Enzyme
    import MySQLdb
    import urllib.request
    import os
    from io import StringIO
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db=biodb) # name of the data base
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

    data = urllib.request.urlopen(enzyme_file).read().decode('utf-8')

    sql1 = 'CREATE TABLE IF NOT EXISTS enzyme_enzymes (enzyme_id INT AUTO_INCREMENT PRIMARY KEY,' \
          ' ec VARCHAR(200));'

    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme_enzymes_dat (enzyme_dat_id INT,' \
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

    for n, data in enumerate(Enzyme.parse(StringIO(data))):
        enzyme = data['ID']
        # insert enzyme id into primary TABLE
        sql = 'INSERT into enzyme_enzymes (ec) values ("%s");' % enzyme

        print(n, sql)

        cursor.execute(sql,)
        conn.commit()

        sql = 'SELECT LAST_INSERT_ID();'

        cursor.execute(sql, )
        id = cursor.fetchall()[0][0]

        # description
        sql = 'INSERT into enzyme_enzymes_dat (enzyme_dat_id, line, value) values (%s, "description", "%s");' % (id, data['DE'])
        cursor.execute(sql,)
        # alternative names
        for i in data['AN']:
            sql = 'INSERT into enzyme_enzymes_dat (enzyme_dat_id,line, value) values(%s, "alternative name", "%s");' % (id, i)
            cursor.execute(sql,)

        # Catalytic activity
        sql = 'INSERT into enzymes_dat (enzyme_dat_id,line, value) values(%s,"catalytic activity", "%s");' % (id, data['CA'])
        cursor.execute(sql,)

        # Cofactors
        sql = 'INSERT into enzyme_enzymes_dat (enzyme_dat_id,line, value) values(%s, "cofactors", "%s");' % (id, data['CF'])
        cursor.execute(sql,)

        # prosite crossref
        for i in data['PR']:
            sql = 'INSERT into enzyme_enzymes_dat (enzyme_dat_id,line, value) values(%s, "prosite", "%s");' % (id, i)
            cursor.execute(sql,)
        # comments
        for i in data['CC']:
            sql = 'INSERT into enzyme_enzymes_dat (enzyme_dat_id,line, value) values(%s, "comment", "%s");' % (id, i)
            cursor.execute(sql,)

        conn.commit()



def get_kegg_pathway_classification():
    """
    get kegg pathway bread classification from we page
    :return: table
    """

    import urllib.request
    import re

    url = "http://www.genome.jp/kegg/pathway.html"

    data = urllib.request.urlopen(url).read().decode('utf-8').split("\n")

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



def get_pathay_table(map2category, 
                     biodb):
    import MySQLdb
    import urllib.request
    import sys
    import os
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db=biodb) # name of the data base
    cursor = conn.cursor()

    sql1 = 'CREATE TABLE IF NOT EXISTS enzyme_kegg_pathway (pathway_id INT AUTO_INCREMENT PRIMARY KEY,' \
           ' pathway_name VARCHAR(200),' \
           ' pathway_category_short VARCHAR(200),' \
           ' pathway_category VARCHAR(200),' \
           ' description LONG,' \
           ' index pn(pathway_name),' \
           ' index pcs(pathway_category_short),' \
           ' index pc(pathway_category));'

    cursor.execute(sql1,)
    conn.commit()

    pathway_file_file = 'http://rest.kegg.jp/list/pathway'
    data = urllib.request.urlopen(pathway_file_file).read().decode('utf-8').split("\n")

    print ('iter pathway list...')
    for line in data:
        # manage empty line(s)
        if len(line) == 0:
            continue
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
        sql = 'INSERT into enzyme_kegg_pathway (pathway_name, description,' \
              ' pathway_category_short, pathway_category) values ("%s", "%s", "%s", "%s");' % (map,
                                                                                               description,
                                                                                               cat_short,
                                                                                               cat)

        cursor.execute(sql,)
        conn.commit()


def get_pathway2ko(ko_accession2ko_id, 
                   biodb):
    '''
    1. get all kegg pathways from API (http://rest.kegg.jp/) => create enzyme.kegg_module table
    2. get all KO associated for each pathway => create enzyme_pathway2ko table
    todo: remove existing tables for uptade if rerun

    :return: nothing
    '''


    import MySQLdb
    import urllib.request
    import re
    import sys
    import os
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db=biodb) # name of the data base
    cursor = conn.cursor()



    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme_pathway2ko (pathway_id INT,' \
           ' ko_id INT,' \
           ' index pathway_id(pathway_id),' \
           ' index ko_id(ko_id));'

    print (sql2)
    cursor.execute(sql2,)
    conn.commit()

    sql = 'select pathway_name, pathway_id from enzyme_kegg_pathway'
    cursor.execute(sql,)
    pathway_name2pathway_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    pathway_file = 'http://rest.kegg.jp/list/pathway'
    data = urllib.request.urlopen(pathway_file).read().decode('utf-8').split("\n")
    for line in data:
        # handle empty line(s)
        if len(line) == 0:
            continue
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
        ko_data = urllib.request.urlopen(ko_numbers_link).read().decode('utf-8').split("\n")

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
                sql = 'INSERT into enzyme_pathway2ko (pathway_id, ko_id) values (%s,%s);' % (pathway_id,
                                                                                       ko_id)
                cursor.execute(sql,)
            except:
                print (line)
                import sys
                sys.exit()
        conn.commit()


def get_module_table(module2category,
                     ko_accession2ko_id,
                     biodb):
    '''
    1. get all kegg pathways from API (http://rest.kegg.jp/) => create enzyme.kegg_module table
    2. get all KO associated for each module => create enzyme.module2ko table
    todo: remove existing tables for uptade if rerun

    :return: nothing
    '''


    import MySQLdb
    import urllib.request
    import re
    import sys
    import os
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db=biodb) # name of the data base
    cursor = conn.cursor()

    sql1 = 'CREATE TABLE IF NOT EXISTS enzyme_kegg_module (module_id INT AUTO_INCREMENT PRIMARY KEY,' \
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

    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme_module2ko (module_id INT,' \
           ' ko_id INT,' \
           ' index ko_id(ko_id),' \
           ' index module_id(module_id));'

    sqlm = 'select module_name from enzyme_kegg_module'
    cursor.execute(sqlm,)
    module_in_db = [i[0] for i in cursor.fetchall()]


    sql3 = ' CONSTRAINT fk_pathway_id' \
           ' FOREIGN KEY(pathway_id) REFERENCES enzyme_kegg_pathway(id)' \
           ' ON DELETE CASCADE,' \
           ' CONSTRAINT fk_ec_id' \
           ' FOREIGN KEY(ec_id) REFERENCES enzyme_enzymes(id)' \
           ' ON DELETE CASCADE);'

    print (sql2)
    cursor.execute(sql2,)

    module_file_file = 'http://rest.kegg.jp/list/module'
    data = urllib.request.urlopen(module_file_file).read().decode('utf-8').split("\n")
    for line in data:
        if len(line) == 0:
            continue
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

        sql = 'INSERT into enzyme_kegg_module (module_name, module_cat,module_sub_cat, module_sub_sub_cat, description) ' \
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
        ko_data = urllib.request.urlopen(ko_numbers_link).read().decode('utf-8').split("\n")

        if ko_data[0] == '\n':
            print ('MODULE MADE OF SUBMODULES----------')
            module_link = "http://rest.kegg.jp/get/%s" % module
            module_data = urllib.request.urlopen(module_link).read().decode('utf-8').split("\n")
            definition = ko2definition(module_data)
            module_list = definition2module_list(definition)
            print ('modules:', module_list)
            for n, m in enumerate(module_list):
                print (n, m)
                ko_numbers_link2 = "http://rest.kegg.jp/link/ko/%s" % m
                ko_data2 = urllib.request.urlopen(ko_numbers_link2).read().decode('utf-8').split("\n")
                get_module2ko_data(cursor, ko_data2, module_id, ko_accession2ko_id)
        else:
            get_module2ko_data(cursor,
                               ko_data,
                               module_id,
                               ko_accession2ko_id)

        conn.commit()



def get_kegg_module_hierarchy():
    """
    get kegg module broad classification from we page
    :return: table
    """
    import urllib.request
    import re
    url = "http://www.genome.jp/kegg-bin/download_htext?htext=ko00002.keg&format=htext&filedir="
    data = urllib.request.urlopen(url).read().decode('utf-8').split("\n")
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
            #print main_cat, sub_cat, sub_sub_cat, mod
            module2category[mod] = [main_cat, sub_cat, sub_sub_cat, description]
         else:
             pass
    return module2category


def get_ec2get_pathway_table(biodb):
    '''
    1. get all kegg pathways from API (http://rest.kegg.jp/) => create enzyme.kegg_pathway table
    2. get all ec associated for each pathway => create enzyme.kegg2ec table
    todo: remove existing tables for uptade if rerun

    :return: nothing
    '''


    import MySQLdb
    import urllib.request
    import sys
    from chlamdb.biosqldb import manipulate_biosqldb
    import os
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db=biodb) # name of the data base
    cursor = conn.cursor()

    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme_kegg2ec (pathway_id INT,' \
           ' ec_id INT,' \
           ' FOREIGN KEY(pathway_id) REFERENCES enzyme_kegg_pathway(pathway_id),' \
           ' FOREIGN KEY(ec_id) REFERENCES enzyme_enzymes(enzyme_id))' \


    sql3 = ' CONSTRAINT fk_pathway_id' \
           ' FOREIGN KEY(pathway_id) REFERENCES enzyme_kegg_pathway(id)' \
           ' ON DELETE CASCADE,' \
           ' CONSTRAINT fk_ec_id' \
           ' FOREIGN KEY(ec_id) REFERENCES enzyme_enzymes(id)' \
           ' ON DELETE CASCADE);'

    print (sql2)
    cursor.execute(sql2,)
    conn.commit()


    sql = 'select pathway_name,pathway_id from enzyme_kegg_pathway;'
    cursor.execute(sql,)
    map2id = manipulate_biosqldb.to_dict(cursor.fetchall())
    for map in map2id:
        id = map2id[map]

        ec_numbers_link = "http://rest.kegg.jp/link/ec/%s" % map
        print (ec_numbers_link)
        ec_data = urllib.request.urlopen(ec_numbers_link).read().decode('utf-8').split("\n")

        for line in ec_data:
            if len(line) == 0:
                continue
            try:
                ec = line.rstrip().split("\t")[1]
            except:
                print ('problem:', line)
                continue
            if '-' in ec:
                continue
            else:
                ec = ec[3:]
                sql_ec = 'select enzyme_id from enzyme_enzymes where ec="%s"' % ec
                cursor.execute(sql_ec, )
                try:
                    ec_id = cursor.fetchall()[0][0]
                except IndexError:
                    print("problem:", sql_ec)
                    sys.stdout.write("trying to get enzyme data from IUBMB: %s...\n" % ec)
                    try:
                        ec_id = int(get_ec_data_from_IUBMB(ec))
                    except urllib.request.HTTPError:
                        sys.stdout.write("NO DATA FOR %s\n" % ec)
                        continue

                sql = 'INSERT into enzyme_kegg2ec (pathway_id, ec_id) values (%s,"%s");' % (id, ec_id)
                cursor.execute(sql,)
        conn.commit()


def definition2module_list(definition_string):
    import re
    return re.findall('M[0-9]+', definition_string)


def get_module2ko_data(cursor, ko_data, module_id, ko_accession2ko_id):
        import urllib.request
        import re

        for line in ko_data:
            if len(line) == 0:
                continue
            try:
                ko = line.rstrip().split("\t")[1][3:]
            except:
                print ('ko not found in line:', line)
                print ('ko data:', ko_data)
                import sys
                sys.exit()
            try:
                ko_url = "http://rest.kegg.jp/get/%s" % ko
                ko_data = urllib.request.urlopen(ko_url).read().decode('utf-8').split("\n")
                definition = ko2definition(ko_data)
                ko_id = ko_accession2ko_id[ko]
                sql = 'INSERT into enzyme_module2ko (module_id, ko_id) values (%s,%s);' % (module_id,
                                                                                    ko_id)
                print (sql)
                cursor.execute(sql,)
            except:
                print (line)
                import sys
                sys.exit()


def get_ec_data_from_IUBMB(ec, 
                           biodb):
    import urllib.request
    import re
    import MySQLdb
    from bs4 import BeautifulSoup
    import os
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost",
                                user="root",
                                passwd=sqlpsw,
                                db=biodb)
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

    adress = "https://www.qmul.ac.uk/sbcs/iubmb/enzyme/EC%s/%s/%s/%s.html" % (ec_sep[0], ec_sep[1], ec_sep[2], ec_sep[3])

    html = urllib.request.urlopen(adress).read().decode('utf-8')

    #html = urllib2.urlopen(adress).read()
    html = re.sub("\&\#", "-", html)
    soup = BeautifulSoup(html, "lxml")
    html = str(soup.encode('utf-8'))#.encode('latin-1') #encode('utf-8') # prettify()

    print(html)
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

    sql_new = 'INSERT into enzyme_enzymes (ec) values ("%s");' % ec
    cursor.execute(sql_new, )
    conn.commit()

    sql_id = 'SELECT LAST_INSERT_ID()'
    cursor.execute(sql_id, )
    id = int(cursor.fetchall()[0][0])
    print (id)

    sql = 'INSERT into enzyme_enzymes_dat (enzyme_dat_id, line, value) values (%s, "description", "%s");' % (id, name)
    cursor.execute(sql,)
    # alternative names
    for i in altname:
        sql = 'INSERT into enzyme_enzymes_dat (enzyme_dat_id,line, value) values(%s, "alternative name", "%s");' % (id, re.sub("<[a-z\/]+>", "", i))
        cursor.execute(sql,)

    # Catalytic activity
    for i in reaction_list:
        sql = 'INSERT into enzyme_enzymes_dat (enzyme_dat_id,line, value) values(%s,"catalytic activity", "%s");' % (id, re.sub("<[a-z\/]+>", "", i))
        cursor.execute(sql,)

    # comments
    for i in cc:
        #i = re.sub("\r\r<b>Reaction:</b> ","",i)
        sql = 'INSERT into enzyme_enzymes_dat (enzyme_dat_id,line, value) values(%s, "comment", "%s");' % (id, re.sub("<[ =\"A-Za-z0-9\/\.]+>", "", i))
        cursor.execute(sql,)

    conn.commit()
    return id


def get_ko2ec(biodb):
    import urllib.request
    import re
    conn, cursor =connect_db(biodb)

    sql = 'CREATE TABLE IF NOT EXISTS enzyme_ko2ec (ko_id VARCHAR(20),' \
           ' enzyme_id INT,' \
           ' index enzyme_id (enzyme_id), ' \
           ' index ko_id (ko_id));'
    cursor.execute(sql)


    url = 'http://rest.kegg.jp/list/ko'
    data = [line for line in urllib.request.urlopen(url).read().decode('utf-8').split("\n")]

    total = len(data)
    #data = ["ko:K00512"]
    for n, ko_data in enumerate(data):
        if len(ko_data) == 0:
            continue
        print ("%s / %s" % (n, total))

        ko_id = ko_data.split('\t')[0][3:]
        print (ko_id)

        url_ko = 'http://rest.kegg.jp/link/ec/%s' % ko_id

        ko_data = [line.rstrip().split('\t') for line in urllib.request.urlopen(url_ko).read().decode('utf-8').split("\n")]

        for one_ec in ko_data:
            try:
                ec_name = one_ec[1][3:]
            except:
                print ('problem with:', ko_id, 'no ec?')
                continue

            try:
                sql_ec_id = 'select enzyme_id from enzyme.enzymes where ec="%s";' % ec_name
                cursor.execute(sql_ec_id,)
                ec_id = cursor.fetchall()[0][0]

                sql = 'INSERT INTO enzyme.ko2ec(ko_id, enzyme_id) VALUES ("%s",%s);' % (one_ec[0][3:], ec_id)
                #print sql
                cursor.execute(sql)

            except IndexError:
                print ('problem getting ec ID for:', ec_name, 'not in biosqldb?')
                print ('trying to add ec data from IUBMB')
                try:
                    ec_id = get_ec_data_from_IUBMB(ec_name)
                    print ("new_ec_id:", ec_id)
                    sql = 'INSERT INTO enzyme.ko2ec(ko_id, enzyme_id) VALUES ("%s",%s);' % (one_ec[0][3:], ec_id)
                    #print sql
                    cursor.execute(sql)
                except:
                    print ('FAIL')
        conn.commit()


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



def get_module_table_legacy(module2category, biodb):
    '''
    1. get all kegg pathways from API (http://rest.kegg.jp/) => create enzyme.kegg_module table
    2. get all KO associated for each module => create enzyme.module2ko table
    todo: remove existing tables for uptade if rerun

    :return: nothing
    '''


    import MySQLdb
    import urllib.request
    import re
    import sys
    import os

    sqlpsw = os.environ['SQLPSW']

    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db=biodb) # name of the data base
    cursor = conn.cursor()

    sql1 = 'CREATE TABLE IF NOT EXISTS enzyme_kegg_module_v1 (module_id INT AUTO_INCREMENT PRIMARY KEY,' \
           ' module_name VARCHAR(200),' \
           ' module_cat VARCHAR(200),' \
           ' module_sub_cat VARCHAR(200),' \
           ' module_sub_sub_cat VARCHAR(200),' \
           ' description LONG);'

    print (sql1)
    cursor.execute(sql1,)

    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme_module2ko_v1 (module_id INT,' \
           ' ko_id VARCHAR(200),' \
           ' ko_description TEXT);'


    sql3 = ' CONSTRAINT fk_pathway_id' \
           ' FOREIGN KEY(pathway_id) REFERENCES enzyme_kegg_pathway_v1(id)' \
           ' ON DELETE CASCADE,' \
           ' CONSTRAINT fk_ec_id' \
           ' FOREIGN KEY(ec_id) REFERENCES enzymes(id)' \
           ' ON DELETE CASCADE);'

    print (sql2)
    cursor.execute(sql2,)



    module_file_file = 'http://rest.kegg.jp/list/module'
    data = urllib.request.urlopen(module_file_file).read().decode('utf-8').split('\n')
    for line in data:
        if len(line) == 0:
            continue
        pathway = line.rstrip().split("\t")
        module = pathway[0][3:]
        description = pathway[1]
        try:
            cat = module2category[module][0]
            sub_cat = module2category[module][1]
            sub_sub_cat = module2category[module][2]
            description = module2category[module][3]
        except:
            cat = 'uncategorized'
            cat_short = 'uncategorized'
            print ('------------------------------------------------')
        sql = 'INSERT into enzyme_kegg_module_v1 (module_name, module_cat,module_sub_cat, module_sub_sub_cat, description) ' \
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
        ko_data = urllib.request.urlopen(ko_numbers_link).read().decode('utf-8').split('\n')

        for line in ko_data:
            if len(line) == 0:
                continue
            ko = line.rstrip().split("\t")[1][3:]
            try:
                ko_url = "http://rest.kegg.jp/get/%s" % ko
                ko_data = urllib.request.urlopen(ko_url).read().decode('utf-8').split('\n')
                definition = ko2definition(ko_data)
                sql = 'INSERT into enzyme_module2ko_v1 (module_id, ko_id, ko_description) values ("%s","%s", "%s");' % (module_id,
                                                                                                              ko,
                                                                                                              re.sub('"',
                                                                                                              '',
                                                                                                              definition))
                print (sql)
                cursor.execute(sql,)
            except:
                print (line)
                import sys
                sys.exit()
        conn.commit()


def get_pathway2ko_legacy(biodb):
    '''
    1. get all kegg pathways from API (http://rest.kegg.jp/) => create enzyme.kegg_module table
    2. get all KO associated for each pathway => create enzyme_pathway2ko table
    todo: remove existing tables for uptade if rerun

    :return: nothing
    '''


    import MySQLdb
    import urllib
    import re
    import sys

    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd="estrella3", # your password
                                db=biodb) # name of the data base
    cursor = conn.cursor()



    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme_pathway2ko_v1 (pathway_id INT,' \
           ' ko_id VARCHAR(200),' \
           ' index pathway_id(pathway_id),' \
           ' index ko_id(ko_id));'



    print (sql2)
    cursor.execute(sql2,)



    pathway_file = 'http://rest.kegg.jp/list/pathway'
    data = urllib.request.urlopen(pathway_file).read().decode('utf-8').split('\n')
    for line in data:
        if len(line) == 0:
            continue
        raw = line.rstrip().split("\t")
        pathway = raw[0][5:]
        description = raw[1]
        print ('pathway', pathway)
        sql = 'select pathway_id from enzyme_kegg_pathway_v1 where pathway_name="%s";' % pathway
        print (sql)

        cursor.execute(sql)

        try:
            pathway_id = cursor.fetchall()[0][0]
            print ('path id', pathway_id)
        except:
            print ('no pathway id for: %s, incomplete pathway table?' %  pathway)


        ko_numbers_link = "http://rest.kegg.jp/link/ko/%s" % pathway
        ko_data = urllib.request.urlopen(ko_numbers_link).read().decode('utf-8').split('\n')


        for line in ko_data:
            if len(line) == 0:
                continue
            try:
                ko = line.rstrip().split("\t")[1][3:]
            except IndexError:
                print ('No Ko for pathway %s?' % pathway)
                continue
            try:
                sql = 'INSERT into enzyme_pathway2ko_v1 (pathway_id, ko_id) values (%s,"%s");' % (pathway_id,
                                                                                       ko)
                cursor.execute(sql,)
            except:
                print (line)
                import sys
                sys.exit()
        conn.commit()


def get_pathay_table_legacy(map2category, biodb):
    import MySQLdb
    import urllib
    import sys

    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd="estrella3", # your password
                                db=biodb) # name of the data base
    cursor = conn.cursor()

    sql1 = 'CREATE TABLE IF NOT EXISTS enzyme_kegg_pathway_v1 (pathway_id INT AUTO_INCREMENT PRIMARY KEY,' \
           ' pathway_name VARCHAR(200),' \
           ' pathway_category_short VARCHAR(200),' \
           ' pathway_category VARCHAR(200),' \
           ' description LONG);'

    print (sql1)
    cursor.execute(sql1,)



    pathway_file_file = 'http://rest.kegg.jp/list/pathway'
    data = urllib.request.urlopen(pathway_file_file).read().decode('utf-8').split('\n')
    for line in data:
        if len(line) == 0:
            continue
        pathway = line.rstrip().split("\t")
        map = pathway[0][5:]
        description = pathway[1]

        try:
            cat = map2category[map][1]
            cat_short = map2category[map][0]
        except:
            cat = 'uncategorized'
            cat_short = 'uncategorized'
        sql = 'INSERT into enzyme_kegg_pathway_v1 (pathway_name, description,pathway_category_short, pathway_category) values ("%s", "%s", "%s", "%s");' % (map,
                                                                                                                                                  description,
                                                                                                                                                  cat_short,
                                                                                                                                                  cat)

        print (sql)
        cursor.execute(sql,)
        conn.commit()


def get_complete_ko_table_legacy(biodb):
    import urllib
    import re
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'CREATE TABLE IF NOT EXISTS enzyme_ko_annotation_v1 (ko_id VARCHAR(20),' \
           ' name varchar(200),' \
           ' definition TEXT,' \
           ' EC TEXT,' \
           ' pathways TEXT,' \
           ' modules TEXT, ' \
           ' dbxrefs TEXT, index ko_id (ko_id));'
    server.adaptor.execute_and_fetchall(sql)

    sql = 'select ko_id from enzyme_ko_annotation_v1'

    ko_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql)]
    print(ko_list)
    url = 'http://rest.kegg.jp/list/ko'
    data = urllib.request.urlopen(url).read().decode('utf-8').split('\n')

    total = len(data)
    for n, line in enumerate(data):
        if len(line) == 0:
            continue
        print ("%s / %s" % (n, total))
        ko = line.rstrip().split('\t')[0][3:]
        print (ko)
        if ko in ko_list:
            continue
        url_ko = 'http://rest.kegg.jp/get/%s' % ko

        ko_data = urllib.request.urlopen(url_ko).read().decode('utf-8').split('\n')

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

        sql = 'insert into enzyme_ko_annotation_v1 (ko_id, name, definition, EC, pathways, modules, dbxrefs)' \
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


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb

    parser = argparse.ArgumentParser()
    parser.add_argument("-u", '--updated', action="store_true", help="Setup updated tables")
    parser.add_argument("-l", '--legacy', action="store_true", help="Setup legacy tables")
    parser.add_argument("-d", '--db_name', help="Biodb name")

    args = parser.parse_args()

    conn, cursor = connect_db(args.db_name)

    print('getting map2category...')
    map2category = get_kegg_pathway_classification()

    if args.updated:
        print('getting complete_ko_table...')
        #get_complete_ko_table(args.db_name)
        print('load_enzyme_nomenclature_table map2category...')
        #load_enzyme_nomenclature_table(args.db_name)
        sql = 'select ko_accession, ko_id from enzyme_ko_annotation'
        cursor.execute(sql,)
        ko_accession2ko_id = manipulate_biosqldb.to_dict(cursor.fetchall())
        print('getting pathway table...')
        #get_pathay_table(map2category, args.db_name)
        print('getting ko2pathway...')
        #get_pathway2ko(ko_accession2ko_id, args.db_name)
        print('getting module2ko...')
        #module_hierarchy = get_kegg_module_hierarchy()
        #get_module_table(module_hierarchy, ko_accession2ko_id, args.db_name)
        #print('getting get_ec2get_pathway_table...')
        #get_ec2get_pathway_table(args.db_name)
        #print('getting get_microbial_metabolism_in_diverse_environments_kegg01120...')
        #get_microbial_metabolism_in_diverse_environments_kegg01120()
        print('getting get_ko2ec...')
        get_ko2ec(args.db_name)

        #get_ec_data_from_IUBMB("3.2.1.196")
        #get_ec_data_from_IUBMB("1.14.13.217")
        #get_ec_data_from_IUBMB("1.1.1.1")
    if args.legacy:
        get_complete_ko_table_legacy(args.db_name)
        get_pathay_table_legacy(map2category)
        get_pathway2ko_legacy()
        get_module_table_legacy(get_kegg_module_hierarchy())
