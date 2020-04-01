#!/usr/bin/env python


def connect_db(biodb):
    import os
    import re
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    try:
        conn.set_character_set('utf8')
        cursor.execute('SET NAMES utf8;')
        cursor.execute('SET CHARACTER SET utf8;')
        cursor.execute('SET character_set_connection=utf8;')
    except:
        cursor.execute('PRAGMA encoding="UTF-8";')


    return conn, cursor


def get_ko2modules():
    import urllib.request
    ko2modules = {}
    module2ko = {}
    url = 'http://rest.kegg.jp/link/module/ko'
    data = urllib.request.urlopen(url).read().decode('utf-8').split('\n') 

    for row in data:
        if len(row) == 0:
            continue
        # ko:K15376	md:M00880
        row_data = row.split("\t")
        ko = row_data[0].split(":")[1]
        module = row_data[1].split(":")[1]
        if ko not in ko2modules:
            ko2modules[ko] = [module]
        else:
            ko2modules[ko].append(module)

        if module not in module2ko:
            module2ko[module] = [ko]
        else:
            module2ko[module].append(ko)

    return ko2modules, module2ko


def get_ko2pathways():
    import urllib.request
    ko2pathways = {}
    pathway2ko = {}
    url = 'http://rest.kegg.jp/link/pathway/ko'
    data = urllib.request.urlopen(url).read().decode('utf-8').split('\n') 

    for row in data:
        if len(row) == 0:
            continue
        # ko:K15376	md:M00880
        row_data = row.split("\t")
        ko = row_data[0].split(":")[1]
        pathway = row_data[1].split(":")[1]
        if 'map' in pathway:
            # redudancy between path:map03018 and path:ko03018
            # only keep ko... syntax
            continue
        if ko not in ko2pathways:
            ko2pathways[ko] = [pathway]
        else:
            ko2pathways[ko].append(pathway)
        
        if pathway not in pathway2ko:
            pathway2ko[pathway] = [ko]
        else:
            pathway2ko[pathway].append(ko)
            
    return ko2pathways, pathway2ko


def get_pathway2ec():
    import urllib.request
    pathway2ec = {}
    url = 'http://rest.kegg.jp/link/pathway/ec'
    data = urllib.request.urlopen(url).read().decode('utf-8').split('\n') 

    for row in data:
        if len(row) == 0:
            continue
        # ko:K15376	md:M00880
        row_data = row.split("\t")
        ec = row_data[0].split(":")[1]
        pathway = row_data[1].split(":")[1]
        if not 'map' in pathway:
            # redudancy between path:map03018 and path:ec03018
            # only keep map... syntax
            continue
        
        if pathway not in pathway2ec:
            pathway2ec[pathway] = [ec]
        else:
            pathway2ec[pathway].append(ec)
            
    return pathway2ec


def get_ko2ec():
    import urllib.request
    ko2ec = {}
    url = 'http://rest.kegg.jp/link/ec/ko'
    data = urllib.request.urlopen(url).read().decode('utf-8').split('\n') 

    for row in data:
        if len(row) == 0:
            continue
        # ko:K15376	md:M00880
        row_data = row.split("\t")
        ko = row_data[0].split(":")[1]
        ec = row_data[1].split(":")[1]
        if ko not in ko2ec:
            ko2ec[ko] = [ec]
        else:
            ko2ec[ko].append(ec)
    return ko2ec


def get_ko2data():
    import urllib.request
    ko2name = {}
    ko2definition = {}

    url = 'http://rest.kegg.jp/list/ko'
    data = urllib.request.urlopen(url).read().decode('utf-8').split('\n')

    for row in data:
        if len(row) == 0:
            continue 
        row_data = row.split("\t")
        ko = row_data[0].split(":")[1]
        name = row_data[1].split(";")[0]
        definition = ';'.join(row_data[1].split(";")[1:])
        ko2name[ko] = name 
        ko2definition[ko] = definition

    return ko2name, ko2definition


def get_complete_ko_table(biodb,
                          ko2name,
                          ko2definition,
                          ko2ec,
                          ko2modules,
                          ko2pathways):
    import urllib.request
    import re
    
    conn, cursor = connect_db(biodb)

    sql = 'CREATE TABLE IF NOT EXISTS enzyme_ko_annotation (ko_id INTEGER PRIMARY KEY, ' \
          ' ko_accession VARCHAR(20),' \
          ' name varchar(60),' \
          ' definition TEXT,' \
          ' EC TEXT,' \
          ' pathways TEXT,' \
          ' modules TEXT, ' \
          ' dbxrefs TEXT);'
          
    cursor.execute(sql)

    sql = 'select ko_accession from enzyme_ko_annotation;'
    
    cursor.execute(sql,)

    ko_already_in_db = [i[0] for i in cursor.fetchall()]

    print('Already into DB:', ko_already_in_db[0:10])

    total = len(ko2name)
    
    for n, ko in enumerate(ko2name):
        if n % 1000 == 0:
            print("%s / %s" % (n, total))

        if ko in ko_already_in_db:
            print ('%s already in DB' % ko)
            continue

        definition = ko2definition[ko]
        name = ko2name[ko]
        
        try:
            ec = ','.join(ko2ec[ko])
        except KeyError:
            ex = '-'
        try:
            modules = ','.join(ko2modules[ko])
        except KeyError:
            modules = '-'

        try:
            pathways = ','.join(ko2pathways[ko])
        except KeyError:
            pathways = '-'

        sql = 'insert into enzyme_ko_annotation (ko_accession, name, definition, EC, pathways, modules, dbxrefs)' \
              ' values ("%s", "%s", "%s","%s", "%s", "%s", "%s")' % (ko,
                                                                     name,
                                                                     re.sub('"','',definition),
                                                                     ec,
                                                                     pathways,
                                                                     modules,
                                                                     '-')

        cursor.execute(sql,)
        conn.commit()

    # add indexes
    sql_index1 = 'create index ekaki on enzyme_ko_annotation(ko_id);'
    sql_index2 = 'create index ekaka on enzyme_ko_annotation(ko_accession);'
    
    cursor.execute(sql_index1)
    cursor.execute(sql_index2)
    conn.commit()
    

def load_enzyme_nomenclature_table(biodb):


    '''

    Download all SIB enzyme nomenclature from FTP (ftp://ftp.expasy.org/databases/enzyme/)
    create the enzyme_enzymes table with the list of all EC with associated description
    create the enzyme_enzymes_dat with detailed information about each EC
    todo: remove existing tables for uptade if rerun

    :return: nothing
    '''

    from Bio.ExPASy import Enzyme
    import urllib.request
    from io import StringIO
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor
    
    try:
        conn.set_character_set('utf8')
        cursor.execute('SET NAMES utf8;')
        cursor.execute('SET CHARACTER SET utf8;')
        cursor.execute('SET character_set_connection=utf8;')
    except:
        cursor.execute('PRAGMA encoding="UTF-8";')

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

    sql1 = 'CREATE TABLE IF NOT EXISTS enzyme_enzymes (enzyme_id INTEGER PRIMARY KEY,' \
          ' ec VARCHAR(200));'

    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme_enzymes_dat (enzyme_dat_id INT,' \
          ' line VARCHAR(20),' \
          ' value LONG,' \
           ' CONSTRAINT fk_enzyme_id' \
           ' FOREIGN KEY(enzyme_dat_id) REFERENCES enzyme_enzymes(enzyme_id)' \
           ' ON DELETE CASCADE);'

    print ('create enzyme table')

    cursor.execute(sql1,)
    print ('create dat table')

    cursor.execute(sql2)

    all_data = [i for i in Enzyme.parse(StringIO(data))]

    for n, data in enumerate(all_data):
        
        if n %  1000 == 0:
            print("%s / %s" % (n, len(all_data)))
        
        enzyme = data['ID']
        # insert enzyme id into primary TABLE
        sql = 'INSERT into enzyme_enzymes (ec) values ("%s");' % enzyme

        cursor.execute(sql,)
        conn.commit()

        id = cursor.lastrowid

        # description
        sql = 'INSERT into enzyme_enzymes_dat (enzyme_dat_id, line, value) values (%s, "description", "%s");' % (id, data['DE'])

        cursor.execute(sql,)
        # alternative names
        for i in data['AN']:
            sql = 'INSERT into enzyme_enzymes_dat (enzyme_dat_id,line, value) values(%s, "alternative name", "%s");' % (id, i)

            cursor.execute(sql,)

        # Catalytic activity
        sql = 'INSERT into enzyme_enzymes_dat (enzyme_dat_id,line, value) values(%s,"catalytic activity", "%s");' % (id, data['CA'])
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
    import urllib.request
    import sys
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor
    
    sql1 = 'CREATE TABLE IF NOT EXISTS enzyme_kegg_pathway (pathway_id INTEGER PRIMARY KEY,' \
           ' pathway_name VARCHAR(200),' \
           ' pathway_category_short VARCHAR(200),' \
           ' pathway_category VARCHAR(200),' \
           ' description LONG);'

    cursor.execute(sql1,)
    conn.commit()

    sql_index1 = 'create index ekpp on enzyme_kegg_pathway(pathway_name);'
    sql_index2 = 'create index ekppcs on enzyme_kegg_pathway(pathway_category_short);'
    sql_index3 = 'create index ekppc on enzyme_kegg_pathway(pathway_category);'
    
    server.adaptor.execute(sql_index1)
    server.adaptor.execute(sql_index2)
    server.adaptor.execute(sql_index3)
    conn.commit()

    pathway_file_file = 'http://rest.kegg.jp/list/pathway'
    data = urllib.request.urlopen(pathway_file_file).read().decode('utf-8').split("\n")

    print ('iter pathway list...')
    pathway_name2pathway_id = {}
    for line in data:
        # manage empty line(s)
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
        sql = 'INSERT into enzyme_kegg_pathway (pathway_name, description,' \
              ' pathway_category_short, pathway_category) values ("%s", "%s", "%s", "%s");' % (map,
                                                                                               description,
                                                                                               cat_short,
                                                                                               cat)

        cursor.execute(sql,)
        pathway_id = cursor.lastrowid
        pathway_name2pathway_id[map] = pathway_id
        
    conn.commit()
    


def get_pathway2ko(biodb,
                   ko_accession2ko_id, 
                   pathway2ko):
    '''
    1. get all kegg pathways from API (http://rest.kegg.jp/) => create enzyme_kegg_module table
    2. get all KO associated for each pathway => create enzyme_pathway2ko table
    todo: remove existing tables for uptade if rerun

    :return: nothing
    '''

    import urllib.request
    import re
    import sys
    from chlamdb.biosqldb import manipulate_biosqldb
    
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme_pathway2ko (pathway_id INT,' \
           ' ko_id INT);'

    cursor.execute(sql2,)
    conn.commit()

    sql = 'select pathway_name, pathway_id from enzyme_kegg_pathway'
    cursor.execute(sql,)
    pathway_name2pathway_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    pathway_file = 'http://rest.kegg.jp/list/pathway'
    data = urllib.request.urlopen(pathway_file).read().decode('utf-8').split("\n")
    
    for line in data:
        if len(line) == 0:
            continue
        raw = line.rstrip().split("\t")
        pathway = raw[0][5:]
        description = raw[1]

        cursor.execute(sql)

        pathway_id = pathway_name2pathway_id[pathway] 

        try:
            ko_list = pathway2ko[re.sub("map", "ko", pathway)]
        except KeyError:
            print("No KO data for pathway %s/%s, skipping..." % (pathway, description))
            continue

        for ko in ko_list:
            ko_id = ko_accession2ko_id[ko]
            sql = 'INSERT into enzyme_pathway2ko (pathway_id, ko_id) values (%s,%s);' % (pathway_id,
                                                                                         ko_id)
            cursor.execute(sql,)
        conn.commit()
        
    sql_index1 = 'create index ekpppid on enzyme_pathway2ko(pathway_id);'
    sql_index2 = 'create index ekpkid on enzyme_pathway2ko(ko_id);'
    
    server.adaptor.execute(sql_index1)
    server.adaptor.execute(sql_index2)
    conn.commit()


def get_module_table(biodb,
                     module2category,
                     ko_accession2ko_id):
    '''
    1. get all kegg pathways from API (http://rest.kegg.jp/) => create enzyme_kegg_module table
    2. get all KO associated for each module => create enzyme_module2ko table
    todo: remove existing tables for uptade if rerun

    :return: nothing
    '''

    import urllib.request
    import re
    import sys
    from chlamdb.biosqldb import manipulate_biosqldb
    
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    sql1 = 'CREATE TABLE IF NOT EXISTS enzyme_kegg_module (module_id INTEGER PRIMARY KEY,' \
           ' module_name VARCHAR(200),' \
           ' module_cat VARCHAR(200),' \
           ' module_sub_cat VARCHAR(200),' \
           ' module_sub_sub_cat VARCHAR(200),' \
           ' description LONG);'

    cursor.execute(sql1,)

    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme_module2ko (module_id INT,' \
           ' ko_id INT);'

    server.adaptor.execute(sql2)
   
    conn.commit()

    sqlm = 'select module_name from enzyme_kegg_module'
    cursor.execute(sqlm,)
    module_in_db = [i[0] for i in cursor.fetchall()]


    sql3 = ' CONSTRAINT fk_pathway_id' \
           ' FOREIGN KEY(pathway_id) REFERENCES enzyme_kegg_pathway(id)' \
           ' ON DELETE CASCADE,' \
           ' CONSTRAINT fk_ec_id' \
           ' FOREIGN KEY(ec_id) REFERENCES enzyme_enzymes(id)' \
           ' ON DELETE CASCADE);'

    cursor.execute(sql2,)

    module_file_file = 'http://rest.kegg.jp/list/module'
    data = urllib.request.urlopen(module_file_file).read().decode('utf-8').split("\n")
    
    
    module_name2module_id = {}
    
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

        cursor.execute(sql,)
        conn.commit()
        
        module_id = cursor.lastrowid
        
        module_name2module_id[module] = module_id
        

    # add indexes

    sql_index1 = 'create index ekmdn on enzyme_kegg_module(module_name);'
    sql_index2 = 'create index ekmdc on enzyme_kegg_module(module_cat);'
    sql_index3 = 'create index ekmdsc on enzyme_kegg_module(module_sub_cat);'
    sql_index4 = 'create index ekmdssc on enzyme_kegg_module(module_sub_sub_cat);'
    
    server.adaptor.execute(sql_index1)
    server.adaptor.execute(sql_index2)
    server.adaptor.execute(sql_index3)
    server.adaptor.execute(sql_index4)
    
    conn.commit()
    
    return module_name2module_id




def definition2module_list(definition_string):
    import re
    return re.findall('M[0-9]+', definition_string)


def ko2definition(ko_record):
    import re

    m = re.compile("^DEFINITION")

    for line in ko_record:
        if re.match(m, line):
            definition=line.rstrip().split('  ')[1]
            return definition


def get_multimodule(module):
    import urllib.request
    
    module_link = "http://rest.kegg.jp/get/%s" % module
    module_data = urllib.request.urlopen(module_link).read().decode('utf-8').split("\n")
    definition = ko2definition(module_data)
    module_list = definition2module_list(definition)
    ko_list = []
    for n, m in enumerate(module_list):
        ko_list += module2ko[m]
    return ko_list



def get_module2ko(biodb,
                  module2ko,
                  module_name2module_id):


    for module in module_name2module_id:        

        try:
            ko_list = module2ko[module]
        except KeyError:
            print ('MODULE MADE OF SUBMODULES----------')
            ko_list = get_multimodule(module)
        
        module_id = module_name2module_id[module]

        for ko in ko_list:

            ko_id = ko_accession2ko_id[ko]          
            
            sql = 'INSERT into enzyme_module2ko (module_id, ko_id) values (%s,%s);' % (module_id,
                                                                                       ko_id)

            cursor.execute(sql,)

        conn.commit()
    
    # add indexes
    sql_index1 = 'create index emkkid on enzyme_module2ko(ko_id);'
    sql_index2 = 'create index emkmid on enzyme_module2ko(module_id);'

    cursor.execute(sql_index1)
    cursor.execute(sql_index2)
    conn.commit()


def get_module2ko_legacy(biodb,
                         module2ko,
                         module_name2module_id,
                         ko2definition):


    for module in module_name2module_id:        

        try:
            ko_list = module2ko[module]
        except KeyError:
            print ('MODULE MADE OF SUBMODULES----------')
            ko_list = get_multimodule(module)

        for ko in ko_list:       
            
            definition = ko2definition[ko]
            
            module_id = module_name2module_id[module]
            
            sql = 'INSERT into enzyme_module2ko_v1 (module_id, ko_id, ko_description) values (%s,"%s", "%s");' % (module_id,
                                                                                                                  ko,
                                                                                                                  definition)

            cursor.execute(sql,)

        conn.commit()
    
    # add indexes
    sql_index1 = 'create index emkkidv1 on enzyme_module2ko_v1(ko_id);'
    sql_index2 = 'create index emkmidv1 on enzyme_module2ko_v1(module_id);'

    cursor.execute(sql_index1)
    cursor.execute(sql_index2)
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


def setup_kegg_pathway2ec(biodb,
                          pathway2ec):
    '''
    2. get all ec associated for each pathway => create enzyme_kegg2ec table
    todo: remove existing tables for uptade if rerun

    :return: nothing
    '''


    import urllib.request
    import sys
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

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

    cursor.execute(sql2,)
    conn.commit()


    sql = 'select pathway_name,pathway_id from enzyme_kegg_pathway;'
    cursor.execute(sql,)
    
    pathay2pathway_id = manipulate_biosqldb.to_dict(cursor.fetchall())
    
    sql = 'select ec,enzyme_id from enzyme_enzymes;'
    cursor.execute(sql,)
    
    ec2ec_id = manipulate_biosqldb.to_dict(cursor.fetchall())    
    
    for map in pathay2pathway_id:
        id = pathay2pathway_id[map]

        try:
            ec_list = pathway2ec[map]
        except KeyError:
            print("No EC for pathway %s, skipping" % map)
            continue
        
        for ec in ec_list:
            try:
                ec_id = ec2ec_id[ec]
            except KeyError:
                # commit to avoid sqlite lock
                conn.commit()
                try:
                    ec_id = get_ec_data_from_IUBMB(ec, biodb)
                    # update dictionnary
                    ec2ec_id[ec_name] = ec_id
                except:
                    print("FAIL ==> might be a wrong EC mumber, slipping")
                    continue
            sql = 'INSERT into enzyme_kegg2ec (pathway_id, ec_id) values (%s,"%s");' % (id, ec_id)
            cursor.execute(sql,)
            
    conn.commit()



def get_ec_data_from_IUBMB(ec, 
                           biodb):
    import urllib.request
    import re
    from bs4 import BeautifulSoup
    from chlamdb.biosqldb import manipulate_biosqldb
    
    server, db = manipulate_biosqldb.load_db(biodb)
    
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    name_m = re.compile(u".*Accepted name.*")
    alname = re.compile(u".*Other name.*")
    reaction = re.compile(u".*Reaction:\<\/b\>.*")
    comments = re.compile(u".*Comments.*")

    try:
        conn.set_character_set('utf8')
        cursor.execute('SET NAMES utf8;')
        cursor.execute('SET CHARACTER SET utf8;')
        cursor.execute('SET character_set_connection=utf8;')
    except:
        cursor.execute('PRAGMA encoding="UTF-8";')

    ec_sep = ec.split('.')

    adress = "https://www.qmul.ac.uk/sbcs/iubmb/enzyme/EC%s/%s/%s/%s.html" % (ec_sep[0], ec_sep[1], ec_sep[2], ec_sep[3])
    print(adress)

    html = urllib.request.urlopen(adress).read().decode('utf-8')

    #html = urllib2.urlopen(adress).read()
    html = re.sub("\&\#", "-", html)
    soup = BeautifulSoup(html, "lxml")
    html = str(soup.encode('utf-8'))#.encode('latin-1') #encode('utf-8') # prettify()

    #all_data = soup.find_all("p")#[i.get_text() for i in soup.find_all("p")]
    # empty list in case no match was found
    altname = []
    reaction_list = []
    cc = []
    
    for i, data in enumerate(list(html.split('<p>'))):
        if re.match(name_m, data):
            name = data.split("</b>")[-1]
            #name = re.sub("&alpha;","", name)
        elif re.match(alname, data):
            print("altname")
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

    sql_new = 'INSERT into enzyme_enzymes (ec) values ("%s");' % ec
    cursor.execute(sql_new, )
    conn.commit()

    id = cursor.lastrowid

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


def get_ko2ec_table(biodb, 
                    ko2ec):
    
    import urllib.request
    import re
    conn, cursor = connect_db(biodb)

    sql = 'CREATE TABLE IF NOT EXISTS enzyme_ko2ec (ko_id VARCHAR(20),' \
           ' enzyme_id INT);'
    cursor.execute(sql)
    
    sql = 'select ec, enzyme_id from enzyme_enzymes;'
    cursor.execute(sql,)
    
    ec2ec_id = manipulate_biosqldb.to_dict(cursor.fetchall())    
       
    for n, ko_id in enumerate(ko2ec):

        ec_list = ko2ec[ko_id]

        for ec_name in ec_list:

            try:
                ec_id = ec2ec_id[ec_name]
            except KeyError:
                print ('problem getting ec ID for:', ec_name, 'not in biosqldb?')
                print ('Trying to add ec data from IUBMB')
                try:
                    ec_id = get_ec_data_from_IUBMB(ec_name, biodb)
                    # update dictionnary
                    ec2ec_id[ec_name] = ec_id
                except:
                    print ('FAIL ==> might be a wrong EC number, skipping')
                    
            sql = 'INSERT INTO enzyme_ko2ec(ko_id, enzyme_id) VALUES ("%s", %s);' % (ko_id, ec_id)

            cursor.execute(sql)
            
        conn.commit()

    sql_index1 = 'create index ekeceid on enzyme_ko2ec(enzyme_id);'
    sql_index2 = 'create index ekeckid on enzyme_ko2ec(ko_id);'

    cursor.execute(sql_index1)
    cursor.execute(sql_index2)
    
    conn.commit()

def get_module_table_legacy(module2category, 
                            biodb):
    '''
    1. get all kegg pathways from API (http://rest.kegg.jp/) => create enzyme_kegg_module table
    2. get all KO associated for each module => create enzyme_module2ko table
    todo: remove existing tables for uptade if rerun

    :return: nothing
    '''

    import urllib.request
    import re
    import sys
    from chlamdb.biosqldb import manipulate_biosqldb
    
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    sql1 = 'CREATE TABLE IF NOT EXISTS enzyme_kegg_module_v1 (module_id INTEGER PRIMARY KEY,' \
           ' module_name VARCHAR(200),' \
           ' module_cat VARCHAR(200),' \
           ' module_sub_cat VARCHAR(200),' \
           ' module_sub_sub_cat VARCHAR(200),' \
           ' description LONG);'

    cursor.execute(sql1,)

    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme_module2ko_v1 (module_id INT,' \
           ' ko_id VARCHAR(200),' \
           ' ko_description TEXT);'


    sql3 = ' CONSTRAINT fk_pathway_id' \
           ' FOREIGN KEY(pathway_id) REFERENCES enzyme_kegg_pathway_v1(id)' \
           ' ON DELETE CASCADE,' \
           ' CONSTRAINT fk_ec_id' \
           ' FOREIGN KEY(ec_id) REFERENCES enzyme_enzymes(id)' \
           ' ON DELETE CASCADE);'

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


        cursor.execute(sql,)
        conn.commit()
        
    conn.commit()


def get_pathway2ko_legacy(biodb,
                          pathway2ko):
    '''
    1. get all kegg pathways from API (http://rest.kegg.jp/) => create enzyme_kegg_module table
    2. get all KO associated for each pathway => create enzyme_pathway2ko table
    todo: remove existing tables for uptade if rerun

    :return: nothing
    '''

    import urllib
    import re
    import sys
    from chlamdb.biosqldb import manipulate_biosqldb
    
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor


    sql2 = 'CREATE TABLE IF NOT EXISTS enzyme_pathway2ko_v1 (pathway_id INT,' \
           ' ko_id VARCHAR(200));'
    cursor.execute(sql2,)
    
    sql = 'select pathway_name, pathway_id from enzyme_kegg_pathway_v1'
    cursor.execute(sql,)
    pathway_name2pathway_id = manipulate_biosqldb.to_dict(cursor.fetchall())
       
    pathway_file = 'http://rest.kegg.jp/list/pathway'
    data = urllib.request.urlopen(pathway_file).read().decode('utf-8').split('\n')
    for line in data:
        if len(line) == 0:
            continue
        raw = line.rstrip().split("\t")
        pathway = raw[0][5:]
        description = raw[1]
        
        pathway_id = pathway_name2pathway_id[pathway]

        try:
            ko_list = pathway2ko[re.sub("map", "ko", pathway)]
        except KeyError:
            print("No KO data for pathway %s/%s, skipping..." % (pathway, 
                                                                 description))
            continue

        for ko in ko_list:
            sql = 'INSERT into enzyme_pathway2ko_v1 (pathway_id, ko_id) values (%s, "%s");' % (pathway_id,
                                                                                                ko)
            cursor.execute(sql,)

        conn.commit()

    sql_index1 = 'create index epakov1pid on enzyme_pathway2ko_v1(pathway_id);'
    sql_index2 = 'create index epakov1kid on enzyme_pathway2ko_v1(ko_id);'

    server.adaptor.execute(sql_index1)
    server.adaptor.execute(sql_index2)
    
    conn.commit()


def get_pathay_table_legacy(map2category, biodb):
    import urllib
    import sys
    from chlamdb.biosqldb import manipulate_biosqldb
    
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    sql1 = 'CREATE TABLE IF NOT EXISTS enzyme_kegg_pathway_v1 (pathway_id INTEGER PRIMARY KEY,' \
           ' pathway_name VARCHAR(200),' \
           ' pathway_category_short VARCHAR(200),' \
           ' pathway_category VARCHAR(200),' \
           ' description LONG);'

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

        cursor.execute(sql,)
    conn.commit()


def get_complete_ko_table_legacy(biodb,
                                 ko2name,
                                 ko2definition,
                                 ko2ec,
                                 ko2modules,
                                 ko2pathways):
    import urllib
    import re
    
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'CREATE TABLE IF NOT EXISTS enzyme_ko_annotation_v1 (ko_id VARCHAR(20),' \
           ' name varchar(200),' \
           ' definition TEXT,' \
           ' EC TEXT,' \
           ' pathways TEXT,' \
           ' modules TEXT, ' \
           ' dbxrefs TEXT);'
           
    server.adaptor.execute(sql)   

    sql = 'select ko_id from enzyme_ko_annotation_v1'

    ko_already_in_db = [i[0] for i in server.adaptor.execute_and_fetchall(sql)]

    total = len(ko2name)
    
    for n, ko in enumerate(ko2name):
        if n % 1000 == 0:
            print("%s / %s" % (n, total))

        if ko in ko_already_in_db:
            print ('%s already in DB' % ko)
            continue

        definition = ko2definition[ko]
        name = ko2name[ko]
        
        try:
            ec = ','.join(ko2ec[ko])
        except KeyError:
            ex = '-'
        try:
            modules = ','.join(ko2modules[ko])
        except KeyError:
            modules = '-'

        try:
            pathways = ','.join(ko2pathways[ko])
        except KeyError:
            pathways = '-'

        sql = 'insert into enzyme_ko_annotation_v1 (ko_id, name, definition, EC, pathways, modules, dbxrefs)' \
              ' values ("%s", "%s", "%s","%s", "%s", "%s", "%s")' % (ko,
                                                                     name,
                                                                     re.sub('"','', definition),
                                                                     ec,
                                                                     pathways,
                                                                     modules,
                                                                     "-")

        server.adaptor.execute(sql,)
        
    server.commit()

    sql_index = 'create index ekav1k on enzyme_ko_annotation_v1(ko_id);'
    
    server.adaptor.execute(sql_index)
    conn.commit()
    

if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb

    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--db_name', help="Biodb name")

    args = parser.parse_args()

    conn, cursor = connect_db(args.db_name)
    
    print('Retrieve ko2name, ko2definition...')
    # complete ko list
    ko2name, ko2description = get_ko2data()

    print('Retrieve ko2ec...')
    # ko2EC
    ko2ec = get_ko2ec() 
    
    print('Retrieve ko2modules...')
    # ko2modules
    ko2modules, module2ko = get_ko2modules()
    
    print('Retrieve ko2pathways...')
    # ko2pathways
    ko2pathways, pathway2ko = get_ko2pathways()
    

    print('Setup ko table...')
    get_complete_ko_table(args.db_name,
                            ko2name,
                            ko2description,
                            ko2ec,
                            ko2modules,
                            ko2pathways)
    
    
    
    print('Setup enzyme nomenclature table...')
    load_enzyme_nomenclature_table(args.db_name)

    sql = 'select ko_accession, ko_id from enzyme_ko_annotation'
    cursor.execute(sql,)
    
    ko_accession2ko_id = manipulate_biosqldb.to_dict(cursor.fetchall())
    
    map2category = get_kegg_pathway_classification()
    
    print('Setup pathway table...')
    get_pathay_table(map2category, 
                        args.db_name)
    
    print('Setup ko2pathway...')
    get_pathway2ko(args.db_name,
                    ko_accession2ko_id, 
                    pathway2ko)       
    
    print('Setup module table...')
    module_hierarchy = get_kegg_module_hierarchy()
    
    module_name2module_id = get_module_table(args.db_name,
                                                module_hierarchy, 
                                                ko_accession2ko_id)
    
    
    print('Setup module2ko...')       
    get_module2ko(args.db_name,
                    module2ko,
                    module_name2module_id)   
    
    print('Setup pathway2ec table...')
    pathway2ec = get_pathway2ec()

    setup_kegg_pathway2ec(args.db_name, 
                            pathway2ec)

    #print('getting get_microbial_metabolism_in_diverse_environments_kegg01120...')
    #get_microbial_metabolism_in_diverse_environments_kegg01120()
    
    #get_ec_data_from_IUBMB("4.6.1.24", args.db_name)
    
    print('Setup get_ko2ec...')
    get_ko2ec_table(args.db_name,
                    ko2ec)


    # TODO: stop using legacy enzyme tables and remove them 
    print('Setup ko legacy table...')
    get_complete_ko_table_legacy(args.db_name,
                                    ko2name,
                                    ko2description,
                                    ko2ec,
                                    ko2modules,
                                    ko2pathways)
    
    
    print('Setup pathway legacy table...')
    map2category = get_kegg_pathway_classification()
            
    map2category = get_kegg_pathway_classification()
    get_pathay_table_legacy(map2category, 
                            args.db_name)
            
    print('Setup pathway2ko legacy table...')
    get_pathway2ko_legacy(args.db_name,
                            pathway2ko)       
    
    print('Setup module legacy table...')
    module_hierarchy = get_kegg_module_hierarchy()
    
    get_module_table_legacy(module_hierarchy,
                            args.db_name)
    
    print('Setup module2ko legacy table...')
    sql = 'select module_name, module_id from enzyme_kegg_module_v1'
    cursor.execute(sql,)
    
    module_name2module_id = manipulate_biosqldb.to_dict(cursor.fetchall())
    
    get_module2ko_legacy(args.db_name,
                            module2ko,
                            module_name2module_id,
                            ko2description)
    

