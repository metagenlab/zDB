#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# create html table from blast tabulated file
# headers: accession	size	gi	n proteins	n contigs	gc 	description
# add 4 columns with links of the form /assets/chlamdb/ffn/ for gbk/faa/ffm/fna
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2015
# ---------------------------------------------------------------------------

def get_interpro_entry_tables(interpro_release= '60.0'):
    '''

    :param interpro_release: par exemple: 60.0
    :return:
    '''
    import urllib2
    import MySQLdb
    import os
    sqlpsw = os.environ['SQLPSW']

    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="interpro") # name of the data base
    cursor = conn.cursor()



    sql = 'CREATE table entry (interpro_id INT AUTO_INCREMENT, name varchar(400), description TEXT, index interpro_id(interpro_id))'

    cursor.execute(sql,)
    conn.commit()
    link = 'ftp://ftp.ebi.ac.uk/pub/databases/interpro/%s/entry.list' % interpro_release
    print (link)
    req = urllib2.Request(link)
    entry_list = urllib2.urlopen(req)
    for line in entry_list:
        if not 'IPR' in line:
            continue
        accession = line.rstrip().split(' ')[0]
        description = ' '.join(line.rstrip().split(' ')[1:])
        sql = 'insert into entry (name, description) values ("%s", "%s")' % (accession, description)
        cursor.execute(sql,)
    conn.commit()

def get_interpro_parent2child_table(interpro_release):



    link = ''

def get_interpro2go_table():

    '''

    :param interpro_release: par exemple: 60.0
    :return:
    '''
    import urllib2
    import MySQLdb
    from chlamdb.biosqldb import manipulate_biosqldb
    import os
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="interpro") # name of the data base
    cursor = conn.cursor()



    sql = 'CREATE table interpro2gene_ontology (interpro_id INT, go_id INT, index interpro_id(interpro_id), index go_id(go_id))'
    #cursor.execute(sql,)
    #conn.commit()

    sql = 'select name,interpro_id from entry'
    cursor.execute(sql,)
    interpro_name2interpro_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    sql = 'select acc,id from gene_ontology.term where acc like "GO%"'
    cursor.execute(sql,)
    go_name2go_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    link = 'ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro2go'

    req = urllib2.Request(link)
    entry_list = urllib2.urlopen(req)
    echec = 0
    for i, line in enumerate(entry_list):
        if line[0] == '!':
            continue
        interpro_name = line.rstrip().split(':')[1].split(' ')[0]
        go_name = line.rstrip().split(' ; ')[1]
        try:
            sql = 'insert into interpro2gene_ontology values(%s, %s)' % (interpro_name2interpro_id[interpro_name], go_name2go_id[go_name])
            cursor.execute(sql,)
            conn.commit()
        except:
            echec+=1
    print (echec, i)


    # InterPro:IPR000003 Retinoid X receptor/HNF4 > GO:DNA binding ; GO:0003677


def interpro2biosqlV2(server,
                    interpro_accession2interpro_id,
                    seqfeature_id2locus_tag,
                    locus_tag2genome_taxon_id,
                    protein_id2genome_taxon_id,
                    locus_tag2seqfeature_id,
                    protein_id2seqfeature_id,
                    seqfeature_id2organism,
                    db_name,
                    *input_files):
    import re
    '''
    1. Protein Accession (e.g. P51587)
    2. Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
    3.. equence Length (e.g. 3418)
    4. Analysis (e.g. Pfam / PRINTS / Gene3D)
    5. Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
    6. Signature Description (e.g. BRCA2 repeat profile)
    7. Start location
    8. Stop location
    9. Score - is the e-value of the match reported by member database method (e.g. 3.1E-52)
    10. Status - is the status of the match (T: true)
    11. Date - is the date of the run
    12. (InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprscan option is switched on)
    13. (InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprscan option is switched on)
    14. (GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
    15. (Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)
    :param input_file:
    :return:
    '''

    sql = 'CREATE TABLE interpro_interpro (accession VARCHAR(100),' \
          ' seqfeature_id INT, ' \
          ' organism VARCHAR(200),  ' \
          ' taxon_id INT,' \
          ' sequence_length INT, ' \
          ' analysis VARCHAR(100) NOT NULL, ' \
          ' signature_accession VARCHAR(100), ' \
          ' signature_description VARCHAR(1000), ' \
          ' start INT, ' \
          ' stop INT, ' \
          ' score VARCHAR(20) NOT NULL, ' \
          ' interpro_id INT, ' \
          ' GO_terms varchar(10000),' \
          ' pathways varchar(10000),' \
          ' INDEX seqfeature_id (seqfeature_id),' \
          ' INDEX interpro_id (interpro_id),' \
          ' INDEX ia (interpro_accession))' % db_name
    try:
        server.adaptor.execute(sql)
    except:
        pass
    for one_interpro_file in input_files:
        print (one_interpro_file)
        from pandas import read_csv
        with open(one_interpro_file, 'r') as f:
            tsvin = read_csv(f, sep='\t', error_bad_lines=False, names=["accession",
                                                                        "MD5",
                                                                        "sequence_length",
                                                                        "analysis",
                                                                        "signature_accession",
                                                                        "signature_description",
                                                                        "start",
                                                                        "stop",
                                                                        "score",
                                                                        "status",
                                                                        "date",
                                                                        "interpro_accession",
                                                                        "interpro_description",
                                                                        "GO_terms",
                                                                        "pathways"])
            tsvin = tsvin.fillna(0)

            for i in range(len(tsvin['accession'])):

                data= list(tsvin.loc[i,:])
                for index, item in enumerate(data):
                    if type(item) == str:
                        data[index] = item.replace('\"','')

                accession = data[0]

                sequence_length = data[2]
                analysis = data[3]
                signature_accession = data[4]
                signature_description = data[5]
                start = data[6]
                stop = data[7]
                score = data[8]


                interpro_accession = data[11]
                interpro_description = data[12]
                GO_terms = data[13]
                pathways = data[14]

                try:
                    taxon_id = protein_id2genome_taxon_id[accession]
                    seqfeature_id = protein_id2seqfeature_id[accession]
                except KeyError:
                    taxon_id = locus_tag2genome_taxon_id[accession]
                    seqfeature_id = locus_tag2seqfeature_id[accession]
                organism = seqfeature_id2organism[str(seqfeature_id)]
                #print organism
                locus_tag = seqfeature_id2locus_tag[str(seqfeature_id)]

                sql = 'INSERT INTO interpro_interpro(accession, locus_tag, organism, taxon_id,' \
                      ' sequence_length, analysis, signature_accession, signature_description, start, ' \
                      ' stop, score, interpro_accession, interpro_description, GO_terms, pathways) ' \
                      ' values ("%s", "%s", "%s", %s, %s, "%s", "%s", "%s", %s, %s, "%s", "%s", "%s", "%s", "%s");' % (db_name,
                                                                                                     accession,
                                                                                                     seqfeature_id,
                                                                                                     organism,
                                                                                                     taxon_id,
                                                                                                       sequence_length,
                                                                                                       analysis,
                                                                                                       signature_accession,
                                                                                                       signature_description,
                                                                                                       int(start),
                                                                                                       int(stop),
                                                                                                       str(score),
                                                                                                       interpro_id,
                                                                                                       GO_terms,
                                                                                                       pathways)
                try:
                    server.adaptor.execute(sql)
                    server.adaptor.commit()
                except:
                    print (sql)
                    print (data)
                    import sys
                    sys.exit()


def update_analysis_dico(server):

    from chlamdb.biosqldb import manipulate_biosqldb

    sql = 'select analysis_name, analysis_id from interpro_analysis'

    analysis_nam2analysis_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    return analysis_nam2analysis_id

def interpro2biosql(server,
                    locus_tag2seqfeature_id,
                    db_name,
                    *input_files):
    import MySQLdb
    '''
    1. Protein Accession (e.g. P51587)
    2. Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
    3.. equence Length (e.g. 3418)
    4. Analysis (e.g. Pfam / PRINTS / Gene3D)
    5. Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
    6. Signature Description (e.g. BRCA2 repeat profile)
    7. Start location
    8. Stop location
    9. Score - is the e-value of the match reported by member database method (e.g. 3.1E-52)
    10. Status - is the status of the match (T: true)
    11. Date - is the date of the run
    12. (InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprscan option is switched on)
    13. (InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprscan option is switched on)
    14. (GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
    15. (Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)
    :param input_file:
    :return:
    '''

    sql = 'CREATE TABLE if not exists interpro_analysis(analysis_id INT AUTO_INCREMENT PRIMARY KEY, ' \
          ' analysis_name varchar(400),' \
          ' index analysis_name(analysis_name))'

    server.adaptor.execute(sql,)

    sql2 = 'CREATE TABLE if not exists interpro_signature (signature_id INT AUTO_INCREMENT PRIMARY KEY, ' \
           ' signature_accession varchar(400),' \
           ' signature_description TEXT,' \
           ' analysis_id INT,' \
           ' interpro_id INT, ' \
           ' GO_terms TEXT,' \
           ' pathways TEXT,' \
           ' INDEX analysis_id (analysis_id),' \
           ' index signature_accession(signature_accession),' \
           ' index interpro_id(interpro_id))' \

    server.adaptor.execute(sql2,)

    sql3 = 'CREATE TABLE if not exists interpro_interpro (seqfeature_id INT,' \
          ' sequence_length INT, ' \
          ' signature_id INT, ' \
          ' start INT, ' \
          ' stop INT, ' \
          ' score TEXT,' \
          ' INDEX signature_id (signature_id),' \
          ' INDEX seqfeature(seqfeature_id))' % db_name

    server.adaptor.execute(sql3)

    analysis2analysis_id = update_analysis_dico(server)
    sql = 'select signature_accession, signature_id from interpro_signature'
    signature2signature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    sql = 'select name, interpro_id from interpro_entry'
    interpro_entry2interpro_entry_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    for one_interpro_file in input_files:
        print (one_interpro_file)
        from pandas import read_csv
        with open(one_interpro_file, 'r') as f:
            tsvin = read_csv(f, sep='\t', error_bad_lines=False, names=["accession",
                                                                        "MD5",
                                                                        "sequence_length",
                                                                        "analysis",
                                                                        "signature_accession",
                                                                        "signature_description",
                                                                        "start",
                                                                        "stop",
                                                                        "score",
                                                                        "status",
                                                                        "date",
                                                                        "interpro_accession",
                                                                        "interpro_description",
                                                                        "GO_terms",
                                                                        "pathways"])
            tsvin = tsvin.fillna(0)

            for i in range(len(tsvin['accession'])):

                data = list(tsvin.loc[i,:])
                for index, item in enumerate(data):
                    if type(item) == str:
                        data[index] = item.replace('\"','')

                accession = data[0]

                sequence_length = data[2]
                analysis = data[3]
                signature_description = data[5]
                interpro_accession = data[11]
                interpro_description = data[12]
                GO_terms = data[13]
                pathways = data[14]

                try:
                    analysis_id = analysis2analysis_id[analysis]
                except KeyError:
                    print ('New analysis:', analysis)
                    sql = 'insert into interpro_analysis (analysis_name) values ("%s")' % analysis
                    server.adaptor.execute(sql)
                    server.adaptor.commit()
                    analysis2analysis_id = update_analysis_dico(server)
                    analysis_id = analysis2analysis_id[analysis]

                signature_accession = data[4]

                #sql = 'select signature_id from interpro_signature where signature_accession="%s"' % signature_accession

                try:
                    signature_id = signature2signature_id[signature_accession]#server.adaptor.execute_and_fetchall(sql,)[0][0]
                except KeyError:
                    print ('New signature', signature_accession, signature_description)
                    #sql1 = 'select interpro_id from interpro_entry where name="%s"' % (interpro_accession)

                    try:
                        interpro_id = interpro_entry2interpro_entry_id[interpro_accession]#server.adaptor.execute_and_fetchall(sql1,)[0][0]
                    except KeyError:
                        if interpro_accession == 0:
                            print ('No interpro-accession for ', signature_accession, signature_description)
                            interpro_id="NULL"
                        else:
                            print ('New Interpro entry', interpro_accession, interpro_description)

                            sql1b = 'insert into interpro_entry(name, description) values("%s","%s")' % (interpro_accession,
                                                                                                         interpro_description)
                            print (sql1b)
                            server.adaptor.execute(sql1b,)
                            server.adaptor.commit()
                            sql1 = 'select interpro_id from interpro_entry where name="%s"' % (interpro_accession)
                            interpro_id = server.adaptor.execute_and_fetchall(sql1,)[0][0]
                            # update dictionnray
                            interpro_entry2interpro_entry_id[interpro_accession] = interpro_id

                    sql2 = 'insert into interpro_signature (signature_accession, signature_description, ' \
                          ' analysis_id, interpro_id, GO_terms, pathways) values ("%s", "%s", %s, %s, "%s", "%s")' % (signature_accession,
                                                                                                     signature_description,
                                                                                                     analysis_id,
                                                                                                     interpro_id,
                                                                                                     GO_terms,
                                                                                                     pathways)

                    server.adaptor.execute(sql2,)
                    server.adaptor.commit()
                    sql = 'select signature_id from interpro_signature where signature_accession="%s"' % signature_accession
                    signature_id = server.adaptor.execute_and_fetchall(sql,)[0][0]
                    # update dictionnary
                    signature2signature_id[signature_accession] = signature_id

                start = data[6]
                stop = data[7]
                score = data[8]

                if analysis in ['Phobius', 'Coils', 'SignalP-TM', 'SignalP_EUK', 'ProSitePatterns', 'SignalP_GRAM_NEGATIVE', 'SignalP_GRAM_POSITIVE']:
                    score = "NULL"

                seqfeature_id = locus_tag2seqfeature_id[accession]

                sql = 'INSERT INTO interpro_interpro(seqfeature_id,' \
                      ' signature_id,' \
                      ' sequence_length, ' \
                      ' start, ' \
                      ' stop, ' \
                      ' score) ' \
                      ' values (%s, %s, %s, %s, %s, "%s");' % (db_name,
                                                             seqfeature_id,
                                                             signature_id,
                                                             sequence_length,
                                                             int(start),
                                                             int(stop),
                                                             score)
                try:
                    server.adaptor.execute(sql)

                except:
                    print (sql)
                    print (data)
                    import sys
                    sys.exit()
            server.adaptor.commit()


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_interpro', type=str, help="input interpro csv file", nargs='+')
    parser.add_argument("-d", '--database_name', type=str, help="database name")
    parser.add_argument("-v2", '--v2_table', action="store_true", help="create V2 table")

    args = parser.parse_args()

    biodb = args.database_name


    server, db = manipulate_biosqldb.load_db(biodb)
    #print "creating locus_tag2seqfeature_id"
    #locus_tag2seqfeature_id = manipulate_biosqldb.locus_tag2seqfeature_id_dict(server, biodb)



    #print "creating protein_id2seqfeature_id"
    #protein_id2seqfeature_id = manipulate_biosqldb.protein_id2seqfeature_id_dict(server, biodb)

    #print "getting seqfeature_id2organism"
    #seqfeature_id2organism = manipulate_biosqldb.seqfeature_id2organism_dico(server, biodb)

    #print "creating locus_tag2taxon_id dictionnary..."
    #locus_tag2genome_taxon_id = manipulate_biosqldb.locus_tag2genome_taxon_id(server, biodb)

    #print "creating protein_id2taxon_id dictionnary..."
    #protein_id2genome_taxon_id = manipulate_biosqldb.protein_id2genome_taxon_id(server, biodb)
    #print "getting seqfeature_id2locus_tag"
    #seqfeature_id2locus_tag = manipulate_biosqldb.seqfeature_id2locus_tag_dico(server, biodb)

    #print seqfeature_id2locus_tag.keys()[1:10]
    if not args.v2_table:

        sql = 'select locus_tag, seqfeature_id from annotation_seqfeature_id2locus' % biodb
        locus_tag2seqfeature_id =manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

        interpro2biosql(server,
                        locus_tag2seqfeature_id,
                        biodb, *args.input_interpro)

        #get_interpro2go_table()
    else:


        sql = 'select name,interpro_id from interpro_entry'

        interpro_accession2interpro_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

        interpro2biosqlV2(server,
                        interpro_accession2interpro_id,
                        seqfeature_id2locus_tag,
                        locus_tag2genome_taxon_id,
                        protein_id2genome_taxon_id,
                        locus_tag2seqfeature_id,
                        protein_id2seqfeature_id,
                        seqfeature_id2organism,
                        biodb, *args.input_interpro)