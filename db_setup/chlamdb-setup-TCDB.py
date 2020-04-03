#!/usr/bin/env python


from chlamdb.eutils import tcdb_utils
import MySQLdb
import re
import os
from chlamdb.eutils import accession2taxon_id
from chlamdb.eutils import sequence_id2scientific_classification
from Bio import SeqIO
from bs4 import BeautifulSoup
from io import StringIO
import urllib.request
from chlamdb.biosqldb import blastswiss2sqltable


class TCDB():
    def __init__(self, biodb):
        from chlamdb.biosqldb import manipulate_biosqldb
        
        server, db = manipulate_biosqldb.load_db(biodb)
        self.conn = server.adaptor.conn
        self.cursor = server.adaptor.cursor

        try:
            self.conn.set_character_set('utf8mb4')
            self.cursor.execute('SET NAMES utf8mb4;')
            self.cursor.execute('SET CHARACTER SET utf8mb4;')
            self.cursor.execute('SET character_set_connection=utf8mb4;')
            self.cursor.execute('SET character_set_database=utf8mb4;')
            self.cursor.execute('SET character_set_server=utf8mb4;')
        except:
            self.cursor.execute('PRAGMA encoding="utf8mb4";')
        self.conn.commit()

        # for transporter entries: eg 2.A.88.1.2 (full path)
        sql = 'CREATE table if not exists transporters_transporters (transporter_id INTEGER primary key, tc_name varchar(200), superfamily INT,' \
              ' family INT, subfamily INT, substrate TEXT, description TEXT, domain varchar(200), superkingdom varchar(200))'

        self.cursor.execute(sql,)
        
        sql_index1 = 'create index ttttc1 on transporters_transporters(tc_name)'
        sql_index2 = 'create index tttspf on transporters_transporters(superfamily)'
        sql_index3 = 'create index tttf on transporters_transporters(family)'
        sql_index4 = 'create index tttsuf on transporters_transporters(subfamily)'
        self.cursor.execute(sql_index1,)
        self.cursor.execute(sql_index2,)
        self.cursor.execute(sql_index3,)
        self.cursor.execute(sql_index4,)

        self.conn.commit()

        # only for categories (class, families, subfamilies )
        # 1.B.1.1
        sql2 = 'CREATE table if not EXISTS transporters_classification (tc_id INTEGER primary key, tc_name varchar(400), description TEXT)'
        self.cursor.execute(sql2,)
        self.conn.commit()

        sql = 'CREATE TABLE if not exists transporters_protein_entry (protein_id INTEGER primary key, protein_accession varchar(400), organism TEXT)'
              
        sql2 = 'CREATE TABLE if not exists transporters_protein_entry2transporter (protein_id INTEGER , transporter_id INTEGER)'
                       
        sql_index1 = 'create index tpetpid on transporters_protein_entry2transporter(protein_id)'
        sql_index2 = 'create index tpettid on transporters_protein_entry2transporter(transporter_id)'
        self.cursor.execute(sql,)
        self.cursor.execute(sql2,)
        self.cursor.execute(sql_index1,)
        self.cursor.execute(sql_index2,)
        self.conn.commit()


    def insert_classification(self,
                              tcid2definition):
        
        for tcid in tcid2definition:
            sql = 'insert into transporters_classification (tc_name, description) values ("%s", "%s")' % (tcid, tcid2definition[tcid])
            self.cursor.execute(sql,)
            
        self.conn.commit()
        
    
    def insert_missing_tc_id(self, tc_id):
        
        sql = 'insert into transporters_classification (tc_name, description) values ("%s", "%s")' % (tc_id, "n/a")
        self.cursor.execute(sql,)
        self.conn.commit()
        id = self.cursor.lastrowid
        return id
        

    def insert_transporters(self, 
                            acc2tcid_list,
                            tcid2data,
                            tcid2substrate):


        from chlamdb.biosqldb import manipulate_biosqldb
        
        sql = 'select tc_name, tc_id from transporters_classification'
        
        tc_classif2classification_id = manipulate_biosqldb.to_dict(self.cursor.execute(sql,).fetchall())
        
        tcid_lists = list(acc2tcid_list.values())
        
        nr_tcid_list = list(set([tcid for tcid_list in tcid_lists for tcid in tcid_list]))

        print("Loading %s transporters into DB..." % (len(nr_tcid_list)))
        for tc_name in nr_tcid_list:

            separated_ids = tc_name.split('.')
            
            if len(separated_ids) != 5:
                print ('WARNING invalid id: %s, skipping...' % tc_name)
                continue
            else:
                id1 = separated_ids[0]
                tc_id_superfamily = '.'.join(separated_ids[0:2])
                tc_id_family = '.'.join(separated_ids[0:3])
                tc_id_subfamily = '.'.join(separated_ids[0:4])

                # deal with eventual missing entries
                try:
                    id_superfamily = tc_classif2classification_id[tc_id_superfamily]
                except KeyError:
                    id_superfamily = self.insert_missing_tc_id(tc_id_superfamily)
                    tc_classif2classification_id[tc_id_superfamily] = id_superfamily
                
                try:
                    id_family = tc_classif2classification_id[tc_id_family]
                except KeyError:
                    id_family = self.insert_missing_tc_id(tc_id_family)
                    tc_classif2classification_id[tc_id_family] = id_family
                
                try:
                    id_subfamily = tc_classif2classification_id[tc_id_subfamily]
                except KeyError:
                    id_subfamily = self.insert_missing_tc_id(tc_id_subfamily)
                    tc_classif2classification_id[tc_id_subfamily] = id_subfamily
                
            try:    
                description = tcid2data[tc_name]["name"]
            except:
                description = 'n/a'
            try:
                domain = tcid2data[tc_name]["domain"]
            except:
                domain = 'n/a'
            try:
                superkingdom = tcid2data[tc_name]["superkingdom"]
            except:
                superkingdom = 'n/a'
            try:
                substrate = tcid2substrate[tc_name]
            except:
                substrate = "n/a"


            sql = 'insert into transporters_transporters (tc_name, superfamily, family, subfamily, substrate, description, domain, superkingdom) values ("%s",' \
                    ' %s,%s,%s,"%s","%s","%s","%s")' % (tc_name, 
                                                        id_superfamily, 
                                                        id_family, 
                                                        id_subfamily, 
                                                        substrate, 
                                                        re.sub('"', "",description), 
                                                        domain, 
                                                        superkingdom)
            #print(sql)
            self.cursor.execute(sql,)
        self.conn.commit()


    def import_annot(self,
                     acc2tcid_list,
                     acc2species):
        from chlamdb.biosqldb import manipulate_biosqldb
        sql = 'select tc_name, transporter_id from transporters_transporters'

        transporter2transporter_id = manipulate_biosqldb.to_dict(self.cursor.execute(sql,).fetchall())
        
        for n, accession in enumerate(acc2tcid_list):
            tcid_list = acc2tcid_list[accession]
            try:
                species = acc2species[accession]
            except:
                species = 'n/a'

            # insert transporter                
            sql = 'insert into transporters_protein_entry (protein_accession, organism) values ("%s", "%s")' % (accession, species)

            self.cursor.execute(sql,)
            
            protein_id = self.cursor.lastrowid
                
            # the same protein can be associated with muitiple tcid 
            for tcid in tcid_list:
                try:
                    transporter_id = transporter2transporter_id[tcid]
                except:
                    print("ERROR: missing id for TCID %s, skipping" % tcid)
                    continue
                    
                sql = 'insert into transporters_protein_entry2transporter (protein_id, transporter_id) values (%s, %s)' % (protein_id, transporter_id)
                #print(sql)
                self.cursor.execute(sql,)
        self.conn.commit()

import re

def cleanhtml(raw_html):
  cleanr = re.compile('<.*?>')
  cleantext = re.sub(cleanr, '', raw_html)
  return cleantext


def get_family2description():
    import urllib.request
    
    families = 'http://www.tcdb.org/cgi-bin/projectv/public/families.py'
    
    data = urllib.request.urlopen(families).read().decode('utf-8').split('\n')
    family2description = {}
    for row in data:
        if len(row) == 0:
            continue
        try:
            family, description = row.split("\t")
        except:
            row_data = row.split("\t")
            family = row_data[0]
            description = ' '.join(row_data[1:])
        family2description[family] = cleanhtml(description)
    
    return family2description


def parse_tcdb_fasta():
    import re
    print("Download fasta db")
    handle = urllib.request.urlopen("http://www.tcdb.org/public/tcdb")
    fasta_st = StringIO(handle.read().decode("UTF-8"))
    records = SeqIO.parse(fasta_st, 'fasta')
    
    acc2species = {}
    missing = []
    for n, record in enumerate(records):
        # gnl|TC-DB|1001796365|4.F.1.1.5
        header_data = record.name.split("|")
        hit_uniprot_accession = header_data[2]
        description = record.description
        try:
            hit_tcid = re.findall("[0-9]+\.[A-Z]\.[0-9]+\.[0-9]+\.[0-9]+", description)[0]
        except:
            print("WARNING malformed entry:", description)
            hit_tcid = description.split("|")[3].split(" ")[0]
            print("==> considering %s as TCID" % hit_tcid)
        # extract species from fasta header
        # 3 formats
        # 1. >gnl|TC-DB|1001796365|4.F.1.1.5 CDP-alcohol phosphatidyltransferase [Marinobacter excellens]
        # 2. Amurin-1 protein OS=Rana amurensis GN=amurin-1 PE=2 SV=1
        # 3. >gnl|TC-DB|A2ARV4|9.B.87.1.1 Low-density lipoprotein receptor-related protein 2 - Mus musculus (Mouse).
        species = re.findall("\[(.*)\]", description)
        if len(species) > 0:
            acc2species[hit_uniprot_accession] = species[0]
            continue 
        species = re.findall("OS=([^\s]+\s+[^\s]+)", description)
        if len(species) > 0:
            acc2species[hit_uniprot_accession] = species[0]
            continue        
        species = re.findall("- (.*$)", description)    
        if len(species) > 0:
            acc2species[hit_uniprot_accession] = species[0]
            continue 
        else:
            print("WARNING: species could not be found:", description)
            acc = description.split("|")[2]
            missing.append(acc)
            acc2species[hit_uniprot_accession] = 'Unknown'

    return acc2species, missing


def get_superfamilies(acc='2.A'):
    
    id2description = {}
    
    import urllib.request
    
    url= 'http://www.tcdb.org/search/result.php?tc=%s' % acc
    
    page = urllib.request.urlopen(url)
    
    soup = BeautifulSoup(page, 'html.parser')
    
    box = soup.find("div", {"class": "next-level-result"})
    
    title = soup.find("div", {"class": "title"})
    
    # extract tcid and description
    tcid = cleanhtml(title.text).split(":")[0]
    description = ':'.join(cleanhtml(title.text).split(":")[1:]).strip()
    id2description[tcid] = description
    
    superfam_list = box.findAll("div", {"class": "next-level-links"})
    
    for superfam in superfam_list:
        # <div class="next-level-links">2.A.132 <a href="/search/result.php?tc=2.A.132">The Ferrous Iron Transporter (IroT/MavN) Family </a></div>
        tcid = re.findall("(^[0-9]+\.[A-Z]) ", superfam.text.strip())[0]
        description = re.findall("%s (.*)" % tcid, superfam.text.strip())[0]
        id2description[tcid] = description
    return id2description

def get_acc2tcid_list():
    import urllib.request
    
    acc2tcid = 'http://www.tcdb.org/cgi-bin/projectv/public/acc2tcid.py'

    data = urllib.request.urlopen(acc2tcid).read().decode('utf-8').split('\n')
    
    acc2tcid = {}
    
    for row in data:
        if len(row) == 0:
            continue
        row_data = [i.strip() for i in row.split("\t") if len(i) != 0]
        if len(row_data) != 2:
            print("WARNING incomplete row: %s" % row)
        else:
            acc, tcid = row_data
            if acc not in acc2tcid:
                acc2tcid[acc] = [tcid]
            else:
                acc2tcid[acc].append(tcid)
    
    return acc2tcid


def get_tcid2substrate():
    import urllib.request
    
    url = 'http://tcdb.org/cgi-bin/projectv/getSubstrates.py'

    data = urllib.request.urlopen(url).read().decode('utf-8').split('\n')
    
    tcid2substrate = {}
    
    for row in data:
        if len(row) == 0:
            continue
        row_data = [i for i in row.split("\t") if len(i) != 0]
        if len(row_data) != 2:
            print("WARNING incomplete row: %s" % row)
        else:
            tcid, substrate = row_data
        
        tcid2substrate[tcid.strip()] = substrate.strip()
    
    return tcid2substrate


def get_protein_description(acc="1.A.1.9"):
    import urllib.request
    from bs4 import UnicodeDammit
    
    url= 'http://www.tcdb.org/search/result.php?tc=%s' % acc
    
    page = urllib.request.urlopen(url).read().decode('Latin1')
    
    soup = BeautifulSoup(page, 'lxml')
    
    #table = soup.find("table", {"id": "result-cluster"})

    superfam2description = {}
    
    # retrieve fam list and description
    famlist = soup.findAll("td", {"id": "right-border"})
    for fam in famlist:
        txt = fam.text.strip()
        try:
            tcid = re.findall("[0-9]+\.[A-Z]\.[0-9]+\.[0-9]+", txt)[0]
            description = txt.split(":")[1].strip()
            superfam2description[tcid] = description
        except:
            pass
        
    rows = soup.findAll('tr')
    
    tcid2data = {}
    #print("Number of transporters:", len(rows))
    for n, row in enumerate(rows):
        cols = row.find_all('td')
        # BeautifulSoup(raw_html, "lxml").text
        cols = [BeautifulSoup(ele.text.strip(), "lxml").text.strip() for ele in cols]
        #print (len(cols))
        #print(cols)
        if len(cols) > 1:
            tcid = re.sub("\*", "", cols[0])
            tcid2data[tcid] = {}
            try:
                tcid2data[tcid]["name"] = cols[1].replace("\r\n", "").replace("\xa0", "")
                tcid2data[tcid]["domain"] = cols[2]
                tcid2data[tcid]["superkingdom"] = cols[3]
            except:
                print("WARNING: problem with row", cols)
                continue
    return tcid2data, superfam2description



def get_family2members(family2description, acc2tcid_list):
    
    family2members = {}
    
    
    for acc in acc2tcid_list:
        # one protein acc can be associated with multiple tcid
        tcid_list = acc2tcid_list[acc]
        
        # 9.A.7.1.3
        for tcid in tcid_list:
            data = tcid.split(".")
            if len(data) != 5:
                print("WARNING unexpected TCID fortmat: %s" % tcid)

            family = '%s.%s.%s' % (data[0],
                                    data[1],
                                    data[2])
            
            if family not in family2description:
                print("WARNING missing family: %s" % family)
            if family not in family2members:
                family2members[family] = [tcid]
            else:
                family2members[family].append(tcid)
                
    return family2members                  
    

if __name__ == '__main__':
    import argparse
    import time
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)

    args = parser.parse_args()

    # family_abreviations = 'http://tcdb.org/cgi-bin/projectv/family_abbreviations.py'
   
    # 1 class ok
    # 1.A superfamily ok
    # 1.A.1 family ok
    # 1.A.1.1 subfamily => no description 
    # 1.A.1.1.3 TCID

    id2description = {}
    
    print("%s Retrieveing family description..." % time.ctime())
    family2description = get_family2description()
    
    id2description.update(family2description)


    print("%s Retrieveing class and superfamily description..." % time.ctime())
    class_list = []
    for family in family2description:
        class_index = family.split(".")[0]
        if class_index not in class_list:
            class_list.append(class_index)
    
    for one_class in class_list:
        print(one_class)
        id2description.update(get_superfamilies(acc=one_class))

    print("%s Retrieveing spacies names from fasta..." % time.ctime())
    acc2tcid_list = get_acc2tcid_list()

    family2members = get_family2members(family2description, acc2tcid_list)
    
    
    print("%s Retrieveing transporters description, domain & superkingdom..." % time.ctime())
    tcid2data_all = {}
    for n, family in enumerate(family2members):
        if n % 100 == 0:
            print("%s %s / %s -- downloading transporter description" % (time.ctime(), n, len(family2members)))
            
        first_member = '.'.join(family2members[family][0].split(".")[0:4])
        tcid2data, superfam2description = get_protein_description(acc=first_member)
        # update names
        id2description.update(superfam2description)
        tcid2data_all.update(tcid2data)

    print("%s Retrieveing transporters substrates..." % time.ctime())
    tcid2substrate = get_tcid2substrate()
  
    acc2species, missing_list = parse_tcdb_fasta()
    
    print("%s Try to search species name on the UNIPROTKB" % time.ctime())
    
    annot_dico = blastswiss2sqltable.get_swissprot_annotation(missing_list)
    # update the dictionnary
    for acc in annot_dico:
        acc2species[acc] = annot_dico[acc][-1]    

    t = TCDB(args.db_name)
    t.insert_classification(id2description)

    t.insert_transporters(acc2tcid_list,
                          tcid2data_all,
                          tcid2substrate)
    
    t.import_annot(acc2tcid_list,
                   acc2species)
   
