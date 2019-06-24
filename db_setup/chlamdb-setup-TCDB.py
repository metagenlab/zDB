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
    def __init__(self):

        sqlpsw = os.environ['SQLPSW']
        self.conn = MySQLdb.connect(host="localhost",
                                    user="root",
                                    passwd=sqlpsw)
        self.cursor = self.conn.cursor()
        sql = 'create database if not exists transporters'
        self.cursor.execute(sql,)
        sql = 'use transporters'
        self.cursor.execute(sql,)
        self.conn.set_character_set('utf8mb4')
        self.cursor.execute('SET NAMES utf8mb4;')
        self.cursor.execute('SET CHARACTER SET utf8mb4;')
        self.cursor.execute('SET character_set_connection=utf8mb4;')
        self.cursor.execute('SET character_set_database=utf8mb4;')
        self.cursor.execute('SET character_set_server=utf8mb4;')
        self.conn.commit()

        sql = 'CREATE table if not exists transporter_table (transporter_id INT primary key, tc1 INT, superfamily INT,' \
              ' family INT, subfamily INT, index tc1 (tc1), index superfamily (superfamily),' \
              ' index family (family), index subfamily (subfamily))'
        self.cursor.execute(sql,)
        self.conn.commit()

        sql2 = 'CREATE table if not EXISTS tc_table (tc_id INT primary key AUTO_INCREMENT, tc_name varchar(400), description TEXT)'
        self.cursor.execute(sql2,)
        self.conn.commit()

        sql = 'CREATE TABLE if not exists uniprot_table (uniprot_id INT primary key AUTO_INCREMENT, uniprot_accession varchar(400)' \
              ' ,substrate TEXT, taxon_id INT, ' \
              ' uniprot_description TEXT, tcdb_description TEXT, organism TEXT, uniprot_gene TEXT, uniprot_annotation_score INT)'
        self.cursor.execute(sql,)
        self.conn.commit()

    def add_one_tc_id(self, tc_id):

        tc_db_id = self.check_if_tc_accession_already_in_db(tc_id)
        if not tc_db_id:
            print ('%s not into database!' % tc_id)
            description = re.sub('"', '', tcdb_utils.accession2family(tc_id))
            sql = 'insert into tc_table(tc_name, description) values (%s, %s)'
            print(sql % (tc_id, description))
            self.cursor.execute(sql, (tc_id, description))
            self.conn.commit()
            sql3 = 'select tc_id from tc_table where tc_name="%s"' % tc_id
            self.cursor.execute(sql3,)
            return self.cursor.fetchall()[0][0]
        else:
            return tc_db_id

    def insert_new_tc_path(self, tc_name):

        separated_ids = tc_name.split('.')
        if len(separated_ids) != 5:
            print ('invalid id!')
            import sys
            sys.exit()
        else:
            id1 = separated_ids[0]
            id_superfamily = '.'.join(separated_ids[0:2])
            id_family = '.'.join(separated_ids[0:3])
            id_subfamily = '.'.join(separated_ids[0:4])
            id_complete = tc_name

            id_id1 = self.add_one_tc_id(id1)
            id_superfamily = self.add_one_tc_id(id_superfamily)
            id_family = self.add_one_tc_id(id_family)
            id_subfamily = self.add_one_tc_id(id_subfamily)
            id_complete = self.add_one_tc_id(id_complete)
        sql = 'select transporter_id from transporter_table where transporter_id=%s' % id_complete
        try:
            self.cursor.execute(sql,)
            id = self.cursor.fetchall()[0][0]
            return id
        except:
            sql = 'insert into transporter_table (transporter_id, tc1,superfamily, family, subfamily) values ("%s",' \
                  ' "%s","%s","%s","%s")' % (id_complete, id_id1, id_superfamily, id_family, id_subfamily)
            self.cursor.execute(sql,)
            self.conn.commit()
            return id_complete

    def check_if_tc_accession_already_in_db(self, tc_name):

        sql3 = 'select tc_id from tc_table where tc_name="%s"' % tc_name
        try:
            self.cursor.execute(sql3,)
            tc_id = self.cursor.fetchall()[0][0]
            return tc_id
        except IndexError:
            return False

    def check_if_uniprot_accession_already_in_db(self, tc_name):
        sql3 = 'select uniprot_id from uniprot_table where uniprot_accession="%s"' % tc_name
        try:
            self.cursor.execute(sql3,)
            uniprot_id = self.cursor.fetchall()[0][0]
            return uniprot_id
        except IndexError:
            return False

    def insert_uniprot_in_db(self, uniprot_accession, tc_id, tcdb_description):

        sql2 = 'select uniprot_id from uniprot_table where uniprot_accession="%s"' % (uniprot_accession)

        try:
            self.cursor.execute(sql2,)
            return self.cursor.fetchall()[0][0]
        except:
            print ('%s not in database' % uniprot_accession)
            annot_dico = blastswiss2sqltable.get_swissprot_annotation([uniprot_accession])
            try:
                # try with uniprot
                taxon_id = annot_dico[uniprot_accession][0]
                uniprot_description = annot_dico[uniprot_accession][2]
                annot_score = annot_dico[uniprot_accession][1].split(' ')[0]
                uniprot_gene = annot_dico[uniprot_accession][3]
                uniprot_organism = annot_dico[uniprot_accession][4]
            except:
                try:
                    # try with NCBI (removed uniprot entries)
                    taxon_id = accession2taxon_id.accession2taxon_id([uniprot_accession])[uniprot_accession]
                    data = accession2taxon_id.accession2description([uniprot_accession])
                    uniprot_description = data[uniprot_accession]['description']
                    uniprot_organism = data[uniprot_accession]['source']
                    annot_score = 0
                    uniprot_gene = '-'
                except:
                    try:
                        # try with tcdb itself
                        data = tcdb_utils.accession2species_and_product(uniprot_accession, tc_id)
                        taxon_id = data[0].split('[')[1].split(']')[0]
                        uniprot_description = data[1]
                        uniprot_organism = data[0].split('[')[0]
                        annot_score = 0
                        uniprot_gene = '-'
                    except:
                        # try to search with genus and species id
                        try:
                            data = tcdb_utils.accession2species_and_product(uniprot_accession, tc_id)
                            uniprot_description = data[1]
                            uniprot_organism = data[0].split('[')[0]
                            genus_and_species = ' '.join(uniprot_organism.split(' ')[0:2])
                            taxon_id = accession2taxon_id.species_name2taxon_id(genus_and_species)
                        # if nothing works, set as "-"
                        except:
                            uniprot_description = "-"
                            uniprot_organism = "-"
                            genus_and_species = "-"
                            taxon_id = None

                        if not taxon_id:
                            # set cellular organism
                            taxon_id = 131567
                        annot_score = 0
                        uniprot_gene = '-'

            uniprot_description = re.sub('"', '', uniprot_description)
            substrate = tcdb_utils.accession2substrate(uniprot_accession, tc_id)
            sql = 'insert into uniprot_table (uniprot_accession,taxon_id, substrate, uniprot_description, tcdb_description, organism, ' \
                  ' uniprot_gene, uniprot_annotation_score) values ("%s", %s, "%s", "%s","%s","%s","%s",%s)' % (uniprot_accession,
                                                                                                           taxon_id,
                                                                                                           substrate,
                                                                                                           uniprot_description,
                                                                                                           tcdb_description,
                                                                                                           uniprot_organism,
                                                                                                           uniprot_gene,
                                                                                                           annot_score)
            self.cursor.execute(sql,)
            self.conn.commit()
            sql2 = 'select uniprot_id from uniprot_table where uniprot_accession="%s"' % (uniprot_accession)
            self.cursor.execute(sql2,)
            return self.cursor.fetchall()[0][0]

    def import_annot(self,
                     update=True,
                     replace=True,
                     db_fasta=False):

        '''

        Tables
        TC avec description, family id, family et substrats
        TCDB hit (accession)
        hits_chlamydia_04_16 avec details du hit (score,...)

        '''
        if not db_fasta:
            print("Download fasta db")
            handle = urllib.request.urlopen("http://www.tcdb.org/public/tcdb")
            fasta_st = StringIO(handle.read().decode("UTF-8"))
            db_accession2record = SeqIO.to_dict(SeqIO.parse(fasta_st, 'fasta'))
        else:
            print("Local fasta db")
            db_accession2record = SeqIO.to_dict(SeqIO.parse(db_fasta, 'fasta'))
        uniprot_nr_list = []
        for n, accession in enumerate(db_accession2record):
            print(n, accession)
            # gnl|TC-DB|1001796365|4.F.1.1.5
            header_data = accession.split("|")
            if len(header_data) == 4:
                hit_tcid = header_data[3]
            else:
                # malformed header, search accession in description
                hit_tcid = re.findall("[0-9]+\.[A-Z]\.[0-9]+\.[0-9]+\.[0-9]+", db_accession2record[accession].description)[0]
            hit_uniprot_accession = header_data[2]

            if hit_uniprot_accession not in uniprot_nr_list:
                uniprot_nr_list.append(hit_uniprot_accession)

            transporter_id = self.insert_new_tc_path(hit_tcid)
            uniprot_id = self.check_if_uniprot_accession_already_in_db(hit_uniprot_accession)
            if not uniprot_id:
                uniprot_id = self.insert_uniprot_in_db(hit_uniprot_accession, hit_tcid, db_accession2record[accession].description)

        self.conn.commit()
        print("number of nr uniprot accessions:", len(uniprot_nr_list))

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-u", '--update', action='store_true', help="update DB (add missing entries)")
    parser.add_argument("-r", '--replace', action='store_true', help="replace existing tables")
    parser.add_argument("-f", '--fasta_file', help="Fasta file (download from TCDB.org if not provided)")

    args = parser.parse_args()

    t = TCDB()
    t.import_annot(update=args.update, replace=args.replace, db_fasta=args.fasta_file)
