#!/usr/bin/env python

from Bio import Entrez


Entrez.email = "trestan.pillonel@unil.ch"


class Uniprot_annot():

    def __init__(self,
                 biodb,
                 uniprot_mapping,
                 uniprot_data,
                 hash_locus_mapping):

        from chlamdb.biosqldb import manipulate_biosqldb
        from Bio import SeqIO
        import http.client
        import sqlite3

        self.sqlite_db_path = '/tmp/uniprot.db'
        self.biodb = biodb
        self.uniprot_mapping = uniprot_mapping
        self.uniprot_data = uniprot_data
        self.hash_locus_mapping = hash_locus_mapping
        self.taxid2most_frequent_proteome = {}

        # create temporary sqlite3 database
        self.sqlite_conn = sqlite3.connect(self.sqlite_db_path)
        self.sqlite_cursor = self.sqlite_conn.cursor()
        self.mysql_server, self.mysql_db = manipulate_biosqldb.load_db(self.biodb)

    def get_whole_db_uniprot_crossref(self):

        import MySQLdb
        from datetime import datetime
        import time
        from chlamdb.biosqldb import manipulate_biosqldb
        import re
        import os

        from tempfile import NamedTemporaryFile

        server, db = manipulate_biosqldb.load_db(self.biodb)
        conn = server.adaptor.conn
        cursor = server.adaptor.cursor


        sql1 = 'CREATE TABLE IF NOT EXISTS custom_tables_uniprot_id2seqfeature_id (seqfeature_id INT UNIQUE, uniprot_id INT AUTO_INCREMENT,' \
               ' uniprot_accession varchar(400), uniprot_status varchar(400), annotation_score INT, proteome varchar(200), insert_date varchar(300), INDEX uniprot_id(uniprot_id))'


        sql5 = 'CREATE TABLE IF NOT EXISTS custom_tables_uniprot_annotation (seqfeature_id INT, comment_function TEXT,' \
               ' ec_number TEXT,comment_similarity TEXT,comment_catalyticactivity TEXT,comment_pathway TEXT,keywords TEXT,' \
               ' comment_subunit TEXT, gene TEXT, recommendedName_fullName TEXT, proteinExistence TEXT, ' \
               ' developmentalstage TEXT, index seqfeature_id(seqfeature_id))'

        #sql6 = 'CREATE TABLE IF NOT EXISTS uniprot_keywords_%s (seqfeature_id INT, uniprot_accession varchar(200), keyword TEXT)' % self.biodb

        cursor.execute(sql1, )
        cursor.execute(sql5, )
        #cursor.execute(sql6, )
        conn.commit()

        sql1 = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id'

        cursor.execute(sql1, )
        locus2seqfeature_id = manipulate_biosqldb.to_dict(cursor.fetchall())

        sql = 'select distinct taxid, t1.locus_tag from locus_tag2hash t1 inner join hash2uniprot_accession t2 on t1.hash=t2.hash inner join taxid2locus_tag t3 on t1.locus_tag=t3.locus_tag;'
        self.sqlite_cursor.execute(sql,)
        taxid_and_locus_list = [i for i in self.sqlite_cursor.fetchall()]

        for i, row in enumerate(taxid_and_locus_list):
            if i % 10 == 0:
                print(i)
                conn.commit()
            taxid, locus_tag = row
            try:
                most_common_proteome = self.taxid2most_frequent_proteome[taxid]
            except KeyError:
                most_common_proteome = False

            if most_common_proteome:
                sql = 'select t4.* from taxid2locus_tag t1 inner join locus_tag2hash t2 on t1.locus_tag=t2.locus_tag inner join hash2uniprot_accession t3 on t2.hash=t3.hash inner join uniprot_data t4 on t3.uniprot_accession=t4.uniprot_accession where t1.locus_tag="%s" and proteome="%s";' % (locus_tag, most_common_proteome)
                try:
                    self.sqlite_cursor.execute(sql,)
                    annot = self.sqlite_cursor.fetchall()[0]
                except IndexError:
                    # no match with most frequent proteome, keep first hit with priority to swissprot entries (ORDER BY alphabetical order)
                    sql = 'select t4.* from taxid2locus_tag t1 inner join locus_tag2hash t2 on t1.locus_tag=t2.locus_tag inner join hash2uniprot_accession t3 on t2.hash=t3.hash inner join uniprot_data t4 on t3.uniprot_accession=t4.uniprot_accession where t1.locus_tag="%s" ORDER by t3.source_db;' % (locus_tag)
                    self.sqlite_cursor.execute(sql,)
                    annot = self.sqlite_cursor.fetchall()[0]

            '''
            0    uniprot_accession
            1    uniprot_score
            2    uniprot_status
            3    proteome
            4    comment_function
            5    ec_number
            6    comment_subunit
            7    gene
            8    recommendedName_fullName
            9    proteinExistence
            10    developmentalstage
            11    comment_similarity
            12    comment_catalyticactivity
            13    comment_pathway
            14    keywords
            '''
            
            seqfeature_id = locus2seqfeature_id[locus_tag]
            comment_function = annot[4]
            ec_number = annot[5]
            comment_similarity = annot[11]
            comment_catalyticactivity = annot[12]
            comment_pathway = annot[13]
            keywords = annot[14]
            keywords_list = annot[14].split(";")
            comment_subunit = annot[6]
            gene = annot[7]
            recommendedName_fullName = annot[8]
            proteinExistence = annot[9]
            developmentalstage = annot[10]

            uniprot_accession = annot[0]
            uniprot_status = annot[2]
            proteome = annot[3]
            # '1 out of 5'
            if not isinstance(annot[1], int):
                uniprot_score = annot[1].split(' ')[0]
            else:
                uniprot_score = annot[1]

            now = datetime.now()
            str_date = "%s-%s-%s" % (now.year, now.month, now.day)

            sql = 'insert into custom_tables_uniprot_id2seqfeature_id'
            sql += '(seqfeature_id, uniprot_accession, uniprot_status, annotation_score, proteome, insert_date) ' \
                   ' values (%s, %s, %s, %s, %s, %s)'
            cursor.execute(sql, (seqfeature_id,
                               uniprot_accession,
                               uniprot_status,
                               uniprot_score,
                               proteome,
                               str_date))

            # add annotation
            sql = 'insert into custom_tables_uniprot_annotation'
            sql += '(seqfeature_id, comment_function,' \
                   ' ec_number,comment_similarity,comment_catalyticactivity,comment_pathway,keywords,' \
                   ' comment_subunit, gene, recommendedName_fullName, proteinExistence,developmentalstage) values' \
                   '  (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'

            cursor.execute(sql, (seqfeature_id,
                                 comment_function,
                                 ec_number,
                                 comment_similarity,
                                 comment_catalyticactivity,
                                 comment_pathway,
                                 keywords,
                                 comment_subunit,
                                 gene,
                                 recommendedName_fullName,
                                 proteinExistence,
                                 developmentalstage))
        conn.commit()

    def get_most_frequent_proteome(self):

        # import taxid2locus_tag
        sql1 = 'create table taxid2locus_tag (taxid INT, locus_tag varchar(200))'
        self.sqlite_cursor.execute(sql1,)
        sql2 = 'select taxon_id, locus_tag from orthology_detail'
        data = self.mysql_server.adaptor.execute_and_fetchall(sql2,)
        sql = 'insert into taxid2locus_tag values (?, ?)'
        for row in data:
            self.sqlite_cursor.execute(sql, row)
        self.sqlite_conn.commit()

        # indexes
        sqlid1 = 'create index taxid on taxid2locus_tag(taxid);'
        sqlid2 = 'create index loc on taxid2locus_tag(locus_tag);'
        self.sqlite_cursor.execute(sqlid1)
        self.sqlite_cursor.execute(sqlid2)

        # import locus_tag2hash
        sql = 'create table locus_tag2hash (locus_tag varchar(200), hash binary)'
        self.sqlite_cursor.execute(sql,)
        with open(self.hash_locus_mapping, 'r') as f:
            sql = 'insert into locus_tag2hash values (?, ?)'
            for row in f:
                data = row.rstrip().split("\t")
                self.sqlite_cursor.execute(sql, data)
        self.sqlite_conn.commit()

        # indexes
        sqlid1 = 'create index lh1 on locus_tag2hash(locus_tag);'
        sqlid2 = 'create index lh2 on locus_tag2hash(hash);'
        self.sqlite_cursor.execute(sqlid1)
        self.sqlite_cursor.execute(sqlid2)

        # import hash2uniprot_accession
        sql = 'create table hash2uniprot_accession (hash binary, uniprot_accession varchar(200), taxon_id INTEGER, description TEXT, source_db varchar(200))'
        self.sqlite_cursor.execute(sql,)
        with open(self.uniprot_mapping, 'r') as f:
            sql = 'insert into hash2uniprot_accession values (?, ?, ?, ?, ?)'
            for n, row in enumerate(f):
                if n == 0:
                    continue
                data = row.rstrip().split("\t")
                data[1] = data[1].split(".")[0]
                self.sqlite_cursor.execute(sql, data)
        # indexes
        sqlid1 = 'create index hu1 on hash2uniprot_accession(hash);'
        sqlid2 = 'create index hu2 on hash2uniprot_accession(uniprot_accession);'
        sqlid3 = 'create index hu3 on hash2uniprot_accession(taxon_id);'
        sqlid4 = 'create index hu4 on hash2uniprot_accession(source_db);'
        self.sqlite_cursor.execute(sqlid1)
        self.sqlite_cursor.execute(sqlid2)
        self.sqlite_cursor.execute(sqlid3)
        self.sqlite_cursor.execute(sqlid4)

        self.sqlite_conn.commit()

        # import uniprot2data
        sql = 'create table uniprot_data (uniprot_accession varchar(200), uniprot_score INTEGER, uniprot_status varchar(200), proteome varchar(200),comment_function TEXT, ec_number varchar(200),comment_subunit TEXT, gene varchar(200), recommendedName_fullName TEXT, proteinExistence TEXT, developmentalstage TEXT,comment_similarity TEXT,comment_catalyticactivity TEXT, comment_pathway TEXT, keywords TEXT)'
        self.sqlite_cursor.execute(sql,)
        with open(self.uniprot_data, 'r') as f:
            sql = 'insert into uniprot_data values (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)'
            for n,row in enumerate(f):
                if n == 0:
                    continue
                data = row.rstrip().split("\t")
                self.sqlite_cursor.execute(sql, data)
        self.sqlite_conn.commit()

        # indexes
        sqlid1 = 'create index up1 on uniprot_data(uniprot_accession);'
        sqlid2 = 'create index up2 on uniprot_data(proteome);'
        self.sqlite_cursor.execute(sqlid1)
        self.sqlite_cursor.execute(sqlid2)

        sql = 'select taxid,proteome, count(*) as n from taxid2locus_tag t1 inner join locus_tag2hash t2 on t1.locus_tag=t2.locus_tag inner join hash2uniprot_accession t3 on t2.hash=t3.hash inner join uniprot_data t4 on t3.uniprot_accession=t4.uniprot_accession group by taxid,proteome order by n DESC;'
        self.sqlite_cursor.execute(sql,)

        for row in self.sqlite_cursor.fetchall():
            if row[0] not in self.taxid2most_frequent_proteome and row[1] != '-':
                self.taxid2most_frequent_proteome[row[0]] = row[1]


if __name__ == '__main__':
    import argparse
    import sys
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-d",'--biodb', type=str, help="biodb_name")
    parser.add_argument("-um",'--uniprot_mapping', type=str, help="uniprot_mapping")
    parser.add_argument("-ud",'--uniprot_data', type=str, help="uniprot_data")
    parser.add_argument("-hm",'--hash_locus_mapping', type=str, help="hash_locus_mapping")

    args = parser.parse_args()

    t = Uniprot_annot(args.biodb,
                      args.uniprot_mapping,
                      args.uniprot_data,
                      args.hash_locus_mapping)
    t.get_most_frequent_proteome()
    t.get_whole_db_uniprot_crossref()
