#!/usr/bin/env python

from chlamdb.biosqldb import manipulate_biosqldb
import sqlite3 
from Bio import SeqIO
from Bio import Medline
from Bio import Entrez
Entrez.email = 'trestan.pillone@chuv.chl'

class StringPMID():
    def __init__(self,
                 string_sqlite,
                 blast_results,
                 db_name,
                 hash2locus_list,
                 query_fasta_file,
                 db_fasta_file):

        self.string_conn = sqlite3.connect(string_sqlite)
        self.string_cursor = self.string_conn.cursor()

        self.server, self.db = manipulate_biosqldb.load_db(db_name)
        self.blast_results = blast_results
        self.hash2locus_list = hash2locus_list
        self.db_name = db_name
        self.query2len = self.parse_fasta(query_fasta_file)
        self.hit2len = self.parse_fasta(db_fasta_file)

        self.missing_pmid_data =[] 

        # create tables if not exists
        self.create_tables()

        
        # retrieve data from string db
        sql = 'select pmid,publication_date,publication_source,linkout_url,authors,title from publications;'
        
        sql_species = 'select species_id,compact_name from species'
        self.pmid2article_data = manipulate_biosqldb.to_dict(self.string_cursor.execute(sql,).fetchall())
        self.species_id2species_name = manipulate_biosqldb.to_dict(self.string_cursor.execute(sql_species).fetchall())
                
        # retrieve protein hash from biosqldb
        sql = f'select distinct hash from string.seqfeature_id2string_protein_mapping t1 inner join annotation.hash2seqfeature_id_{db_name} t2 on t1.seqfeature_id=t2.seqfeature_id;'
        self.hash_in_db =set([i[0] for i in self.server.adaptor.execute_and_fetchall(sql,)])
        
        # check if some data were already inserted
        sql = 'select pmid from string.pmid2data_stringdb'
        self.pmid_in_db = [i[0] for i in self.server.adaptor.execute_and_fetchall(sql,)]
        
        sql = 'select accession,id from string.string_protein_entry'
        self.string_proteins2string_protein_id = manipulate_biosqldb.to_dict(self.server.adaptor.execute_and_fetchall(sql,))

        

    def create_tables(self):
        sql1 = 'create table if not exists string.seqfeature_id2string_protein_mapping (seqfeature_id INT, string_protein_id INT, query_cov FLOAT, hit_cov FLOAT, identity FLOAT, evalue FLOAT, score FLOAT)'
        sql2 = 'create table if not exists string.string_protein2pmid (string_protein_id INT, PMID BIGINT)'
        sql3 = 'create table if not exists string.pmid2data_stringdb (pmid BIGINT, title TEXT, journal TEXT, year varchar(200), linkout_url TEXT, authors TEXT)'
        sql4 = 'create table if not exists string.string_protein_entry (id INT AUTO_INCREMENT PRIMARY KEY, accession varchar(200), organism TEXT, description TEXT, preferred_name varchar(400))'
        self.server.adaptor.execute(sql1)
        self.server.adaptor.execute(sql2)
        self.server.adaptor.execute(sql3)
        self.server.adaptor.execute(sql4)

    def parse_fasta(self, fasta_file):
        accession2len ={}
        records = SeqIO.parse(fasta_file, "fasta")
        for record in records:
            accession2len[record.id] = len(record.seq) 
        return accession2len


    def load_string_data(self):

        db_name = self.db_name
                
        sql = f'select locus_tag, seqfeature_id from custom_tables.locus2seqfeature_id_{db_name}'
        locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(self.server.adaptor.execute_and_fetchall(sql,))

        check_if_new = True
        for blast_file in self.blast_results:
            with open(blast_file, "r") as f:
                print("Reading blast file...")
                rows = [i.rstrip().split("\t") for i in f]
                print("OK")
                count = 1
                for n, row in enumerate(rows):
                    # skip already inserted entries
                    if check_if_new:
                        if row[0] in self.hash_in_db:
                            if n % 100000 == 0:
                                print(n)
                            continue
                        else:
                            # stop checking if the entry exists
                            self.hash_in_db = [] 
                            check_if_new = False

                    if n % 1000 == 0:
                        print(n)
                        self.server.adaptor.commit()
                    if n != 0:
                       # print(row[0], rows[n-1][0])
                        if row[0] == rows[n-1][0]:
                            count += 1
                        else:
                            print("%s\t%s\t%s" % (n, rows[n-1][0], count))
                            count = 1
                    # only keep the top 20 hits for each protein
                    if count > 20:
                        continue
                    query_hash = row[0]
                    hit_accession = row[1]
                    evalue = row[10]
                    score = row[11]
                    identity = row[2]
                    query_start = int(row[6])
                    query_end = int(row[7])
                    hit_start = int(row[8])
                    hit_end =  int(row[9])
                    query_cov = round((query_end-query_start + 1) / self.query2len[row[0]]*100, 2)
                    hit_cov = round((hit_end-hit_start + 1)/self.hit2len[row[1]]*100, 2)

                    # check if hit id already into db 
                    try:
                        string_protein_id = self.string_proteins2string_protein_id[hit_accession]
                    except KeyError:
                        string_protein_id = False
                   
                    # if not create a new entry with all associated PMID 
                    if not string_protein_id:
                        # new string entry
                        
                        sql_protein_data = f'select protein_external_id,annotation,preferred_name,pmid from proteins t1 inner join protein_id2pmid t2 on t1.protein_id=t2.protein_id where protein_external_id="{hit_accession}";'
                    
                        protein_data = list(list(i) for i in self.string_cursor.execute(sql_protein_data,).fetchall())
                        nr_pmid_list = [str(i[3]) for i in protein_data]
                        
                        protein_external_id, annotation, preferred_name, pmid = protein_data[0]
                        species_id = protein_external_id.split(".")[0]
                        species_name = self.species_id2species_name[species_id]
                        
                        sql = f'insert into string.string_protein_entry (accession, organism, description, preferred_name) values (%s, %s, %s, %s)'
                        self.server.adaptor.execute(sql, (protein_external_id, species_name, annotation, preferred_name))
                        
                        # add new id to dictionnary
                        string_protein_id = self.server.adaptor.last_id("string.string_protein_entry")
                        self.string_proteins2string_protein_id[hit_accession] = string_protein_id
                        
                        for pmid in nr_pmid_list:
                            sql = f'insert into string.string_protein2pmid (string_protein_id, pmid) values ({string_protein_id}, {pmid})'
                            self.server.adaptor.execute(sql,)
                            # insert data if not already present in db 
                            if pmid not in self.pmid_in_db:
                                article_data = self.pmid2article_data[pmid]
                                publication_date, publication_source, linkout_url, authors, title = article_data
                                sql2 = f'insert into string.pmid2data_stringdb (pmid, title, journal, year, linkout_url, authors) values (%s, %s, %s, %s, %s, %s)'
                                self.server.adaptor.execute(sql2, (pmid, title, publication_source, publication_date, linkout_url, authors))
                                self.pmid_in_db.append(pmid)

                    # insert blast result data
                    for locus_tag in self.hash2locus_list[query_hash]:
                        seqfeature_id = locus_tag2seqfeature_id[locus_tag] 
                        # insert Match 
                        sql = f'insert into string.seqfeature_id2string_protein_mapping (seqfeature_id, string_protein_id, query_cov, hit_cov, identity, evalue, score)' \
                              f' values ({seqfeature_id}, {string_protein_id}, {query_cov}, {hit_cov}, {identity}, {evalue} ,{score})'
                        self.server.adaptor.execute(sql,)
        self.server.adaptor.commit()  


if __name__ == '__main__':
    import argparse
    import chlamdb_setup_utils
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--blast_output', type=str, help="BLAST output file(s)", required=True, nargs="+")
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)
    parser.add_argument("-qf", '--query_fasta', type=str, help="db name", required=True)
    parser.add_argument("-df", '--db_fasta', type=str, help="db name", required=True)
    parser.add_argument("-p", '--string_db_name', type=str, help="string sqlite3 db_name", required=True)
    parser.add_argument("-c", '--corresp_table',  type=str, help="hash to locus correspondance table")



    args = parser.parse_args()
    
    hash2locus_list = chlamdb_setup_utils.get_hash2locus_list(args.corresp_table)

    pb = StringPMID(args.string_db_name,
                    args.blast_output,
                    args.db_name,
                    hash2locus_list,
                    args.query_fasta,
                    args.db_fasta)
    pb.load_string_data()
