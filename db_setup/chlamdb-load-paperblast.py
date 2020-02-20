#!/usr/bin/env python

from chlamdb.biosqldb import manipulate_biosqldb
import sqlite3 
from Bio import SeqIO
from Bio import Medline
from Bio import Entrez
Entrez.email = 'trestan.pillone@chuv.chl'

class PaperBlast():
    def __init__(self,
                 paperblast_sqlite,
                 blast_results,
                 db_name,
                 hash2locus_list,
                 query_fasta_file,
                 db_fasta_file):

        self.paperblast_conn = sqlite3.connect(paperblast_sqlite)
        self.paperblast_cursor = self.paperblast_conn.cursor()

        self.server, self.db = manipulate_biosqldb.load_db(db_name)
        self.blast_results = blast_results
        self.hash2locus_list = hash2locus_list
        self.db_name = db_name
        self.query2len = self.parse_fasta(query_fasta_file)
        self.hit2len = self.parse_fasta(db_fasta_file)

        self.missing_pmid_data =[] 

        # create tables if not exists
        self.create_tables()

        # check if some data were already inserted
        sql = 'select pmId,doi,title,journal,year from (select distinct pmId,doi,title,journal,year from GenePaper) A;'
        self.pmid2article_data = manipulate_biosqldb.to_dict(self.paperblast_cursor.execute(sql,).fetchall())
        sql = f'select distinct hash from string_seqfeature_id2paperblast t1 inner join annotation.hash2seqfeature_id_{db_name} t2 on t1.seqfeature_id=t2.seqfeature_id;'
        self.hash_in_db =set([i[0] for i in self.server.adaptor.execute_and_fetchall(sql,)])
        sql = 'select pmid from string_pmid2data_paperblast'
        self.pmid_in_db = [i[0] for i in self.server.adaptor.execute_and_fetchall(sql,)]
        sql = 'select pmid from string_pmid2gene_rif_comment'
        self.gene_rif_in_db = [i[0] for i in self.server.adaptor.execute_and_fetchall(sql,)]
        sql = 'select accession,id from string_paperblast_entry'
        self.paperblast_accession2paperblast_id = manipulate_biosqldb.to_dict(self.server.adaptor.execute_and_fetchall(sql,))
        sql_desc_gene_curated = 'select protId,organism,desc,comment from CuratedGene'
        sql_desc_gene = 'select geneId,organism,desc from Gene'
        self.paperblast_accession2desc =  manipulate_biosqldb.to_dict(self.paperblast_cursor.execute(sql_desc_gene_curated,).fetchall())
        self.paperblast_accession2desc.update(manipulate_biosqldb.to_dict(self.paperblast_cursor.execute(sql_desc_gene).fetchall()))
        

    def create_tables(self):
        sql1 = 'create table if not exists string_seqfeature_id2paperblast (seqfeature_id INT, paperblast_id INT, query_cov FLOAT, hit_cov FLOAT, identity FLOAT, evalue FLOAT, score FLOAT)'
        sql2 = 'create table if not exists string_paperblast2pmid (paperblast_id INT, PMID BIGINT)'
        sql3 = 'create table if not exists string_pmid2data_paperblast(pmid BIGINT, title TEXT, journal TEXT, year varchar(200))'
        sql4 = 'create table if not exists string_paperblast_entry (id INT AUTO_INCREMENT PRIMARY KEY, accession varchar(200), db_name varchar(200), organism TEXT, description TEXT, comment TEXT)'
        sql5 = 'create table if not exists string_pmid2gene_rif_comment (pmid BIGINT, comment TEXT)'
        self.server.adaptor.execute(sql1)
        self.server.adaptor.execute(sql2)
        self.server.adaptor.execute(sql3)
        self.server.adaptor.execute(sql4)
        self.server.adaptor.execute(sql5)


    def parse_fasta(self, fasta_file):

        accession2len ={}
        records = SeqIO.parse(fasta_file, "fasta")
        for record in records:
            accession2len[record.id] = len(record.seq) 
        return accession2len


    def retrieve_paperblast_data(self,
                                 accession):

        sql_synonymous = f'select duplicate_id from SeqToDuplicate where sequence_id="{accession}"'
        
        duplicate_list = [i[0] for i in self.paperblast_cursor.execute(sql_synonymous,).fetchall()]
        
        acc_list = [accession] + duplicate_list

        pmid_list = []
        gene_RIF = []
        for accession in acc_list:
            #print("ACC", accession)
            if '::' in accession:
                accession = accession.split('::')[1]

            sql_gene = f'select pmId from GenePaper where geneId="{accession}"'
            sql_prot = f'select pmId from CuratedPaper where protId="{accession}"'
            sql_geneRIF = f'select pmId,comment from GeneRIF where geneId="{accession}";'

            data_gene = [i[0] for i in self.paperblast_cursor.execute(sql_gene,).fetchall()]
            data_prot = [i[0] for i in self.paperblast_cursor.execute(sql_prot,).fetchall()]
            data_geneRIF = [i for i in self.paperblast_cursor.execute(sql_geneRIF,).fetchall()]

            pmid_list += data_gene
            pmid_list += data_prot
            gene_RIF += data_geneRIF

        nr_pmid_list = list(set(pmid_list))
        nr_gene_RIF = list(set(gene_RIF))

        return nr_pmid_list, gene_RIF



    def pmid2abstract_info(self, pmid_list):
        

        # make sure that pmid are strings
        pmid_list = [str(i) for i in pmid_list]

        try:
            handle = Entrez.efetch(db="pubmed", id=','.join(pmid_list), rettype="medline", retmode="text")
            records = Medline.parse(handle)
        except:
            print("FAIL:", pmid_list)
            return None

        pmid2data = {}
        for record in records:
            try:
                pmid = record["PMID"]
            except:
                print(record)
                #{'id:': ['696885 Error occurred: PMID 28696885 is a duplicate of PMID 17633143']}
                if 'duplicate' in record['id:']:
                    duplicate = record['id:'].split(' ')[0]
                    correct = record['id:'].split(' ')[-1]
                    print("removing duplicated PMID... %s --> %s" % (duplicate, correct))
                    # remove duplicate from list
                    pmid_list.remove(duplicate)
                    return self.pmid2abstract_info(pmid_list)

            pmid2data[pmid] = {}
            pmid2data[pmid]["title"] = record.get("TI", "?")
            pmid2data[pmid]["authors"] = record.get("AU", "?")
            pmid2data[pmid]["source"] = record.get("SO", "?")
            pmid2data[pmid]["abstract"] = record.get("AB", "?")
            pmid2data[pmid]["journal"] = record.get("TA", "?")
            pmid2data[pmid]["year"] = record.get("DP", "?")
            pmid2data[pmid]["pmid"] = pmid

        return pmid2data

    def add_pmid_metadata(self, pmid_list):
        data = self.pmid2abstract_info(list(set(pmid_list)))
        for pmid in pmid_list:
            try:
                title = data[pmid]["title"]
                journal = data[pmid]["journal"]
                year = data[pmid]["year"]
                sql2 = f'insert into string_pmid2data_paperblast (pmid, title, journal, year) values (%s, %s, %s, %s)'
                self.server.adaptor.execute(sql2, (pmid,title,journal, year))
                self.pmid_in_db.append(pmid)
            except KeyError:
                print("Missing pmid from PubMed:", pmid)
                            
    def load_paperblast_data(self):

        
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
                    # insert pmid and reset list 
                    if len(self.missing_pmid_data) > 70:
                        print("Searching metadata on pubmed:", self.missing_pmid_data[0:5], "...")
                        self.add_pmid_metadata(self.missing_pmid_data)
                        self.missing_pmid_data =[] 
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
                    if '::' in hit_accession:
                        db, hit_accession = hit_accession.split("::")
                    else:
                        db = '-'
                    evalue = row[10]
                    score = row[11]
                    identity = row[2]
                    query_start = int(row[6])
                    query_end = int(row[7])
                    hit_start = int(row[8])
                    hit_end =  int(row[9])
                    query_cov = round((query_end-query_start + 1) / self.query2len[row[0]]*100, 2)
                    hit_cov = round((hit_end-hit_start + 1)/self.hit2len[row[1]]*100, 2)

                    nr_pmid_list, gene_RIF = self.retrieve_paperblast_data(hit_accession)
                    gene_rif_pmid_list = [i[0] for i in gene_RIF]
                    nr_pmid_list+=gene_rif_pmid_list

                    # add RIF comments 
                    for pmid, comment in gene_RIF:
                        if pmid == "":
                            continue
                        if pmid not in self.gene_rif_in_db:
                            sql = 'insert into string_pmid2gene_rif_comment (pmid, comment) values (%s, %s)'
                            self.server.adaptor.execute(sql, (pmid, comment))
                            self.gene_rif_in_db.append(pmid)

                    # check if hit id already into db 
                    try:
                        paperblast_id = self.paperblast_accession2paperblast_id[hit_accession]
                    except KeyError:
                        paperblast_id = False
                   
                    # if not create a new entry with all associated PMID 
                    if not paperblast_id:
                        # new paperblast entry 
                        paperblast_data = self.paperblast_accession2desc[hit_accession] 
                        if len(paperblast_data) == 3:
                            organism, desc, comment = paperblast_data
                        else:
                            organism, desc = paperblast_data
                            comment = '-'
                        sql = f'insert into string_paperblast_entry (accession, db_name, organism, description, comment) values (%s, %s, %s, %s, %s)'
                        self.server.adaptor.execute(sql, (hit_accession, db, organism, desc, comment))
                        paperblast_id = self.server.adaptor.last_id("string_paperblast_entry")
                        self.paperblast_accession2paperblast_id[hit_accession] = paperblast_id
                        for pmid in nr_pmid_list:
                            # skip missing pmid 
                            if pmid == "":
                                continue
                            sql = f'insert into string_paperblast2pmid (paperblast_id, pmid) values ({paperblast_id}, {pmid})'
                            self.server.adaptor.execute(sql,)
                            # insert data if not already present in db 
                            if pmid not in self.pmid_in_db:
                                #  pmId,doi,title,journal,year
                                try:
                                    article_data = self.pmid2article_data[pmid]
                                    doi,title,journal,year = article_data
                                    sql2 = f'insert into string_pmid2data_paperblast (pmid, title, journal, year, comment) values (%s, %s, %s, %s, "-")'
                                    self.server.adaptor.execute(sql2, (pmid,title,journal, year))
                                    self.pmid_in_db.append(pmid)
                                except:
                                    article_data = False
                                if not article_data:
                                    if pmid not in self.missing_pmid_data:
                                        self.missing_pmid_data.append(pmid)

                    for locus_tag in self.hash2locus_list[query_hash]:
                        seqfeature_id = locus_tag2seqfeature_id[locus_tag] 
                        # insert Match 
                        sql = f'insert into string_seqfeature_id2paperblast (seqfeature_id, paperblast_id, query_cov, hit_cov, identity, evalue, score)' \
                              f' values ({seqfeature_id}, {paperblast_id}, {query_cov}, {hit_cov}, {identity}, {evalue} ,{score})'
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
    parser.add_argument("-p", '--paperblast_db_name', type=str, help="paperblast db_name", required=True)
    parser.add_argument("-c", '--corresp_table',  type=str, help="hash to locus correspondance table")



    args = parser.parse_args()
    
    hash2locus_list = chlamdb_setup_utils.get_hash2locus_list(args.corresp_table)

    pb = PaperBlast(args.paperblast_db_name,
                    args.blast_output,
                    args.db_name,
                    hash2locus_list,
                    args.query_fasta,
                    args.db_fasta)
    pb.load_paperblast_data()
