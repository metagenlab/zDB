#!/usr/bin/env python

from chlamdb.biosqldb import manipulate_biosqldb
import sqlite3 

class PaperBlast():
    def __init__(self,
                 paperblast_sqlite,
                 blast_results,
                 db_name,
                 hash2locus_list):

        self.paperblast_conn = sqlite3.connect(paperblast_sqlite)
        self.paperblast_cursor = self.paperblast_conn.cursor()

        self.server, self.db = manipulate_biosqldb.load_db(db_name)
        self.blast_results = blast_results
        self.hash2locus_list = hash2locus_list
        self.db_name = db_name

        sql = 'select pmId,doi,title,journal,year from (select distinct pmId,doi,title,journal,year from GenePaper) A;'
        self.pmid2article_data = manipulate_biosqldb.to_dict(self.paperblast_cursor.execute(sql,).fetchall())

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
        #print(nr_pmid_list)
        #print(gene_RIF)

    def load_paperblast_data(self):

        
        db_name = self.db_name
                
        sql = f'select locus_tag, seqfeature_id from custom_tables.locus2seqfeature_id_{db_name}'
        locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(self.server.adaptor.execute_and_fetchall(sql,))
        
        '''
        sql = f'create table if not exists custom_tables.seqfeature_id2paperblast_{db_name} (seqfeature_id INT, pdb_accession varchar(200), header TEXT, compound TEXT, source TEXT, resolution FLOAT, method TEXT, identity FLOAT, score FLOAT, evalue FLOAT)'
        self.server.adaptor.execute(sql,)
        
        sql_template = f'insert into custom_tables.seqfeature_id2paperblast_{db_name} values (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'
        '''

        for blast_file in self.blast_results:
            with open(blast_file, "r") as f:
                print("Reading blast file...")
                rows = [i.rstrip().split("\t") for i in f]
                print("OK")
                count = 1
                for n, row in enumerate(f):
                    if n !=0:
                        if row[0] == rows[n-1][0]:
                            count += 1
                        else:
                            count = 1
                    # only keep the top 20 hits for each protein
                    if count > 20:
                        continue
                    if n%1000 == 0:
                        print(n)
                    query_hash = row[0]
                    hit_accession = row[1]
                    if '::' in hit_accession:
                        hit_accession = hit_accession.split("::")[1]
                    e_value = row[10]
                    score = row[11]
                    identity = row[2]

                    self.retrieve_paperblast_data(hit_accession)
                    
                    '''
                    for locus in self.hash2locus_list[query_hash]:
                        self.server.adaptor.execute(sql_template, (locus_tag2seqfeature_id[locus],
                                                    hit_accession,
                                                    header,
                                                    compound,
                                                    source,
                                                    resolution,
                                                    method,
                                                    identity,
                                                    score,
                                                    e_value))
                    '''
         
        #sql = f'create index sfp on custom_tables.seqfeature_id2paperblast_{db_name}(seqfeature_id)'                 
        #self.server.adaptor.execute(sql,)
        #self.server.commit()


if __name__ == '__main__':
    import argparse
    import chlamdb_setup_utils
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--blast_output', type=str, help="BLAST output file(s)", required=True, nargs="+")
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)
    parser.add_argument("-p", '--paperblast_db_name', type=str, help="paperblast db_name", required=True)
    parser.add_argument("-c", '--corresp_table',  type=str, help="hash to locus correspondance table")



    args = parser.parse_args()
    
    hash2locus_list = chlamdb_setup_utils.get_hash2locus_list(args.corresp_table)

    pb = PaperBlast(args.paperblast_db_name,
                    args.blast_output,
                    args.db_name,
                    hash2locus_list)
    pb.load_paperblast_data()
