import os
import sys

from BioSQL import BioSeqDatabase
from Bio import SeqIO


# This file defines a class DB, that encapsulates all the SQL requests
# necessary to create the chlamdb database.
# In the future, the goal is to import all database queries needed by the
# chlamdb website as methods of this class.
#
# This improves code readability by removing SQL queries from the main python
# code and more importantly, it would allow to change of database without having
# to modify chlamdb's code.


# to litteral
# encases the string into quotes
def quote(v):
    return f"\"{v}\""

class DB:

    def __init__(self, server, db_name):
        self.server = server
        self.db_name = db_name

    def create_indices_on_cds(self):
        sql_index1 = 'create index ftgcga on feature_tables_genomes_cds(genome_accession)'    
        sql_index2 = 'create index ftgctx on feature_tables_genomes_cds(taxon_id)'
        sql_index3 = 'create index ftgrga on feature_tables_genomes_rrna(taxon_id)'
        sql_index4 = 'create index ftgrtx on feature_tables_genomes_rrna(genome_accession)'    
        sql_index5 = 'create index ftcain on feature_tables_cds_accessions(id_name)'
        sql_index6 = 'create index ftcait on feature_tables_cds_accessions(id_type)'
        
        self.server.adaptor.execute(sql_index1,)
        self.server.adaptor.execute(sql_index2,)
        self.server.adaptor.execute(sql_index3,)
        self.server.adaptor.execute(sql_index4,)
        self.server.adaptor.execute(sql_index5,)
        self.server.adaptor.execute(sql_index6,)

    def create_cds_tables(self):
        sql_cds = 'CREATE table IF NOT EXISTS feature_tables_genomes_cds' \
           ' (prot_primary_id INT AUTO_INCREMENT PRIMARY KEY,' \
           ' taxon_id INT,' \
           ' genome_accession VARCHAR(40),' \
           ' start INT,' \
           ' end INT,' \
           ' strand INT,' \
           ' gene varchar(20),' \
           ' product TEXT,' \
           ' translation TEXT)'
        sql_rrna = 'CREATE table IF NOT EXISTS feature_tables_genomes_rrna' \
            ' (rrna_primary_id INT AUTO_INCREMENT PRIMARY KEY,' \
            ' taxon_id INT,' \
            ' genome_accession VARCHAR(40),' \
            ' start INT,' \
            ' end INT,' \
            ' strand INT,' \
            ' product TEXT)'
        sql_synonyms = 'CREATE table IF NOT EXISTS feature_tables_cds_accessions' \
            ' (prot_primary_id INT,' \
            ' id_type varchar(40),' \
            ' id_name varchar(40))' \
        #' FOREIGN KEY (prot_primary_id) REFERENCES feature_tables_genomes_cds(prot_primary_id))'
        self.server.adaptor.execute(sql_cds,)
        self.server.adaptor.execute(sql_rrna,)
        self.server.adaptor.execute(sql_synonyms,)

    def get_taxid_from_accession(self, accession):
        sql = (
            f"SELECT taxon_id "
            f"FROM bioentry t1 INNER JOIN biodatabase t2 "
            f"ON t1.biodatabase_id = t2.biodatabase_id "
            f"WHERE t2.name={quote(self.db_name)} AND t1.accession={quote(accession)}"
        )
        return self.server.adaptor.execute_and_fetchall(sql,)[0][0]

    def insert_rRNA_in_feature_table(self, args):
        taxon_id = args["taxon_id"]
        accession = args["accession"]
        start = args["start"]
        end = args["end"]
        strand = args["strand"]
        product = args["product"]
        sql = (
            f"INSERT INTO feature_tables_genomes_rrna"
            f"(taxon_id, genome_accession, start, end, strand, product)"
            f"VALUES"
            f"({taxon_id}, {quote(accession)}, {start}, {end}, {strand}, {quote(product)})"
        )
        print(sql)
        self.server.adaptor.execute(sql,)

    def insert_cds(self, args):
        taxon_id = args["taxon_id"]
        accession = args["accession"]
        start = args["start"]
        end = args["end"]
        strand = args["strand"]
        gene = args["gene"]
        product = args["product"]
        translation = args["translation"]
        old_locus_tag = args["old_locus_tag"]
        taxon_id = args["taxon_id"]

        sql1 = (
            f"INSERT INTO feature_tables_genomes_cds"
            f"(taxon_id, genome_accession, start, end, strand, gene, product, translation)"
            f"VALUES ({taxon_id}, {quote(accession)}, {start}, {end}, {strand},"
            f"        {quote(gene)}, {quote(product)}, {quote(translation)})"
        )
        self.server.adaptor.execute(sql1,)
        cds_id = self.server.adaptor.cursor.lastrowid
        sql2 = 'INSERT into feature_tables_cds_accessions(prot_primary_id, id_type, id_name) values (' \
               ' %s, "%s", "%s")'
        self.server.adaptor.execute(sql2 % (cds_id, 'protein_id', args["protein_id"]),)
        self.server.adaptor.execute(sql2 % (cds_id, 'locus_tag', args["locus_tag"]),)
        if old_locus_tag:
            self.server.adaptor.execute(sql2 % (cds_id, 'old_locus_tag', old_locus_tag),)

    def set_status_in_config_table(self, status_name, status_val):
        sql = f"update biodb_config set status={status_val} where name={quote(status_name)};"
        self.server.adaptor.execute(sql,)

    def create_biosql_database(self, args):
        self.server.new_database(self.db_name)

    # wrapper methods
    def commit(self):
        self.server.commit()

    def load_gbk_wrapper(self, records):
        self.server[self.db_name].load(records)

    # Maybe return different instance of a subclass depending on the type
    # of database? Would allow to avoid code duplication if several database
    # types are to be included.
    def load_db(params):
        sqlpsw = os.environ['SQLPSW']
        db_type = params["chlamdb.db_type"]
        db_name = params["chlamdb.db_name"]

        if db_type != "sqlite":
            server = BioSeqDatabase.open_database(driver="MySQLdb", 
                                                  user="root",
                                                  passwd = sqlpsw, 
                                                  host = "127.0.0.1", 
                                                  db=db_name, 
                                                  charset='utf8',
                                                  use_unicode=True)
        else:
            server = BioSeqDatabase.open_database(driver="sqlite3", 
                                                  user="root",
                                                  passwd = sqlpsw, 
                                                  host = "127.0.0.1",
                                                  db=f"{db_name}")
        return DB(server, db_name)
