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

    def add_orthogroups_to_seq(self, hsh_locus_to_group, hsh_locus_to_feature_id, orthogroup_term_id):
        rank = 1
        for locus, group_id in hsh_locus_to_group.items():
            feature_id = hsh_locus_to_feature_id[locus]
            insert = (
                f"INSERT INTO seqfeature_qualifier_value "
                f"VALUES ({feature_id}, {orthogroup_term_id}, {rank}, {quote(group_id)})"
            )
            self.server.adaptor.execute(insert)

    def get_hsh_locus_to_seqfeature_id(self):
        query = (
            "SELECT t2.value, t1.seqfeature_id "
            "FROM seqfeature as t1 "
            "INNER JOIN seqfeature_qualifier_value AS t2 ON t1.seqfeature_id=t2.seqfeature_id "
            "INNER JOIN term as t3 on t3.term_id=t2.term_id AND t3.name=\"locus_tag\""
        )
        results = self.server.adaptor.execute_and_fetchall(query,)
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = int(line[1])
        return hsh_results

    def get_locus_to_og(self):
        query = (
            f"SELECT locus.value, og.value"
            f"FROM seqfeature_qualifier_value as og "
            f"INNER JOIN term AS term_og ON term_og.term_id=og.term_id AND term_og.name=\"orthogroup\" "
            f"INNER JOIN seqfeature_qualifier_value as locus ON locus.seqfeature_id=og.seqfeature_id "
            f"INNER JOIN term AS term_locus ON term_locus.term_id=locus.term_id "
            f"  AND term_locus.name=\"locus_tag\";"
        )
        query_results = self.server.adaptor.execute_and_fetchall(query,)
        locus_to_og = {}
        for locus, og in query_results:
            locus_to_og[locus] = int(og)
        return locus_to_og

    def create_new_og_matrix(self):
        query = (
            f"CREATE TABLE orthology_identity ( "
            f"orthogroup INT, "
            f"id_1 INT, id_2 INT, identity FLOAT(5), "
            f"PRIMARY KEY(orthogroup, id_1, id_2), "
            f"FOREIGN KEY(orthogroup) REFERENCES orthology_orthogroup(orthogroup_id)"
            f"FOREIGN KEY(id_1) REFERENCES seqfeature(seqfeature_id)"
            f"FOREIGN KEY(id_2) REFERENCES seqfeature(seqfeature_id))"
        )
        self.server.adaptor.execute(query, )

    # For each taxon, get the number of protein present in each of the existing
    # ortholog group
    # Returns a table of maps, organized as follows:
    #  group 0: hsh taxid -> number of proteins of taxid in group 0
    #  group 1: hsh taxid -> number of proteins of taxid in group 1
    #  ... 
    #  group N: hsh taxid -> number of proteins of taxid in group N
    def get_orthogroup_count_table(self):
        query = (
            "SELECT value, taxon_id, count(*) "
            "FROM seqfeature_qualifier_value as ortho_table "
            "INNER JOIN term ON ortho_table.term_id = term.term_id AND term.name=\"orthogroup\" "
            "INNER JOIN seqfeature AS feature ON feature.seqfeature_id=ortho_table.seqfeature_id "
            "INNER JOIN bioentry AS entry ON entry.bioentry_id=feature.bioentry_id "
            "GROUP BY value, taxon_id; "
        )
        results = self.server.adaptor.execute_and_fetchall(query,)
        size = max([int(group) for group, taxid, cnt in results])
        arr_group = [None]*(size+1)
        for line in results:
            group, taxid, cnt = int(line[0]), line[1], int(line[2])
            hsh_taxid_to_count = arr_group[group]
            if hsh_taxid_to_count==None:
                hsh_taxid_to_count = {}
                arr_group[group] = hsh_taxid_to_count
            hsh_taxid_to_count[taxid] = cnt
        return arr_group

    def load_og_matrix(self, matrix):
        for orthogroup_id, id_1, id_2, identity in matrix:
            query = (
                f"INSERT INTO orthology_identity "
                f"VALUES ({orthogroup_id}, {id_1}, {id_2}, {identity});"
            )
            print(query)
            self.server.adaptor.execute(query,)

    def load_og_averages(self, averages):
        for og, average in averages:
            sql = (
                f"UPDATE orthology_orthogroup "
                f"SET average_identity = {average} "
                f"WHERE orthogroup_id = {og};"
            )
            self.server.adaptor.execute(sql,)

    def create_orthology_table(self, arr_cnt_tables):
        sql = (
            f"CREATE TABLE IF NOT EXISTS orthology_orthogroup( "
            f"orthogroup_id INTEGER, orthogroup_size INT, n_genomes INT, average_identity FLOAT(5), "
            f"PRIMARY KEY(orthogroup_id)) "
        )
        self.server.adaptor.execute(sql,)

        for group, hsh_cnt_tables in enumerate(arr_cnt_tables):
            size = 0
            n_genomes = 0
            for taxon_id, cnt in hsh_cnt_tables.items():
                size += cnt
                if cnt>0:
                    n_genomes += 1
            sql = (
                f"INSERT INTO orthology_orthogroup VALUES("
                f"{group}, {size}, {n_genomes}, NULL"
                f")"
            )
            self.server.adaptor.execute(sql,)

    # NOTE: may be more efficient to create indices on a combination of keys, depending
    # on how the table is used.
    def create_og_matrix_indices(self):
        sql_1 = "CREATE INDEX oio ON orthology_identity(orthogroup);"
        sql_2 = "CREATE INDEX oiid_1 ON orthology_identity(id_1);"
        sql_3 = "CREATE INDEX oiid_2 ON orthology_identity(id_2);"
        self.server.adaptor.execute(sql_1,)
        self.server.adaptor.execute(sql_2,)
        self.server.adaptor.execute(sql_3,)

    def load_orthology_table(self, arr_count_tables):
        create_table_sql = (
            f"CREATE TABLE comparative_tables_orthology (orthogroup INT, "
            f" taxid INT, count INT"
            f", PRIMARY KEY(orthogroup, taxid)"
            f", FOREIGN KEY(orthogroup) REFERENCES orthology_orthogroup(orthogroup_id)"
            f", FOREIGN KEY(taxid) REFERENCES taxon(taxon_id)"
            f")"
        )
        self.server.adaptor.execute(create_table_sql)

        for group, hsh_taxid_to_count in enumerate(arr_count_tables):
            for taxid, count in hsh_taxid_to_count.items():
                sql = (
                    f"INSERT INTO comparative_tables_orthology "
                    f"VALUES ({group}, {taxid}, {count})" 
                )
                self.server.adaptor.execute(sql,)

    def taxon_ids(self):
        query = "SELECT taxon_id FROM taxon";
        result = self.server.adaptor.execute_and_fetchall(query,)
        arr_taxon_ids = []
        for line in result:
            arr_taxon_ids.append(int(line[0]))
        return arr_taxon_ids

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

    def setup_orthology_table(self):
        sql1 = 'SELECT ontology_id FROM ontology WHERE name="SeqFeature Keys"'
        ontology_id = self.server.adaptor.execute_and_fetchall(sql1)[0][0]
        sql2 = f"INSERT INTO term (name, ontology_id) VALUES (\"orthogroup\", {ontology_id});"
        self.server.adaptor.execute(sql2)
        return self.server.adaptor.cursor.lastrowid

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
