import os
import sys

from BioSQL import BioSeqDatabase
from Bio import SeqIO

import sqlite3


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

    def create_hash_table(self, to_load):
        # 16 is for the length of the CRC hash converted to string
        # 4 is for the potential suffix if two hashes are identical
        # -> for now a "-n" is added at the end of the hash. Adding
        # 4 of length accounts for the "-" and a pessimistic 3 digits in case
        # there are lots of collisions.
        #
        # TODO: to optimize disk space usage, may use the max length of the hashes
        # instead, as collisions are unlikely.
        hash_len = len("CRC-") + 16 + 4
        sql = (
            f"CREATE TABLE sequence_hash_dictionnary (seqfeature_id, hash VARCHAR({hash_len})," 
            " sequence TEXT, "
            " PRIMARY KEY(seqfeature_id), "
            " FOREIGN KEY(seqfeature_id) REFERENCES seqfeature(seqfeature_id));"
        )
        self.server.adaptor.execute(sql,)
        self.load_data_into_table("sequence_hash_dictionnary", to_load)

        sql = "CREATE INDEX shd_id on sequence_hash_dictionnary (seqfeature_id)"
        self.server.adaptor.execute(sql)
        sql = "CREATE INDEX shd_hsh on sequence_hash_dictionnary (hash)"
        self.server.adaptor.execute(sql)

    # Returns a hash that maps the accesion the orthogroup it was linked to
    # Note: this function already filters the top n hits
    def get_non_PVC_refseq_matches(self, params):
        pvc = params["refseq_diamond_BBH_phylogeny_phylum_filter"]
        max_hit_count = params["refseq_diamond_BBH_phylogeny_top_n_hits"]
        pvc_fmt_string = ",".join([f"\"{phylum}\"" for phylum in pvc])
        query = (
            "SELECT refseq_hit.qseqid, refseq_hit_id.accession, feature.value, refseq_hit.hit_count "
            "FROM diamond_refseq AS refseq_hit "
            "INNER JOIN diamond_refseq_match_id AS refseq_hit_id "
            "   ON refseq_hit.sseqid=refseq_hit_id.match_id "
            " INNER JOIN seqfeature_qualifier_value as feature "
            "   ON feature.seqfeature_id=refseq_hit.qseqid "
            "INNER JOIN term as t ON t.term_id = feature.term_id AND t.name=\"orthogroup\" "
            "INNER JOIN refseq_hits_taxonomy AS taxo ON taxo.taxid=refseq_hit_id.taxid "
            "INNER JOIN taxonomy_mapping taxo_name "
            f"   ON taxo.pylum=taxo_name.taxid AND taxo_name.value NOT IN ({pvc_fmt_string}) "
            "ORDER BY refseq_hit.qseqid ASC, hit_count ASC;"
        )
        results = self.server.adaptor.execute_and_fetchall(query, )

        # The following code relies on the results being grouped by qseqid
        # and ordered by hit count by the SQL query
        hsh_results = {}
        curr_tab = []
        curr_qseq_id = None
        curr_count = 0

        # The accession not in curr_tab are not really efficient, particularly
        # since they compare strings and perform in O(n). May need to improve
        # on it lots of sequences are added for each orthogroup.
        for result in results:
            qseqid, accession, orthogroup = result[0], result[1], result[2]
            if qseqid != curr_qseq_id:
                curr_tab = hsh_results.get(orthogroup, [])
                hsh_results[orthogroup] = curr_tab
                if accession not in curr_tab:
                    curr_tab.append(accession)
                curr_count = 1
                curr_qseq_id = qseqid
            elif curr_count < max_hit_count:
                if accession not in curr_tab:
                    curr_tab.append(accession)
                curr_count += 1
        return hsh_results

    def get_all_sequences_for_orthogroup(self, orthogroup):
        query = (
            "SELECT locus.value, seq.sequence "
            "FROM seqfeature_qualifier_value AS ortho "
            "INNER JOIN term as ortho_term "
            "   ON ortho_term.term_id=ortho.term_id AND ortho_term.name=\"orthogroup\" "
            "INNER JOIN sequence_hash_dictionnary as seq ON seq.seqfeature_id = ortho.seqfeature_id "
            "INNER JOIN seqfeature_qualifier_value locus ON locus.seqfeature_id=ortho.seqfeature_id "
            "INNER JOIN term as locus_term "
            "   ON locus_term.term_id=locus.term_id AND locus_term.name=\"locus_tag\" "
            f"WHERE ortho.value={orthogroup};"
        )
        results = self.server.adaptor.execute_and_fetchall(query,)
        tab = []
        for result in results:
            tab.append([result[0], result[1]])
        return tab

    def create_diamond_refseq_match_id(self):
        query = (
            "CREATE TABLE diamond_refseq_match_id ( "
            "match_id INT, accession VARCHAR(200), taxid INT, description TEXT, length INT, "
            "PRIMARY KEY(match_id));"
        )
        self.server.adaptor.execute(query,)

    def load_diamond_refseq_match_id(self, data):
        self.load_data_into_table("diamond_refseq_match_id", data)

    def create_diamond_refseq_match_id(self):
        query = (
            "CREATE INDEX drmii ON diamond_refseq_match_id(match_id);"
        )
        self.server.adaptor.execute(query,)

    def create_refseq_hits_taxonomy(self):
        sql = (
            "CREATE TABLE IF NOT EXISTS refseq_hits_taxonomy(taxid INT, superkingdom INT, phylum INT, "
            "class INT, order_id INT, family INT, genus INT, specie INT, PRIMARY KEY(taxid));"
        )
        self.server.adaptor.execute(sql)

    def create_refseq_hits_taxonomy_indices(self):
        sql = (
            "CREATE INDEX rhti ON refseq_hits_taxonomy(taxid)"
        )
        self.server.adaptor.execute(sql)

    def create_taxonomy_mapping(self, hsh):
        sql = (
            "CREATE TABLE taxonomy_mapping(taxid INT, "
            "rank TEXT, value TEXT, PRIMARY KEY(taxid));"
        )
        self.server.adaptor.execute(sql)
        lst_values = []
        for key, (rank, value) in hsh.items():
            lst_values.append( (key, rank, value) )
        self.load_data_into_table("taxonomy_mapping", lst_values)

    def load_data_into_table(self, table, data):
        fmt_string = ", ".join(["?"] * len(data[0]))
        sql_string = f"INSERT into {table} VALUES ({fmt_string});"
        self.server.adaptor.executemany(sql_string, data)

    def load_refseq_hits_taxonomy(self, data):
        self.load_data_into_table("refseq_hits_taxonomy", data)

    def get_all_taxids(self):
        query = (
            "SELECT DISTINCT taxid from diamond_refseq_match_id;"
        )
        results = self.server.adaptor.execute_and_fetchall(query,)
        taxids = []
        for line in results:
            taxids.append(int(line[0]))
        return taxids

    hsh_taxo_key = {
            "superkingdom" : (1, 2),
            "phylum" :       (3, 4),
            "class" :        (5, 6),
            "order" :        (7, 8),
            "family" :       (9, 10),
            "genus" :        (11, 12),
            "species" :      (13, 14)
    }
    def parse_taxo(self, results):
        # taxid
        lst = [ results[0] ]
        for (idx_name, idx_taxid) in self.hsh_taxo_key.values():
            lst.append( results[idx_name] )
            lst.append( results[idx_taxid] )
        return lst

    def get_linear_taxonomy(self, args, taxids):
        conn_refseq = sqlite3.connect(args["databases_dir"] + "/ncbi-taxonomy/linear_taxonomy.db")
        cursor = conn_refseq.cursor()

        query_string = ",".join([str(i) for i in taxids])
        query = (
            "SELECT tax_id, `superkingdom`, superkingdom_taxid, "
            " `phylum`, phylum_taxid, `class`, class_taxid, "
            " `order`, order_taxid, `family`, family_taxid, "
            " `genus`, genus_taxid, `species`, species_taxid "
            f"FROM ncbi_taxonomy WHERE tax_id IN ({query_string});"
        )
        results = cursor.execute(query, ).fetchall()
        lst_results = []
        for line in results:
            lst_results.append(self.parse_taxo(line))
        return lst_results

    def get_accession_to_taxid(self, accession, params):
        conn_refseq = sqlite3.connect(params["databases_dir"] + "/ncbi-taxonomy/prot_accession2taxid.db")
        cursor = conn_refseq.cursor()

        query_string = ",".join(["\"" + a + "\"" for a in accession])
        query = (
            f"SELECT accession, taxid FROM accession2taxid "
            f"WHERE accession IN ({query_string});" 
        )
        results = cursor.execute(query,).fetchall()
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = int(line[1])
        return hsh_results

    def get_accession_to_prot(self, accession, params):
        conn_refseq = sqlite3.connect(params["databases_dir"] + "/refseq/merged_refseq.db")
        cursor = conn_refseq.cursor()

        fmt_string = ",".join(["\"" + a + "\"" for a in accession])
        query = (
            f"SELECT accession, description, sequence_length FROM refseq "
            f"WHERE accession IN ({fmt_string});"
        )
        results = cursor.execute(query,).fetchall()
        hsh_results = {}
        for line in results:
            hsh_results[line[0]] = (line[1], int(line[2]))
        return hsh_results

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

    def create_refseq_hits_table(self):
        sql = (
            f"CREATE TABLE IF NOT EXISTS diamond_refseq(hit_count INT, qseqid INT, "
            f"sseqid INT, pident FLOAT, length INT, mismatch INT, "
            f"gapopen INT, qstart INT, qend INT, sstart INT, "
            f"send INT, evalue FLOAT, bitscore FLOAT, "
            f"PRIMARY KEY (hit_count, qseqid), "
            f"FOREIGN KEY(qseqid) REFERENCES seqfeature(seqfeature_id));"
        )
        self.server.adaptor.execute(sql,)

    def hash_to_seqfeature(self):
        query = (
            "SELECT seqfeature_id, hash FROM sequence_hash_dictionnary;"
        )
        results = self.server.adaptor.execute_and_fetchall(query,)
        hsh = {}
        for line in results:
            seqfeature_id = int(line[0])
            crc_hsh = line[1]
            hsh[crc_hsh] = seqfeature_id
        return hsh

    def load_refseq_hits(self, data):
        self.load_data_into_table("diamond_refseq", data)

    def setup_orthology_table(self):
        sql1 = 'SELECT ontology_id FROM ontology WHERE name="SeqFeature Keys"'
        ontology_id = self.server.adaptor.execute_and_fetchall(sql1)[0][0]
        sql2 = f"INSERT INTO term (name, ontology_id) VALUES (\"orthogroup\", {ontology_id});"
        self.server.adaptor.execute(sql2)
        return self.server.adaptor.cursor.lastrowid

    def load_reference_phylogeny(self, tree):
        sql = "CREATE TABLE reference_phylogeny (tree TEXT);"
        self.server.adaptor.execute(sql,)

        sql = "INSERT INTO reference_phylogeny VALUES (?);"
        self.server.adaptor.execute(sql, tree)

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
