import os
import sys

# to litteral
# encases the string into quotes
def quote(v):
    return f"\"{to_l}\""

class DB:

    def Get_DB(db_type):
        return None

    def insert_cds(self, args):
        taxon_id = args["taxon_id"]
        accession = args["accession"]
        sql1 = (
            f"INSERT INTO feature_tables.genomes_cds(taxon_id,"
            f"genome_accession, start, end, strand, gene, product, translation)"
            f"values ({taxon_id}, {quote(accession)}, {start}, {end}, {strand}"
            f"        {quote(gene)}, {quote(product)}, {quote(translation)})"
        )
        server.adaptor.execute(sql1,)
        server.commit()


def create_data_table(kwargs):
    db_type = kwargs["chlamdb.db_type"]
    db_name = kwargs["chlamdb.db_name"]
    if db_type=="sqlite":
        import sqlite3
        conn = sqlite3.connect(db_name)
        cursor = conn.cursor()
    else:
        import os
        import MySQLdb
        sqlpsw = os.environ['SQLPSW']

        conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                    user="root", # your username
                                    passwd=sqlpsw, # your password
                                    db=db_name) # name of the data base
        cursor = conn.cursor()


    entry_list = [
        ("gbk_files", "mandatory", False),
        ("orthology_data", "mandatory", False),
        ("orthology_comparative", "mandatory", False),
        ("orthology_comparative_accession", "mandatory", False),
        ("orthology_consensus_annotation", "mandatory", False),
        ("orthogroup_alignments", "mandatory", False),
        ("old_locus_table", "mandatory", False),
        ("reference_phylogeny", "mandatory", False),
        ("taxonomy_table", "mandatory", False),
        ("genome_statistics", "mandatory", False),
        ("BLAST_database", "optional", False),
        ("gene_phylogenies", "optional", False),
        ("interpro_data", "optional", False),
        ("interpro_comparative", "optional", False),
        ("interpro_comparative_accession", "optional", False),
        ("priam_data", "optional", False),
        ("priam_comparative", "optional", False),
        ("priam_comparative_accession", "optional", False),
        ("COG_data", "optional", False),
        ("COG_comparative", "optional", False),
        ("COG_comparative_accession", "optional", False),
        ("KEGG_data", "optional", False),
        ("KEGG_comparative", "optional", False),
        ("KEGG_comparative_accession", "optional", False),
        ("pfam_comparative", "optional", False),
        ("pfam_comparative_accession", "optional", False),       
        ("TCDB_data", "optional", False),
        ("psortb_data", "optional", False),
        ("T3SS_data", "optional", False),
        ("PDB_data", "optional", False),
        ("BLAST_refseq", "optional", False),
        ("BLAST_swissprot", "optional", False),
        ("BBH_phylogenies", "optional", False),
        ("GC_statistics", "optional", False),
        ("gene_clusters", "optional", False),
        ("phylogenetic_profile", "optional", False),
        ("synonymous_table", "optional", False),
        ("interpro_taxonomy", "optional", False), # interpro taxnonomy statistics
        ("pfam_taxonomy", "optional", False), #  taxnonomy statistics
        ("COG_taxonomy", "optional", False) # COG taxnonomy statistics
    ]
    
    sql = 'create table biodb_config (name varchar(200), type varchar(200), status BOOLEAN)'
    
    cursor.execute(sql)
    conn.commit()
    
    sql = 'insert into biodb_config values ("%s", "%s", %s)'
    for row in entry_list:
        cursor.execute(sql % (row[0], row[1], row[2]),)
    conn.commit()
    
def setup_biodb(kwargs):
    sqlpsw = os.environ['SQLPSW']
    db_type = kwargs["chlamdb.db_type"]
    db_name = kwargs["chlamdb.db_name"]
    schema_dir = kwargs["chlamdb.biosql_schema_dir"]
    err_code = 0

    if db_type=="sqlite":
        import sqlite3
        conn = sqlite3.connect(db_name)
        cursor = conn.cursor()
        url_biosql_scheme = 'biosqldb-sqlite.sql'
        err_code = os.system(f"sqlite3 {db_name} < {schema_dir}/{url_biosql_scheme}")
    else:
        import MySQLdb
        conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                    user="root", # your username
                                    passwd=sqlpsw) # name of the data base
        cursor = conn.cursor()
        sql_db = f'CREATE DATABASE IF NOT EXISTS {db_name};'
        cursor.execute(sql_db,)
        conn.commit()
        cursor.execute(f"use {db_name};",)
        url_biosql_scheme = 'biosqldb-mysql.sql'
        err_code = os.system(f"mysql -uroot -p{sqlpsw} {db_name} < {schema_dir}/{url_biosql_scheme}")

    if err_code != 0:
        raise IOError("Problem loading sql schema:", err_code)
    
def setup_chlamdb(**kwargs):
    setup_biodb(kwargs)
    create_data_table(kwargs)

def load_gbk(gbks, args):
    db_name = args["chlamdb.db_name"]

    for gbk in gbks:
        records = [i for i in SeqIO.parse(gbk_name, 'genbank')]
        assert( len(records) == 1)
        db_name.load(records)
        server_name.adaptor.commit()
