#!/usr/bin/env python


def create_data_table(biodb, sqlite):

    if sqlite:
        import sqlite3

        conn = sqlite3.connect(sqlite)
        cursor = conn.cursor()
    
    else:
        import os
        import MySQLdb
        sqlpsw = os.environ['SQLPSW']

        conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                    user="root", # your username
                                    passwd=sqlpsw, # your password
                                    db=biodb) # name of the data base
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
        ("interpro_taxonomy", "optional", False) # interpro taxnonomy statistics
        ("pfam_taxonomy", "optional", False) #  taxnonomy statistics
        ("COG_taxonomy", "optional", False) # COG taxnonomy statistics
    ]
    
    sql = 'create table biodb_config (name varchar(200), type varchar(200), status BOOLEAN)'
    
    cursor.execute(sql)
    conn.commit()
    
    sql = 'insert into biodb_config values ("%s", "%s", %s)'
    for row in entry_list:
        
        cursor.execute(sql % (row[0], row[1], row[2]),)
    
    conn.commit()  
    
    
    



def setup_biodb(biodb_name,
                sqlite=False):
    import urllib.request
    import sys
    import os
    from subprocess import Popen, PIPE
    
    if sqlite:
        import sqlite3

        conn = sqlite3.connect(sqlite)
        cursor = conn.cursor()
    
    else:
        import MySQLdb
        
        sqlpsw = os.environ['SQLPSW']

        conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                    user="root", # your username
                                    passwd=sqlpsw) # name of the data base
        cursor = conn.cursor()

        sys.stdout.write("Creating mysql database...\n")

        sql_db = f'CREATE DATABASE IF NOT EXISTS {biodb_name};'
        cursor.execute(sql_db,)
        conn.commit()
        cursor.execute(f"use {biodb_name};",)


    if not sqlite:
        url_biosql_scheme = 'https://raw.githubusercontent.com/biosql/biosql/master/sql/biosqldb-mysql.sql'
    else:
        url_biosql_scheme = 'https://raw.githubusercontent.com/biosql/biosql/master/sql/biosqldb-sqlite.sql'

    sys.stdout.write('Downloading Biosql scheme from %s ...\n' % url_biosql_scheme)
    request = urllib.request.Request(url_biosql_scheme)
    page = urllib.request.urlopen(request)
    
    with open("/tmp/biosql.sql", "wb") as f:
        content = page.read()
        f.write(content)

    sys.stdout.write("Importing Biosql schema...\n")
    if not sqlite:
        err_code = os.system(f"mysql -uroot -p{sqlpsw} {biodb_name} < /tmp/biosql.sql")
    else:
        err_code = os.system(f"sqlite3 {sqlite} < /tmp/biosql.sql")
    if err_code == 0:
        sys.stdout.write("OK")
    else:
        raise IOError("Problem loading sql schema:", err_code)
    
    
        
if __name__ == '__main__':
    import argparse
    import os 
        
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)
    parser.add_argument("-s", '--sqlite', action='store_true', help="Use sqlite rather than MySQL")


    args = parser.parse_args()

    if args.sqlite:
        home_folder = os.environ['HOME']
        sqlite_path = f"{home_folder}/{args.db_name}.db"
    else:
        sqlite_path = False

    setup_biodb(args.db_name, 
                sqlite_path)

    create_data_table(args.db_name, 
                      sqlite_path)