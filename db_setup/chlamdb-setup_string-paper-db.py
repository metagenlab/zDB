#!/usr/bin/env python

import psycopg2
import os 
import sqlite3



def setup_sqlitedb(sqlite_db_path):

    pwd = os.environ['SQLPSW']

    pg_connection = psycopg2.connect(user = "tpillone",host = "127.0.0.1",port = "5432",database = "string", password=pwd)
    pg_cursor = pg_connection.cursor()

    sqlite_conn = sqlite3.connect(sqlite_db_path)
    sqlite_cursor = sqlite_conn.cursor()

    sqlite_cursor.execute("PRAGMA synchronous = OFF")
    sqlite_cursor.execute("BEGIN TRANSACTION")


    # setup species table 
    sql_species = 'select species_id,compact_name,kingdom,protein_count from items.species;'
    sql_species_create = 'create table species (species_id INT, compact_name TEXT, kingdom varchar(200), protein_count INT)'
    sqlite_cursor.execute(sql_species_create)
    pg_cursor.execute(sql_species)

    sql_template = 'insert into species values (?,?,?,?)'
    for row in pg_cursor.fetchall():
        sqlite_cursor.execute(sql_template, row)
    sqlite_conn.commit()
    print("ok")
    
    # -- items.proteins --
    # protein_id 
    # protein_external_id 
    # species_id 
    # protein_checksum 
    # protein_size 
    # annotation
    # preferred_name 
    # annotation_word_vectors
    
    sql_proteins = 'select protein_id,protein_external_id,annotation,preferred_name from items.proteins;'
    sql_proteins_create = 'create table proteins (protein_id INT, protein_external_id varchar(400), annotation TEXT, preferred_name varchar(400))'
    sqlite_cursor.execute(sql_proteins_create)
    pg_cursor.execute(sql_proteins)

    n = 0
    entry_list = []
    sql_template = 'insert into proteins values (?,?,?,?)'
    for row in pg_cursor.fetchall():
        n+=1
        entry_list.append(row)
        if n % 100000 == 0:
            print(n)
            sqlite_cursor.executemany(sql_template, entry_list)
            entry_list = []
    sqlite_cursor.executemany(sql_template, entry_list)
    sqlite_conn.commit()
    print("ok")
 
    # -- evidence.publications --
    # publication_id --> PMID:24616884
    # publication_date 
    # publication_source 
    # linkout_url 
    # authors  
    # title 
    # abstract - empty 

    sql_publications = 'select publication_id,publication_date,publication_source,linkout_url,authors,title from evidence.publications;'
    sql_publications_create = 'create table publications (pmid INT, publication_date varchar(400), publication_source TEXT, linkout_url TEXT, authors TEXT, title TEXT)'
    sqlite_cursor.execute(sql_publications_create)
    pg_cursor.execute(sql_publications)

    n = 0
    entry_list = []
    sql_template = 'insert into publications values (?,?,?,?,?,?)'
    for row in pg_cursor.fetchall():
        n+=1
        row = list(row)
        row[0] = row[0].split(":")[1]
        entry_list.append(row)
        if n % 100000 == 0:
            print(n)
            sqlite_cursor.executemany(sql_template, entry_list)
            entry_list = []
    sqlite_cursor.executemany(sql_template, entry_list)

    sqlite_conn.commit()

    # -- evidence.items_publications --
    # item_id 
    # publication_id  PMID:8820243
    # name_shown - putative DNA helicase
    # paragraph_number 
    # sentence_number_in_paragraph 
    # start_position 
    # end_position 
    # go_id
    
    sql_items = 'select item_id,publication_id from evidence.items_publications group by item_id,publication_id;'
    sql_items_create = 'create table protein_id2pmid (protein_id INT, pmid INT);'
    sqlite_cursor.execute(sql_items_create)
    pg_cursor.execute(sql_items)
    
    n = 0
    entry_list = []
    sql_template = 'insert into protein_id2pmid values (?,?)'
    for row in pg_cursor.fetchall():
        n+=1
        row = list(row)
        row[1] = row[1].split(":")[1]
        entry_list.append(row)
        if n % 1000000 == 0:
            print(n)
            sqlite_cursor.executemany(sql_template, entry_list)
            sqlite_conn.commit()
            entry_list = []
    sqlite_cursor.executemany(sql_template, entry_list)

    sqlite_conn.commit()
    

    sql1 = 'create index prid1 on proteins(protein_id);'
    sql2 = 'create index prid2 on proteins(protein_external_id);'
    sql3 = 'create index pmid1 on publications(pmid);'
    sql4 = 'create index prid3 on protein_id2pmid(protein_id);'
    sql5 = 'create index prid4 on protein_id2pmid(pmid);'

    sqlite_cursor.execute(sql1)
    sqlite_cursor.execute(sql2)
    sqlite_cursor.execute(sql3)
    sqlite_cursor.execute(sql4)
    sqlite_cursor.execute(sql5) 

    sqlite_conn.commit()



if __name__ == '__main__':
    import argparse
    import chlamdb_setup_utils
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", '--sqlite_db', type=str, help="sqlite db path", required=True)

    args = parser.parse_args()
    
    setup_sqlitedb(args.sqlite_db)
