#!/usr/bin/env python
import MySQLdb
import os
from csv import DictReader
import gzip
import re


class MySQLDB():
    def __init__(self, biodb):

        from chlamdb.biosqldb import manipulate_biosqldb
        
        server, db = manipulate_biosqldb.load_db(biodb)
        self.mysql_conn = server.adaptor.conn
        self.mysql_cursor = server.adaptor.cursor

        self.biodb = biodb

    def importFromCSV(self, 
                      csvfilename, 
                      tablename, 
                      separator=","):
        with gzip.open(csvfilename, 'rt') as fh:
            dr = DictReader(fh, delimiter=separator)
            fieldlist_header = ''
            fieldlist = ''
            for n, field in enumerate(dr.fieldnames):
                field = re.sub(' ', '_', field)
                fieldlist += '`%s`,' % field
                if ("taxid" in field or "tax_id" in field):
                    fieldlist_header += ('taxon_id INTEGER,')
                elif (field == "order"):
                    fieldlist_header += ('`%s` varchar(400),' % field)
                else:
                    fieldlist_header += ('%s varchar(400),' % field)
            fieldlist = fieldlist[0:-1]

            ph = ("%s,"*len(dr.fieldnames))[:-1]
            self.db.execute("DROP TABLE IF EXISTS %s" % tablename)
            print("CREATE TABLE %s(%s)" % (tablename, fieldlist_header[0:-1]))
            self.db.execute("CREATE TABLE %s(%s)" % (tablename, fieldlist_header[0:-1]))
            ins = "insert into %s (%s) values (%s)" % (tablename, fieldlist, ph)
            for line in dr:
                v = []
                for k in dr.fieldnames: v.append(line[k])
                self.db.execute(ins, v)

        self.db.commit()
        sql2 = 'create index ncbit on blastrn_ncbi_taxonomy(tax_id);'
        self.db.execute(sql2)
        self.db.commit()


    def import_from_sqlite3(self, 
                            sqlite_db_path, 
                            table_name,
                            sqlite_format=False):

        import sqlite3

        sqlite_conn = sqlite3.connect(sqlite_db_path)
        sqlite_cursor = sqlite_conn.cursor()

        sql_header = 'PRAGMA table_info(ncbi_taxonomy);'
        sqlite_cursor.execute(sql_header)
        columns = sqlite_cursor.fetchall()
        column_index = []
        columns_def = []
        for n, col in enumerate(columns):
            col = list(col)
            if col[1] in ['cohort', 'cohort_taxid', 'parvorder_taxid', 'parvorder', 'infraorder', 'infraorder_taxid', 'forma_taxid', 'forma', 'infraclass', 'infraclass_taxid', 'series', 'series_taxid', 'species_group', 'species_group_taxid', 'species1', 'species1_taxid', 'species_subgroup', 'species_subgroup_taxid', 'subclass', 'subclass_taxid', 'subcohort', 'subcohort_taxid', 'subfamily', 'subfamily_taxid', 'subgenus', 'subgenus_taxid', 'subkingdom', 'subkingdom_taxid', 'suborder', 'suborder_taxid', 'subphylum', 'subphylum_taxid', 'subsection', 'subsection_taxid', 'subspecies', 'subspecies_taxid', 'subtribe', 'subtribe_taxid', 'superclass', 'superclass_taxid', 'tribe', 'tribe_taxid', 'varietas', 'varietas_taxid', 'section', 'section_taxid', 'superphamily', 'superfamily_taxid', 'superorder', 'superorder_taxid']:
                continue
            else:
                column_index.append("`%s`" % col[1])
                if n == 0:
                    if col[1] != 'taxon_id':
                        col[1] = 'taxon_id'
                columns_def.append("`%s` %s" % (col[1], col[2]))

        sql = 'select %s from ncbi_taxonomy' % ','.join(column_index)
        print(sql)
        sqlite_cursor.execute(sql)
        taxonomy_data = sqlite_cursor.fetchall()

        column_index[0] = 'taxon_id'

        sql_header_crate = 'create table if not exists blastnr_blastnr_taxonomy (%s)' % ','.join(columns_def)
        print(sql_header_crate)
        self.mysql_cursor.execute(sql_header_crate)

        sql = 'insert into blastnr_blastnr_taxonomy (%s) values (' % ','.join(column_index)
        if not sqlite_format:
            sql += ','.join(["%s"]*len(column_index))
        else:
            # use sqlite syntax
            sql += ','.join(["?"]*len(column_index))
        sql += ')'
        for n, row in enumerate(taxonomy_data):
            if n % 10000 == 0:
                print(n)
            row = list(row)
            row = [i if (i != '') else None for i in row]
            self.mysql_cursor.execute(sql, row)
        self.mysql_conn.commit()
        sql1 = 'CREATE INDEX taxid ON blastnr_blastnr_taxonomy(`taxon_id`);'
        # _mysql_exceptions.OperationalError: (1170, "BLOB/TEXT column 'phylum' used in key specification without a key length")
        # specifiy a max key length
        if not sqlite_format:
            # need to constrain index size for mysql
            sql2 = 'CREATE INDEX ph ON blastnr_blastnr_taxonomy(`phylum`(255));'
        else:
            sql2 = 'CREATE INDEX ph ON blastnr_blastnr_taxonomy(phylum);'
        sql3 = 'CREATE INDEX phid ON blastnr_blastnr_taxonomy(`phylum_taxid`);'

        self.mysql_cursor.execute(sql1)
        self.mysql_cursor.execute(sql2)
        self.mysql_cursor.execute(sql3)
        self.mysql_conn.commit()


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--linear_taxonomy_file', type=str, help='Linear taxonomy file')
    parser.add_argument('-s', '--linear_taxonomy_sqlite', type=str, help='Linear taxonomy sqlite3 file')
    parser.add_argument('-d', '--db_name', type=str, help='DB name', default='Biodb name')
    parser.add_argument('-sf', '--sqlitef', action='store_true', help='Data stored in sqlite rather than MySQL (need to adapt inserts)')

    args = parser.parse_args()




    if args.linear_taxonomy_file:
        db = MySQLDB(args.db_name)
        db.importFromCSV(args.linear_taxonomy, args.db_name)
        
        manipulate_biosqldb.update_config_table(args.db_name, "taxonomy_table")  
           
            
    if args.linear_taxonomy_sqlite:
        db = MySQLDB(args.db_name)
        db.import_from_sqlite3(args.linear_taxonomy_sqlite, 
                               args.db_name, 
                               args.sqlitef)
    
        manipulate_biosqldb.update_config_table(args.db_name, "taxonomy_table")  
    
