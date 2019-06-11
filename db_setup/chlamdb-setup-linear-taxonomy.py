#!/usr/bin/env python
import MySQLdb
import os
from csv import DictReader
import gzip
import re


class SQLiteDB():
    def __init__(self, dbname='blastnr'):

        sqlpsw = os.environ['SQLPSW']
        conn = MySQLdb.connect(host="localhost",
                               user="root",
                               passwd=sqlpsw)
        cursor = conn.cursor()
        sql_db = 'CREATE DATABASE IF NOT EXISTS %s;' % dbname
        cursor.execute(sql_db,)
        conn.commit()
        cursor.execute("use %s;" % dbname,)

    def importFromCSV(self, csvfilename, tablename, separator=","):
        with gzip.open(csvfilename, 'rt') as fh:
            dr = DictReader(fh, delimiter=separator)
            fieldlist_header = ''
            fieldlist = ''
            for n, field in enumerate(dr.fieldnames):
                field = re.sub(' ', '_', field)
                fieldlist += '`%s`,' % field
                if ("taxid" in field or "tax_id" in field):
                    fieldlist_header += ('%s INTEGER,' % field)
                elif (field == "order"):
                    fieldlist_header += ('`%s` TEXT,' % field)
                else:
                    fieldlist_header += ('%s TEXT,' % field)
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
        sql2 = 'create index ncbit on ncbi_taxonomy(tax_id);'
        self.db.execute(sql2)
        self.db.commit()


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--linear_taxonomy_file', type=str, help='Linear taxonomy file')
    parser.add_argument('-d', '--db_name', type=str, help='DB name', default='blastnr')

    args = parser.parse_args()

    db = SQLiteDB(args.db_name)

    db.importFromCSV(args.linear_taxonomy, args.db_name)
