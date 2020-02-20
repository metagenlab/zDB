#!/usr/bin/env python



def setup_biodb(biodb_name):
    import urllib.request
    import sys
    import MySQLdb
    import os

    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost",
                        user="root",
                        passwd=sqlpsw)
    cursor = conn.cursor()

    url_mysql_biosql_scheme = 'https://raw.githubusercontent.com/biosql/biosql/master/sql/biosqldb-mysql.sql'

    sys.stdout.write('Downloading Biosql scheme from %s ...\n' % url_mysql_biosql_scheme)
    request = urllib.request.Request(url_mysql_biosql_scheme)
    page = urllib.request.urlopen(request)

    sys.stdout.write("Creating mysql database...\n")

    sql_db = f'CREATE DATABASE IF NOT EXISTS {biodb_name};'
    cursor.execute(sql_db,)
    conn.commit()
    cursor.execute("use {biodb_name};",)

    sys.stdout.write("Importing Biosql schema...\n")


    for line in page.read().decode('unicode-escape').split('\n'):
        if line == '':
            continue
        if line[0] == '-':
            continue
        cursor.execute("%s\n" % line)
    conn.commit()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)


    args = parser.parse_args()

    setup_biodb(args.db_name)