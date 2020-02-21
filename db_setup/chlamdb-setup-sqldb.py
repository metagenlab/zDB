#!/usr/bin/env python



def setup_biodb(biodb_name):
    import urllib.request
    import sys
    import MySQLdb
    import os
    from subprocess import Popen, PIPE

    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost",
                        user="root",
                        passwd=sqlpsw)
    cursor = conn.cursor()

    sys.stdout.write("Creating mysql database...\n")

    sql_db = f'CREATE DATABASE IF NOT EXISTS {biodb_name};'
    cursor.execute(sql_db,)
    conn.commit()
    cursor.execute(f"use {biodb_name};",)

    url_mysql_biosql_scheme = 'https://raw.githubusercontent.com/biosql/biosql/master/sql/biosqldb-mysql.sql'

    sys.stdout.write('Downloading Biosql scheme from %s ...\n' % url_mysql_biosql_scheme)
    request = urllib.request.Request(url_mysql_biosql_scheme)
    page = urllib.request.urlopen(request)
    
    with open("/tmp/biosqldb-mysql.sql", "wb") as f:
        content = page.read()
        f.write(content)

    sys.stdout.write("Importing Biosql schema...\n")
    err_code = os.system(f"mysql -uroot -p{sqlpsw} {biodb_name} < /tmp/biosqldb-mysql.sql")
    if err_code == 0:
        sys.stdout.write("OK")
    else:
        raise IOError("Problem loading sql schema:", err_code)
        
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)


    args = parser.parse_args()

    setup_biodb(args.db_name)