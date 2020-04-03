#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


def load_interpro_html(biodb,
                       html_files_tar_gz):
    import tarfile

    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor
    sql1 = 'create table if not exists interpro.hash2interpro_html(hash varchar(300) PRIMARY KEY, interpro_html LONGBLOB)'
    cursor.execute(sql1,)
    for html_file_tar_gz in html_files_tar_gz:
        tar = tarfile.open(html_file_tar_gz, 'r:gz')
        for member in tar.getmembers():
            if "resources" in member.name:
                continue
            print(member.name)
            html_file = tar.extractfile(member)
            content = html_file.read()
            hash = os.path.basename(member.name).split(".")[0]
            print("hash", hash)
            sql = 'insert into interpro.hash2interpro_html values (%s, %s)'
            try:
                cursor.execute(sql, [hash,  content])
            except:
                print("Alredy into db?")
                continue
            conn.commit()
        conn.commit()

if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_interpro_tag_gz', type=str, help="input interpro html tag.gz file(s)", nargs='+')
    parser.add_argument("-d", '--biodb', type=str, help="biodb name")

    args = parser.parse_args()

    load_interpro_html(args.biodb,
                       args.input_interpro_tag_gz)

    manipulate_biosqldb.update_config_table(args.db_name, "interpro_data")