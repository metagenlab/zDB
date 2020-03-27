#!/usr/bin/env python

def gbk2taxid(gbk_files, db_name):
    from Bio import SeqIO
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(db_name)
    file_name2taxid = {}

    for gbk in gbk_files:
        records = [i for i in SeqIO.parse(gbk, 'genbank')]
        sql = 'select taxon_id from bioentry t1 ' \
              ' inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id where t1.accession="%s" and t2.name="%s";' % (records[0].name, db_name)
        print(sql)
        file_name2taxid[gbk.split("/")[-1].split(".")[0]] = int(server.adaptor.execute_and_fetchall(sql,)[0][0])
    return file_name2taxid


def load_reference_phylogeny(db_name, newick_file, gbk_files):
    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.plots import parse_newick_tree

    file_name2taxid = gbk2taxid(gbk_files, db_name)
    print(file_name2taxid)

    server, db = manipulate_biosqldb.load_db(db_name)

    sql = 'select biodatabase_id from biodatabase where name="%s"' % (db_name)

    db_id = server.adaptor.execute_and_fetchall(sql,)[0][0]

    sql = 'create table if not exists reference_phylogeny (biodatabase_id INTEGER, tree TEXT)'

    server.adaptor.execute(sql,)

    with open(newick_file, 'r') as f:
        rows = [i for i in f]
        if len(rows) > 1:
            raise IOError("wrong tree format (need newick)")
        newick_string = rows[0].rstrip()
        new_tree = parse_newick_tree.convert_terminal_node_names(newick_string, file_name2taxid, False)

    sql = 'insert into reference_phylogeny values( %s, "%s")' % (db_id, new_tree.write())
    server.adaptor.execute(sql,)
    server.commit()

if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb

    parser = argparse.ArgumentParser()

    parser.add_argument("-d", '--database_name', type=str, help="Database name")
    parser.add_argument("-r", '--reference_phylogeny', type=str, help="Reference phylogeny file (newick, taxids as labels)")
    parser.add_argument("-g", '--gbk_files', type=str, help="gbk files (leaf labels)", nargs='+')

    args = parser.parse_args()
    load_reference_phylogeny(args.database_name, args.reference_phylogeny, args.gbk_files)
    
    
    manipulate_biosqldb.update_config_table(args.database_name, "reference_phylogeny")