#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


def load_hmm_data(biodb, database_name, table_name, hmm_tab_files):

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table %s.%s (locus_tag varchar(400), taxon_id INT, hmm_accession varchar(400));' % (database_name,
                                                                                                      table_name)

    server.adaptor.execute_and_fetchall(sql,)

    sql = 'select locus_tag, taxon_id from orthology_detail' % biodb

    locus2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    for hmm_file in hmm_tab_files:
        print hmm_file
        with open(hmm_file, 'r') as f:
            for row in f:
                if row[0] == '#':
                    continue
                else:
                    data = row.rstrip().split()
                    locus_tag = data[0]
                    hit_accession = data[3]
                    print '%s\t%s\t%s' % (hit_accession, locus_tag, locus2taxon_id[locus_tag])
                    sql = 'insert into %s.%s values ("%s", %s, "%s")' % (database_name, table_name,locus_tag, locus2taxon_id[locus_tag], hit_accession)
                    server.adaptor.execute(sql,)
        server.adaptor.commit()

def plot_profile(biodb, database_name, table_name):
    from chlamdb.biosqldb import manipulate_biosqldb
    import phylo_tree_bar

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select tree from reference_phylogeny t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id' \
          ' where t2.name="%s";' % biodb

    phylogeny = server.adaptor.execute_and_fetchall(sql,)[0][0]

    sql = 'select hmm_accession, taxon_id from %s.%s group by hmm_accession, taxon_id' % (database_name, table_name)

    data = server.adaptor.execute_and_fetchall(sql,)

    header_list = []
    taxon2hmm_accession2count = {}
    for row in data:

        if row[0] not in taxon2hmm_accession2count:
            taxon2hmm_accession2count[row[0]] = {}
            taxon2hmm_accession2count[row[0]][row[1]] = 1
        else:
            taxon2hmm_accession2count[row[0]][row[1]] = 1
        if row[0] not in header_list:
            header_list.append(row[0])

    tree, style = phylo_tree_bar.plot_tree_stacked_barplot(phylogeny,
                                                           taxon2value_list_barplot=False,
                                                           header_list=False,  # header stackedbarplots
                                                           taxon2set2value_heatmap=taxon2hmm_accession2count,
                                                           taxon2label=False,
                                                           header_list2=header_list,  # header counts columns
                                                           biodb=biodb,
                                                           column_scale=True,
                                                           general_max=False,
                                                           header_list3=False,
                                                           set2taxon2value_list_simple_barplot=False,
                                                           set2taxon2value_list_simple_barplot_counts=False,
                                                           taxon2description=False,
                                                           rotate=True)

    path = '/home/trestan/hmm_tree.svg'
    style.rotation = 90
    tree.render(path, dpi=800, h=600, tree_style=style)

if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    import parse_priam_EC

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_hmm_tab_files', type=str, help="input hmm tab file", nargs='+', default=False)
    parser.add_argument("-b", '--biodb', type=str, help="biodb_name", default=False)
    parser.add_argument("-d", '--database_name', type=str, help="database name")
    parser.add_argument("-t", '--table_name', help="table name")


    args = parser.parse_args()
    #load_hmm_data(args.biodb, args.database_name, args.table_name, args.input_hmm_tab_files)
    plot_profile(args.biodb, args.database_name, args.table_name)
