#!/usr/bin/python


def get_single_set_data(biodb, set_name,column="bitscore", bitscore_cutoff=0, query_coverage_cutoff=0):

    '''

    :param biodb:
    :param set_name:
    :param score: should one of those:  bitscore, bias, evalue, query_coverage
    :return:
    '''
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select  from hmm_sets t1 inner join hmm_sets_entry t2 on t1.set_id=t2.set_id ' \
          ' inner join hmm_hits_annotated_genome_chlamydia_04_16 t3 ' \
          ' on t2.hmm_id=t3.hmm_id inner join custom_tables.locus2seqfeature_id_chlamydia_04_16 t4 ' \
          ' on t3.seqfeature_id=t4.seqfeature_id where t1.name="T3SS"'

    sql = 'select t3.taxon_id,t4.name,%s from hmm.hmm_sets t1 ' \
          ' inner join hmm.hmm_sets_entry t2 on t1.set_id=t2.set_id ' \
          ' inner join hmm_hmm_hits_annotated_genome t3 on t2.hmm_id=t3.hmm_id ' \
          ' inner join hmm.hmm_profiles t4 on t2.hmm_id=t4.hmm_id ' \
          ' inner join custom_tables_locus2seqfeature_id t5 on t3.seqfeature_id=t5.seqfeature_id' \
          ' where t1.name="%s" and bitscore>=%s and query_coverage>=%s order by bitscore;' % (column,
                                                                            biodb,
                                                                            biodb,
                                                                            set_name,
                                                                            bitscore_cutoff,
                                                                            query_coverage_cutoff)
    gene2taxon2score = {}
    gene_list = []
    data = server.adaptor.execute_and_fetchall(sql,)
    for row in data:
        if row[1] not in gene_list:
            gene_list.append(row[1])
        if row[1] not in gene2taxon2score:
            gene2taxon2score[row[1]] = {}
            gene2taxon2score[row[1]][str(row[0])] = row[2]
        else:
            if str(row[0]) not in gene2taxon2score[row[1]]:
                gene2taxon2score[row[1]][str(row[0])] = row[2]
            else:
                print ('already present!!!!!!!')
                continue
    return gene2taxon2score, gene_list


def get_set_data(biodb,
                 set_list_restrict=[],
                 frequency=False,
                 six_frame_translation=False,
                 return_lists=False,
                 score_cutoff=0):
    from chlamdb.biosqldb import manipulate_biosqldb
    '''

    :param biodb:
    :param set_list: restrict analysis to specific sets (empty list mean all sets)
    :param frequency: return ratio n genes identified/n genes iun the set
    :param cutoff_percent: onle return presence absence data (1 and 0) given a cutoff percentage of the genes identified/genes in the set
    :param six_frame_translation: get data from the six fram translation analysis
    :return: dictionnary taxon2list of values OR taxon2set2value dictionnary
    '''

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select biodatabase_id from biodatabase where name="%s"' % biodb

    db_id = server.adaptor.execute_and_fetchall(sql,)[0][0]

    if six_frame_translation:
        hmm_table = 'hmm_hits_six_frame_genome'
    else:
        hmm_table = 'hmm_hits_annotated_genome'

    sql = 'select taxon_id,set_id, count(*) from ' \
          ' (select t1.*,t2.set_id from hmm.%s_%s t1 ' \
          ' inner join hmm.hmm_sets_entry t2 on t1.hmm_id=t2.hmm_id where t1.bitscore>%s' \
          ' group by taxon_id,set_id,t1.hmm_id) A group by taxon_id,set_id;' % (hmm_table, biodb, score_cutoff)

    data = server.adaptor.execute_and_fetchall(sql,)

    if frequency:
        sql = 'select taxon_id,count(*) as n from COG_locus_tag2gi_hit t1 ' \
              ' inner join COG_cog_names_2014 t2 on t1.COG_id=t2.COG_id ' \
              ' inner join biosqldb.bioentry as t3 on t1.accession=t3.accession ' \
              ' where biodatabase_id=%s group by taxon_id;' % (biodb, db_id)
        taxon_id2count = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select * from hmm.hmm_sets'
    set_id2description = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    set2taxon2count = {}
    taxon2list = {}
    set_list = []
    for row in data:
        if row[0] not in taxon2list:
            taxon2list[row[0]] = [row[1]]
        else:
            taxon2list[row[0]].append(row[1])

        set = set_id2description[str(row[1])]
        if set not in set_list:
            set_list.append(set)
        if len(set_list_restrict) > 0:
            if set not in set_list_restrict:
                continue
        if set not in set2taxon2count:
            set2taxon2count[set] = {}

            if frequency:
                freq = round((float(row[2])/float(taxon_id2count[str(row[0])]))*100,2)
                set2taxon2count[set][str(row[0])] = freq
            else:
                set2taxon2count[set][str(row[0])] = int(row[2])
        else:
            if frequency:
                freq = round((float(row[2])/float(taxon_id2count[str(row[0])]))*100,2)
                set2taxon2count[set][str(row[0])] = freq
            else:
                set2taxon2count[set][str(row[0])] = int(row[2])

    if not return_lists:
        return set2taxon2count, set_list
    else:
        return taxon2list, set_list




def plot_hmm_heatmap(biodb, ref_tree, taxon_id_list=[], frequency=False, six_frame_translation=False):
    from chlamdb.biosqldb import manipulate_biosqldb
    import ete_motifs

    code2taxon2count, set_list = get_set_data(biodb)

    tree2 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                set_list,
                                                code2taxon2count,
                                                show_labels=True,
                                                column_scale=True,
                                                tree=ref_tree,
                                                as_float=frequency)
    return tree2


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d",'--biodb',type=str,help="database_name")

    args = parser.parse_args()

    taxon_list = ["67",
    "1279767",
    "1279774",
    "1279496",
    "48",
    "46",
    "55",
    "87925",
    "1279815",
    "62",
    "1279822",
    "66",
    "52",
    "49",
    "64",
    "60",
    "804807",
    "886707",
    "283",
    "314",
    "1069693",
    "1069694",
    "1137444",
    "1143376",
    "313",
    "1172027",
    "1172028",
    "1035343",
    "307",
    "293",
    "1279839",
    "1279497"]


    biodb = args.biodb #'chlamydia_04_16'
    from chlamdb.biosqldb import manipulate_biosqldb
    from ete2 import Tree
    import ete_motifs
    import cog_heatmap
    server, db = manipulate_biosqldb.load_db(biodb)
    sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

    tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]

    t1 = Tree(tree)

    taxon_list = []

    tree, style = plot_hmm_heatmap(biodb,
                           tree,
                           taxon_list,
                           False
                           )

    tree.render('test_count.svg', dpi=800, h=600,tree_style=style)
