#!/usr/bin/python

def plot_cog_eatmap(biodb, ref_tree, taxon_id_list=[], frequency=False, group_by_cog_id=False):
    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.phylo_tree_display import ete_motifs

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select biodatabase_id from biodatabase where name="%s"' % biodb

    db_id = server.adaptor.execute_and_fetchall(sql,)[0][0]
    
    # RESTRICT TO AS SUBSET OF THE TAXON AVAILABLE

    sql = ''

    if len(taxon_id_list) > 0:
        filter = ',' .join(taxon_id_list)

        sql = 'select taxon_id, code, count(*) as n from COG_seqfeature_id2best_COG_hit t1 ' \
              ' inner join biosqldb.bioentry t2 on t1.bioentry_id=t2.bioentry_id' \
              ' inner join COG_cog_id2cog_category t3 on t1.hit_cog_id=t3.COG_id ' \
              ' inner join COG_code2category t4 on t3.category_id=t4.category_id ' \
              ' where t2.biodatabase_id=%s and taxon_id in (%s)' \
              ' group by taxon_id, code;' % (biodb,
                db_id,
                filter)

        print (sql)
    else:
        if not group_by_cog_id:
            sql = 'select taxon_id,functon,count(*) as n ' \
                  ' from COG_locus_tag2gi_hit t1 ' \
                  ' inner join COG_cog_names_2014 t2 on t1.COG_id=t2.COG_id ' \
                  ' inner join biosqldb.bioentry as t3 on t1.accession=t3.accession ' \
                  ' where biodatabase_id=%s group by taxon_id,functon' % (biodb, db_id)
        else:
            sql = ' select A.taxon_id,B.functon,count(*) from (select t1.COG_id, t3.taxon_id from COG_locus_tag2gi_hit t1 ' \
                  ' inner join orthology_detail t3 on t1.locus_tag=t3.locus_tag ' \
                  ' group by taxon_id,t1.COG_id) A inner join COG_cog_names_2014 B on A.COG_id=B.COG_id ' \
                  ' group by A.taxon_id,B.functon;' % (biodb, biodb)

    data = server.adaptor.execute_and_fetchall(sql,)

    if frequency:
        '''
        ATTENTION: based on total annotated with COG and not genome size
        
        '''
        sql = 'select taxon_id, count(*) as n from COG_seqfeature_id2best_COG_hit t1' \
              ' inner join biosqldb.bioentry t2 on t1.bioentry_id=t2.bioentry_id' \
              ' where t2.biodatabase_id=%s group by taxon_id;' % (biodb, db_id)
        taxon_id2count = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

        code2taxon2count = {}
        cog_list = []

    else:
        sql = 'select taxon_id, count(*) from orthology_detail t1 left join COG_locus_tag2gi_hit t2 ' \
              ' on t1.locus_tag=t2.locus_tag where COG_id is NULL group by t1.taxon_id;' % (biodb,  biodb)

        taxon2count_no_GOG = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

        sql = 'select taxon_id, count(*) from orthology_detail group by taxon_id' % biodb

        taxon2proteome_size = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

        code2taxon2count = {}
        code2taxon2count['-'] = {}
        code2taxon2count['TOTAL'] = {}
        for taxon in taxon2count_no_GOG:
            if taxon in taxon_id_list:
                code2taxon2count['-'][taxon] = int(taxon2count_no_GOG[taxon])
                code2taxon2count['TOTAL'][taxon] = int(taxon2proteome_size[taxon])

        cog_list = ['TOTAL', '-']

    sql = 'select code, description from COG_code2category;'
    code2description = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    for row in data:
        descr = "%s (%s)" % (code2description[row[1]], row[1])
        if descr not in cog_list:
            cog_list.append(descr)
        if descr not in code2taxon2count:
            code2taxon2count[descr] = {}
            if frequency:
                code2taxon2count[descr][str(row[0])] = round((float(row[2])/float(taxon_id2count[str(row[0])]))*100,2)
            else:
                code2taxon2count[descr][str(row[0])] = int(row[2])
        else:
            if frequency:
                code2taxon2count[descr][str(row[0])] = round((float(row[2])/float(taxon_id2count[str(row[0])]))*100,2)
            else:
                code2taxon2count[descr][str(row[0])] = int(row[2])



    tree2 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                cog_list,
                                                code2taxon2count,
                                                show_labels=True,
                                                column_scale=True,
                                                tree=ref_tree,
                                                as_float=frequency)
    return tree2


'''
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

from chlamdb.biosqldb import manipulate_biosqldb
from ete3 import Tree
server, db = manipulate_biosqldb.load_db('chlamydia_04_16')

sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % 'chlamydia_04_16'

tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]

t1 = Tree(tree)

tree, style = plot_cog_eatmap('chlamydia_04_16',
                       tree,#'/home/trestan/work/projets/rhabdo/core_phylo_chlam_staleyi_all_single_copy_07_16/core_chlamydia.tree',
                       [],#taxon_list,
                       True
                       )



tree.render('test_percent.svg', dpi=800, h=600, tree_style=style)
'''