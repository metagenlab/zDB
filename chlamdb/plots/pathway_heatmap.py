#!/usr/bin/python





def pathway_list2profile_dico(biodb, pathway_list, taxon_id_list=[], group_by_KO=True):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select biodatabase_id from biodatabase where name="%s"' % biodb

    db_id = server.adaptor.execute_and_fetchall(sql,)[0][0]

    filter_2 = '"'+'","'.join(pathway_list)+'"'

    # RESTRICT TO AS SUBSET OF THE TAXON AVAILABLE
    if len(taxon_id_list) > 0:
        filter = ',' .join(taxon_id_list)


        sql = 'select taxon_id, description, count(*) from (select distinct taxon_id,description,t3.pathway_id,t1.ko_id ' \
              ' from enzyme_locus2ko t1 inner join enzyme_pathway2ko t2 on t1.ko_id=t2.ko_id ' \
              ' inner join enzyme_kegg_pathway as t3 on t2.pathway_id=t3.pathway_id ' \
              ' where t3.pathway_name in (%s) and taxon_id in (%s)) A group by taxon_id,pathway_id;' % (biodb,
                                                                                           filter_2,
                                                                                           filter)
        print (sql)
    else:
        if not group_by_KO:
            print ('not grouping')
            sql = 'select taxon_id,description,count(*) from (select distinct taxon_id,description,t3.pathway_id,t1.ko_id ' \
                  ' from enzyme_locus2ko t1 inner join enzyme_pathway2ko t2 on t1.ko_id=t2.ko_id ' \
                  ' inner join enzyme_kegg_pathway as t3 on t2.pathway_id=t3.pathway_id ' \
                  ' where t3.pathway_name in (%s)) A group by taxon_id,pathway_id;' % (biodb,
                                                                                               filter_2)
        else:
            print ('grouping')
            sql = 'select taxon_id, description, count(*) from (select distinct taxon_id,description,t4.ko_accession ' \
                  ' from enzyme_seqfeature_id2ko t1 ' \
                  ' inner join annotation_seqfeature_id2locus tb on t1.seqfeature_id=tb.seqfeature_id ' \
                  ' inner join enzyme_pathway2ko t2 on t1.ko_id=t2.ko_id  inner join enzyme_kegg_pathway as t3 on t2.pathway_id=t3.pathway_id' \
                  ' inner join enzyme_ko_annotation t4 on t1.ko_id=t4.ko_id where pathway_name in (%s) ) A ' \
                  ' group by A.taxon_id, A.description;;' % (biodb,
                                                             biodb,
                                                             filter_2)
            print (sql)
    #print sql
    data = server.adaptor.execute_and_fetchall(sql,)



    code2taxon2count = {}
    pathway_list = []
    for row in data:
        row = list(row)
        #print row
        if row[1] not in pathway_list:
            pathway_list.append(row[1])
        if row[1] not in code2taxon2count:
            code2taxon2count[row[1]] = {}
            code2taxon2count[row[1]][str(row[0])] = int(row[2])
        else:
            code2taxon2count[row[1]][str(row[0])] = int(row[2])
    return pathway_list, code2taxon2count

def plot_module_and_pathway_combinaison_heatmap(biodb,
                                                ref_tree,
                                                pathway_list,
                                                module_list,
                                                taxon_id_list=[],
                                                group_by_KO=True,
                                                rotate=False):
    from chlamdb.biosqldb import manipulate_biosqldb
    import ete_motifs
    import module_heatmap

    pathway_list, code2taxon2count_pathway = pathway_list2profile_dico(biodb, pathway_list, taxon_id_list=taxon_id_list,
                                                                       group_by_KO=group_by_KO)
    module_list, code2taxon2count_modules = module_heatmap.module_list2profile_dico(biodb, module_list, taxon_id_list=taxon_id_list)

    merged_list = sorted(pathway_list + module_list)
    code2taxon2count_pathway.update(code2taxon2count_modules)


    tree2 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                merged_list,
                                                code2taxon2count_pathway,
                                                show_labels=True,
                                                column_scale=True,
                                                tree=ref_tree,
                                                as_float=False,
                                                rotate=rotate)
    return tree2




def plot_pathway_heatmap(biodb, ref_tree, pathway_list, taxon_id_list=[], rotate=False, group_by_KO=True):
    from chlamdb.biosqldb import manipulate_biosqldb
    import ete_motifs

    pathway_list, code2taxon2count = pathway_list2profile_dico(biodb, pathway_list, taxon_id_list=taxon_id_list, group_by_KO=group_by_KO)

    tree2 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                pathway_list,
                                                code2taxon2count,
                                                show_labels=True,
                                                column_scale=True,
                                                tree=ref_tree,
                                                as_float=False,
                                                rotate=rotate)
    return tree2



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

cofactors_and_vitamins =[

"map00130",
"map00670",
"map00730",
"map00740",
"map00750",
"map00760",
"map00770",
"map00780",
"map00785",
"map00790",
"map00830",
"map00860"
]

aa_sythesis = [
"map00220",
"map00250",
"map00260",
"map00270",
"map00280",
"map00290",
"map00300",
"map00310",
"map00330",
"map00340",
"map00350",
"map00360",
"map00380",
"map00400"
]
'''
from chlamdb.biosqldb import manipulate_biosqldb
from ete3 import Tree
server, db = manipulate_biosqldb.load_db('chlamydia_04_16')

sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % 'chlamydia_04_16'

tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]
#print tree
t1 = Tree(tree)
tree, style = plot_pathway_heatmap('chlamydia_04_16',
                            tree, #'/home/trestan/work/projets/rhabdo/core_phylo_chlam_staleyi_all_single_copy_07_16/core_chlamydia.tree',
                            cofactors_and_vitamins,
                            [], #taxon_list
                            )
tree.render('test_pathways2.svg', dpi=800, h=600, tree_style=style)
'''