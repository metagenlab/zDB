#!/usr/bin/python

def module_list2profile_dico(biodb, module_list, taxon_id_list=[]):

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select biodatabase_id from biodatabase where name="%s"' % biodb

    db_id = server.adaptor.execute_and_fetchall(sql,)[0][0]

    filter_2 = '"'+'","'.join(module_list)+'"'

    # RESTRICT TO AS SUBSET OF THE TAXON AVAILABLE
    if len(taxon_id_list) > 0:
        filter = ',' .join(taxon_id_list)


        sql = 'select taxon_id, description, count(*) from (select distinct taxon_id,t1.ko_id,module_name,description ' \
              ' from enzyme_seqfeature_id2ko t1' \
              ' inner join enzyme_module2ko t2 on t1.ko_id=t2.ko_id ' \
              ' inner join enzyme_kegg_module as t3 on t2.module_id=t3.module_id ' \
              ' inner join annotation_seqfeature_id2locus t4 ' \
              ' on t1.seqfeature_id=t4.seqfeature_id' \
              ' where module_name in (%s) and taxon_id in (%s)) A group by A.taxon_id,A.module_name;' % (filter_2,
                                                                                                         filter)
    else:
        sql = 'select taxon_id, description, count(*) from (select distinct taxon_id,t1.ko_id,module_name,description ' \
              ' from enzyme_seqfeature_id2ko t1' \
              ' inner join enzyme_module2ko t2 on t1.ko_id=t2.ko_id ' \
              ' inner join enzyme_kegg_module as t3 on t2.module_id=t3.module_id ' \
              ' inner join annotation_seqfeature_id2locus t4 ' \
              ' on t1.seqfeature_id=t4.seqfeature_id' \
              ' where module_name in (%s)) A group by A.taxon_id,A.module_name;' % (filter_2)
    print (sql)
    data = server.adaptor.execute_and_fetchall(sql,)



    code2taxon2count = {}
    module_list = []
    for row in data:
        row = list(row)
        '''
        split_name = row[1].split('=>')
        if len(split_name[0]) >0:
            row[1] = split_name[0]
        '''
        #else:

        if row[1] not in module_list:
            module_list.append(row[1])
        if row[1] not in code2taxon2count:
            code2taxon2count[row[1]] = {}
            code2taxon2count[row[1]][str(row[0])] = int(row[2])
        else:
            code2taxon2count[row[1]][str(row[0])] = int(row[2])
    return module_list, code2taxon2count

def plot_module_heatmap(biodb, ref_tree, module_list, taxon_id_list=[], rotate=False):

    from chlamdb.phylo_tree_display import ete_motifs

    module_list, code2taxon2count = module_list2profile_dico(biodb,
                                                             module_list,
                                                             taxon_id_list=taxon_id_list)

    tree2 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                module_list,
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

module_list = ['M00015', 'M00028', 'M00016']

cofactors_and_vitamins =[

"M00112",
#"M00114", ascorbate
"M00115",
"M00116",
"M00117",
"M00118",
"M00119",
"M00120",
"M00121",
#"M00122",Cobalamin biosynthesis
"M00123",
#"M00124",Pyridoxal biosynthesis
"M00125",
"M00126",
"M00127",
"M00128",
"M00129",
#"M00140",C1-unit interconversion, prokaryotes
#"M00141",
"M00550",
"M00572",
"M00573",
"M00577",
"M00622"
]

aa_synthesis = ["M00015",
                "M00016",
                "M00017","M00018","M00019","M00020","M00021","M00022","M00023","M00024","M00025","M00026","M00027",
                "M00028","M00029","M00030","M00031","M00032","M00033","M00034","M00035","M00036","M00037","M00038",
                "M00040","M00042","M00043","M00044","M00045","M00046","M00047","M00048","M00049","M00050","M00051","M00052","M00053"]

atp_synthesis = [
"M00142",
"M00143",
"M00144",
"M00145",
"M00146",
"M00147",
"M00148",
"M00149",
"M00150",
"M00151",
"M00152",
"M00153",
"M00154",
"M00155",
"M00156",
"M00157",
"M00158",
"M00159",
"M00160",
"M00162",
"M00416",
"M00417"
]

'''
tree, style = plot_module_heatmap('chlamydia_04_16',
                            '/home/trestan/work/projets/rhabdo/core_phylo_chlam_staleyi_all_single_copy_07_16/core_chlamydia.tree',
                            atp_synthesis,
                            taxon_list
                            )
tree.render('atp_synthesis.svg', dpi=800, h=600, tree_style=style)
'''
