#!/usr/bin/python

def plot_cog_eatmap(biodb, ref_tree, taxon_id_list=[], frequency=False):
    import manipulate_biosqldb
    import ete_motifs

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select biodatabase_id from biodatabase where name="%s"' % biodb

    db_id = server.adaptor.execute_and_fetchall(sql,)[0][0]
    
    # RESTRICT TO AS SUBSET OF THE TAXON AVAILABLE
    if len(taxon_id_list) > 0:
        filter = ',' .join(taxon_id_list)

        sql = 'select taxon_id,functon,count(*) as n ' \
              ' from COG.locus_tag2gi_hit_%s t1 ' \
              ' inner join COG.cog_names_2014 t2 on t1.COG_id=t2.COG_id ' \
              ' inner join biosqldb.bioentry as t3 on t1.accession=t3.accession ' \
              ' where biodatabase_id=%s and taxon_id in (%s) group by taxon_id,functon' % (biodb,db_id, filter)
        '''
        sql = 'select taxon_id,functon,count(*) as n from (select taxon_id,functon,t1.COG_id ' \
              ' from COG.locus_tag2gi_hit_%s t1 ' \
              ' inner join COG.cog_names_2014 t2 on t1.COG_id=t2.COG_id ' \
              ' inner join biosqldb.bioentry as t3 on t1.accession=t3.accession ' \
              ' where biodatabase_id=102 and taxon_id in (%s) group by taxon_id, t1.COG_id) A group by taxon_id,functon' % (biodb, filter)
        '''
    else:
        sql = 'select taxon_id,functon,count(*) as n ' \
              ' from COG.locus_tag2gi_hit_%s t1 ' \
              ' inner join COG.cog_names_2014 t2 on t1.COG_id=t2.COG_id ' \
              ' inner join biosqldb.bioentry as t3 on t1.accession=t3.accession ' \
              ' where biodatabase_id=%s group by taxon_id,functon' % (biodb, db_id)
    print sql    
    data = server.adaptor.execute_and_fetchall(sql,)

    if frequency:
        sql = 'select taxon_id,count(*) as n from COG.locus_tag2gi_hit_%s t1 ' \
              ' inner join COG.cog_names_2014 t2 on t1.COG_id=t2.COG_id ' \
              ' inner join biosqldb.bioentry as t3 on t1.accession=t3.accession ' \
              ' where biodatabase_id=%s group by taxon_id;' % (biodb, db_id)
        taxon_id2count = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    code2taxon2count = {}
    cog_list = []
    for row in data:
        if row[1] not in cog_list:
            cog_list.append(row[1])
        if row[1] not in code2taxon2count:
            code2taxon2count[row[1]] = {}
            if frequency:
                code2taxon2count[row[1]][str(row[0])] = round((float(row[2])/float(taxon_id2count[str(row[0])]))*100,2)
            else:
                code2taxon2count[row[1]][str(row[0])] = int(row[2])
        else:
            if frequency:
                code2taxon2count[row[1]][str(row[0])] = round((float(row[2])/float(taxon_id2count[str(row[0])]))*100,2)
            else:
                code2taxon2count[row[1]][str(row[0])] = int(row[2])

    tree2 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                cog_list,
                                                code2taxon2count,
                                                show_labels=True,
                                                column_scale=True,
                                                tree=ref_tree,
                                                as_float=frequency)
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

tree = plot_cog_eatmap('chlamydia_04_16',
                       '/home/trestan/work/projets/rhabdo/core_phylo_chlam_staleyi_all_single_copy_07_16/core_chlamydia.tree',
                       taxon_list,
                       True
                       )
tree.render('test_percent.svg', dpi=800, h=600)
