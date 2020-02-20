#!/usr/bin/env python


def get_module_count_all_db(biodb, category=False):
    '''

    :param biodb: <biodatabase name>
    :param category: KEGG module category (optional)
    :return: for each module, return the total count from KEGG, and the total count of KO present in the <biodb>
    '''

    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

    database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

    if category:

        sql_pathway_count = 'select BB.module_name,count_all,count_db,count_db/count_all from (select module_id, count(*) ' \
                            ' as count_db from (select distinct ko_id from enzyme_locus2ko) as t1' \
                            ' inner join enzyme_module2ko as t2 on t1.ko_id=t2.ko_id group by module_id) AA ' \
                            ' right join (select t1.module_id,module_name, count_all from (select module_id, count(*) as count_all ' \
                            'from enzyme_module2ko group by module_id) t1 inner join enzyme_kegg_module as t2 ' \
                            'on t1.module_id=t2.module_id where module_sub_cat="%s")BB on AA.module_id=BB.module_id;' % (biodb, category) # where pathway_category!="1.0 Global and overview maps"
    else:
        # select distinct KO
        # join with module
        sql_pathway_count = 'select BB.module_name,count_all,count_db,count_db/count_all from (select module_id, count(*) ' \
                            ' as count_db from (select distinct ko_id from enzyme_locus2ko) as t1' \
                            ' inner join enzyme_module2ko as t2 on t1.ko_id=t2.ko_id group by module_id) AA ' \
                            ' right join (select t1.module_id,module_name, count_all from (select module_id, count(*) as count_all ' \
                            'from enzyme_module2ko group by module_id) t1 inner join enzyme_kegg_module as t2 ' \
                            'on t1.module_id=t2.module_id)BB on AA.module_id=BB.module_id;' % (biodb) # where pathway_category!="1.0 Global and overview maps"


    map2count = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_pathway_count,))
    return map2count

def get_module_count_per_genome(biodb, category=False):

    '''

    :return: dictionnary with module category as key
                with nested dictionnary with taxon id as key
                    with nested dictionnary with module id as key and the corresponding list as values (list):
                        1. KO count for this module
                        2. module description


    '''

    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    if category:
        # C.pathway_category,taxon_id, A.pathway_name,A.n_enzymes, C.description
        # select distinct KO id
        # join with module to get counts/modules
        sql = 'select B.module_sub_cat,A.taxon_id,B.module_name,A.n,B.description from ' \
                            ' (select taxon_id, module_id, count(*) as n from ' \
                            ' (select distinct taxon_id,ko_id from enzyme_locus2ko) t1 ' \
                            ' inner join enzyme_module2ko as t2 on t1.ko_id=t2.ko_id group by taxon_id, module_id) A ' \
                            ' inner join enzyme_kegg_module as B on A.module_id=B.module_id where module_sub_cat="%s";' % (biodb, category)
    else:
        sql = 'select B.module_sub_cat,A.taxon_id,B.module_name,A.n,B.description from ' \
                            ' (select taxon_id, module_id, count(*) as n from ' \
                            ' (select distinct taxon_id,ko_id from enzyme_locus2ko) t1 ' \
                            ' inner join enzyme_module2ko as t2 on t1.ko_id=t2.ko_id group by taxon_id, module_id) A ' \
                            ' inner join enzyme_kegg_module as B on A.module_id=B.module_id;' % (biodb)

    pathway_data = server.adaptor.execute_and_fetchall(sql,)
    all_maps = []
    category2maps = {}
    # pathway cat 2 taxon_id 2 pathway_map 2 [count, pathway description]
    pathway_category2taxon2map = {}
    for one_row in pathway_data:
        # first pathway category
        if one_row[0] not in pathway_category2taxon2map:
            category2maps[one_row[0]] = [[one_row[2],one_row[4]]]
            all_maps.append(one_row[2])
            pathway_category2taxon2map[one_row[0]] = {}
            pathway_category2taxon2map[one_row[0]][one_row[1]] = {}
            pathway_category2taxon2map[one_row[0]][one_row[1]][one_row[2]] = one_row[3:]
        else:

            if one_row[2] not in all_maps:
                category2maps[one_row[0]].append([one_row[2],one_row[4]])
                all_maps.append(one_row[2])
            # if new taxon
            if one_row[1] not in pathway_category2taxon2map[one_row[0]]:
                pathway_category2taxon2map[one_row[0]][one_row[1]] = {}

                pathway_category2taxon2map[one_row[0]][one_row[1]][one_row[2]] = one_row[3:]
            # if new map for existing taxon
            else:

                pathway_category2taxon2map[one_row[0]][one_row[1]][one_row[2]] = one_row[3:]
    return pathway_category2taxon2map

def taxon2module2count(biodb, category=False):

    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    if category:
        # C.pathway_category,taxon_id, A.pathway_name,A.n_enzymes, C.description
        # select distinct KO id
        # join with module to get counts/modules
        sql = 'select B.module_sub_cat,A.taxon_id,B.module_name,A.n,B.description from ' \
                            ' (select taxon_id, module_id, count(*) as n from ' \
                            ' (select distinct taxon_id,ko_id from enzyme_locus2ko) t1 ' \
                            ' inner join enzyme_module2ko as t2 on t1.ko_id=t2.ko_id group by taxon_id, module_id) A ' \
                            ' inner join enzyme_kegg_module as B on A.module_id=B.module_id where module_sub_cat="%s";' % (biodb, category)
    else:
        sql = 'select B.module_sub_cat,A.taxon_id,B.module_name,A.n,B.description from ' \
                            ' (select taxon_id, module_id, count(*) as n from ' \
                            ' (select distinct taxon_id,ko_id from enzyme_locus2ko) t1 ' \
                            ' inner join enzyme_module2ko as t2 on t1.ko_id=t2.ko_id group by taxon_id, module_id) A ' \
                            ' inner join enzyme_kegg_module as B on A.module_id=B.module_id;' % (biodb)

    pathway_data = server.adaptor.execute_and_fetchall(sql,)
    taxon2module2c = {}
    # 0 sub cat
    # 1 taxon_id
    # 2 module name
    # 3 count
    # 4 description
    for row in pathway_data:
        if row[1] not in taxon2module2c:
            taxon2module2c[row[1]] = {}
            taxon2module2c[row[1]][row[2]] = row[3]
        else:
            taxon2module2c[row[1]][row[2]] = row[3]
    return taxon2module2c

module2count = get_module_count_all_db('chlam_metagenom_11_16_all')
#category2taxon2module = get_module_count_per_genome('chlam_metagenom_11_16_all')
taxon2module2c = taxon2module2count('chlam_metagenom_11_16_all')

with open('module_freq.tab', 'w') as g:
    g.write('\t"' + '"\t"'.join([str(i) for i in taxon2module2c.keys()])+'"\n')
    for module in module2count:
        one_module_list = []
        for taxon in taxon2module2c:
            try:
                one_module_list.append(str(taxon2module2c[taxon][module]/float(module2count[module][1])))
            except KeyError:
                one_module_list.append('0')
        if sum([float(i) for i in one_module_list]) == 0:
            continue
        else:
            g.write(module + '\t' + '\t'.join(one_module_list)+'\n')

with open('module_counts.tab', 'w') as f:
    f.write('\t"' + '"\t"'.join([str(i) for i in taxon2module2c.keys()])+'"\n')
    for module in module2count:
        one_module_list = []
        for taxon in taxon2module2c:
            try:
                one_module_list.append(str(taxon2module2c[taxon][module]))
            except KeyError:
                one_module_list.append('0')
        if sum([float(i) for i in one_module_list]) == 0:
            continue
        else:
            f.write(module + '\t' + '\t'.join(one_module_list)+'\n')




