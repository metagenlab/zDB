#!/usr/bin/python

# As a COG can have several functions, it is necessary to post-process the data
# for those COGs that have several functions :
#  remove the entries and add one to the count of all functions
#  e.g. if COG has functions ECL, with a count of 3, remove the entry
#  and add 3 to entry E, entry C and entry L.
def simplify(hsh_data):
    hsh_simplified = {}
    for bioentry, hsh_funct_to_count in hsh_data.items():
        hsh_entry = {}
        for function, count in hsh_funct_to_count.items():
            for i in range(0, len(function)):
                curr = hsh_entry.get(function[i], 0)
                hsh_entry[function[i]] = curr + count
        hsh_simplified[bioentry] = hsh_entry
    return hsh_simplified

def plot_cog_heatmap(db, ref_tree, bioentry_ids=None, taxon_id_list=None, frequency=False, group_by_cog_id=False):
    from chlamdb.phylo_tree_display import ete_motifs

    if False:
        # NOTE: need to check with Trestan whether this code is still necessary
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

    hsh_data = db.get_COG_counts(bioentry_ids=bioentry_ids, taxon_ids=taxon_id_list)
    hsh_data = simplify(hsh_data)

    if frequency:
        # based on total annotated with COG and not genome size

        bioentry_to_count = {}
        for bioentry, hsh_func_to_count in hsh_data.items():
            bioentry_to_count[bioentry] = sum(hsh_func_to_count.values())
        code2taxon2count = {}
        cog_list = []
    else:
        bioentryid2count_no_GOG = db.get_n_prot_without_cog(taxons=taxon_id_list, bioentries=bioentry_ids)
        bioentryid2proteome_size = db.n_CDS(taxons=taxon_id_list, bioentries=bioentry_ids)

        code2taxon2count = {}
        code2taxon2count['-'] = {}
        code2taxon2count['TOTAL'] = {}
        for bioentry, cnt in bioentryid2count_no_GOG.items():
            if (bioentry_ids != None and bioentry in bioentry_ids):
                code2taxon2count['-'][bioentry] = cnt
                code2taxon2count['TOTAL'][bioentry] = bioentryid2proteome_size[bioentry]

        cog_list = ['TOTAL', '-']

    code2description = db.get_cog_code_description()
    for bioentry, hsh_count in hsh_data.items():
        for funct, count in hsh_count.items():
            descr = f"{code2description[funct]} ({funct})"
            if descr not in cog_list:
                cog_list.append(descr)
            if descr not in code2taxon2count:
                code2taxon2count[descr] = {}
            if frequency:
                code2taxon2count[descr][str(bioentry)] = round((float(count)/(bioentry_to_count[bioentry]))*100,2)
            else:
                code2taxon2count[descr][str(bioentry)] = count

    leaf_to_name = db.get_genomes_description()
    tree2 = ete_motifs.multiple_profiles_heatmap(None,
                                                cog_list,
                                                code2taxon2count,
                                                show_labels=True,
                                                column_scale=True,
                                                tree=ref_tree,
                                                as_float=frequency,
                                                leaf_to_name=leaf_to_name)
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
