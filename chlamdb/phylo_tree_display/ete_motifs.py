#!/usr/bin/env python

import sys
from random import sample
from random import randint
from ete3 import Tree, SeqMotifFace, TreeStyle, add_face_to_node, TextFace

import numpy as np
import colorsys
import matplotlib.colors as pcolors
from django.conf import settings

def _get_colors(num_colors):
    colors=[]
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i/360.
        lightness = (50 + np.random.rand() * 10)/100.
        saturation = (90 + np.random.rand() * 10)/100.
        rgb = colorsys.hls_to_rgb(hue, lightness, saturation)

        colors.append(pcolors.rgb2hex(rgb))
    return colors


def get_interpro2taxon_id2count(biodb, 
                             id_list, 
                             analysis="Pfam", 
                             taxon_filter=False):


    '''
    get presence/absence of pfam domain(s) in all organisms of database "biodb"
    return it as a dictionnary taxon[domain_id] --> count
    :param biodb:
    :param id_list: lit of COgs/interpro/pfam identifiers
    :param type: COG/Pfam/interpro or any comparative table stored in "comparative_tables" mysql database
    :return:
    '''

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    id_filter = '"' + '","'.join(id_list)+ '"'
    sql = f'select taxon_id,signature_accession, count(*) as n from' \
          f' (select distinct t1.seqfeature_id,taxon_id,signature_accession,signature_description ' \
          f' from interpro_interpro t1 ' \
          f' inner join interpro_signature t2 on t1.signature_id=t2.signature_id ' \
          f' inner join interpro_analysis t3 on t2.analysis_id=t3.analysis_id ' \
          f' inner join annotation_seqfeature_id2locus t4 on t1.seqfeature_id=t4.seqfeature_id' \
          f' where analysis_name="{analysis}" and signature_accession in ({id_filter})) A ' \
          f' group by taxon_id,signature_accession,signature_description;'

    data = server.adaptor.execute_and_fetchall(sql,)

    interpro2taxon_id2count = {}
    for i, row in enumerate(data):
        taxon_id = str(row[0])
        signature_accession = row[1]
        n = row[2] 
        if signature_accession not in interpro2taxon_id2count:
            interpro2taxon_id2count[signature_accession] = {}
        
        interpro2taxon_id2count[signature_accession][taxon_id] = n
    #print("interpro2taxon_id2count",interpro2taxon_id2count)
    return interpro2taxon_id2count



def get_taxon2name2count(biodb, 
                         id_list, 
                         type="COG", 
                         taxon_filter=False):

    '''
    get presence/absence of pfam domain(s) in all organisms of database "biodb"
    return it as a dictionnary taxon[pfam_id] --> count
    :param biodb:
    :param id_list: lit of COgs/interpro/pfam identifiers
    :param type: COG/Pfam/interpro or any comparative table stored in "comparative_tables" mysql database
    :return:
    '''

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db =manipulate_biosqldb.load_db(biodb)

    ortho_sql = '"' + '","'.join(id_list) + '"'

    if type !='orthogroup':
        if taxon_filter:
            col_filter = 'id,`' + '`,`'.join(taxon_filter) + '`'
            ordered_taxons = taxon_filter
        else:
            col_filter = '*'
            if db_driver == 'mysql':
                sql = 'show columns from comparative_tables_%s' % (type)
            elif db_driver == 'sqlite':
                sql = 'PRAGMA table_info(comparative_tables_%s)' % (type)
            
            ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)][1:]
        sql = 'select %s from comparative_tables_%s where id in (%s)' % (col_filter, type, ortho_sql)
    else:
        if taxon_filter:
            col_filter = 'orthogroup,`' + '`,`'.join(taxon_filter) + '`'
            ordered_taxons = taxon_filter
        else:
            col_filter = '*'
            sql = 'show columns from comparative_tables_orthology'
            if db_driver == 'mysql':
                sql = 'show columns from comparative_tables_orthology'
            elif db_driver == 'sqlite':
                sql = 'PRAGMA table_info(comparative_tables_orthology)'
            ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)][1:]
        sql = 'select %s from comparative_tables_orthology where orthogroup in (%s)' % (col_filter, ortho_sql)
    #print sql
    profile_tuples = list(server.adaptor.execute_and_fetchall(sql,))

    taxon2group2n_homologs = {}
    for i, tuple in enumerate(profile_tuples):
        taxon2group2n_homologs[tuple[0]] = {}
        for i, taxon in enumerate(ordered_taxons):
            taxon2group2n_homologs[tuple[0]][str(taxon)] = tuple[i+1]

    return taxon2group2n_homologs


def get_taxon2orthogroup2count(biodb, orthogroup_id_list):

    '''
    get presence/absence of pfam domain(s) in all organisms of database "biodb"
    return it as a dictionnary taxon[pfam_id] --> count
    :param biodb:
    :param id_list: lit of COgs/interpro/pfam identifiers
    :param type: COG/Pfam/interpro or any comparative table stored in "comparative_tables" mysql database
    :return:
    '''

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db =manipulate_biosqldb.load_db(biodb)

    if db_driver == 'mysql':
        sql = 'show columns from comparative_tables_orthology'
    elif db_driver == 'sqlite':
        sql = 'PRAGMA table_info(comparative_tables_orthology);'
    else:
        raise IOError("Unknown db driver %s" % db_driver)

    ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)][1:]

    ortho_sql = '"' + '","'.join(orthogroup_id_list) + '"'

    sql = 'select * from comparative_tables_orthology where orthogroup in (%s)' % (ortho_sql)

    profile_tuples = list(server.adaptor.execute_and_fetchall(sql,))

    taxon2group2n_homologs = {}
    for i, tuple in enumerate(profile_tuples):
        taxon2group2n_homologs[tuple[0]] = {}
        for i, taxon in enumerate(ordered_taxons):
            taxon2group2n_homologs[tuple[0]][str(taxon)] = tuple[i+1]

    return taxon2group2n_homologs


def get_locus2taxon2identity(biodb, locus_tag_list):


    '''
    get presence/absence of pfam domain(s) in all organisms of database "biodb"
    return it as a dictionnary taxon[pfam_id] --> count
    :param biodb:
    :param id_list: lit of COgs/interpro/pfam identifiers
    :param type: COG/Pfam/interpro or any comparative table stored in "comparative_tables" mysql database
    :return:
    '''

    from chlamdb.biosqldb import manipulate_biosqldb
    from django.conf import settings
    
    db_driver = settings.DB_DRIVER

    server, db =manipulate_biosqldb.load_db(biodb)

    if db_driver == 'mysql':
        sql = 'show columns from comparative_tables_orthology'
    elif db_driver == 'sqlite':
        sql = 'PRAGMA table_info(comparative_tables_orthology);'
    else:
        raise IOError("Unknown db driver %s" % db_driver)

    ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)][1:]

    taxon2locus2identity_closest = {}
    for locus in locus_tag_list:
        taxon2locus2identity_closest[locus] = {}

    '''
    select B.locus_tag,C.locus_tag,identity from (select * from comparative_tables.identity_closest_homolog2_chlamydia_04_16 where locus_1 in (17323786,17275213,17295471,17323784,17129361,17219096,17194936,17100717,17219094,17129363,17100715,17222853,17067094,17335066,17047317,17215018,17176196,17294883,17047319,17013495)) A left join custom_tables.locus2seqfeature_id_chlamydia_04_16 B on A.locus_1=B.seqfeature_id left join custom_tables.locus2seqfeature_id_chlamydia_04_16 C on A.locus_2=C.seqfeature_id;
    '''
    #print 'getting dico'
    sql = 'select seqfeature_id,locus_tag from custom_tables_locus2seqfeature_id'

    seqfeature_id2locus = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    #print 'ok'
    all_id = []
    for locus in locus_tag_list:
        #print 'locus', locus
        seqfeatre_id_sql = 'select seqfeature_id from custom_tables_locus2seqfeature_id where locus_tag="%s";' % (locus)
        #print seqfeatre_id_sql
        try:
            seqfeatre_id = server.adaptor.execute_and_fetchall(seqfeatre_id_sql,)[0][0]
            all_id.append(str(seqfeatre_id))
        except:
            # pseudogene
            continue
    filter = ','.join(all_id)
    sql = 'select taxon_2,locus_1,identity from comparative_tables_identity_closest_homolog2 where locus_1 in (%s) ;' % (filter)

    identity_data = server.adaptor.execute_and_fetchall(sql, )
    taxon2identity_closest = {}

    for row in identity_data:
        #taxon2identity_closest[str(row[0])] = row[2]
        #taxon2locus2identity_closest[seqfeature_id2locus[str(row[1])]][str(row[0])] = row[2]
        #print taxon2locus2identity_closest
        taxon2locus2identity_closest[seqfeature_id2locus[str(row[1])]][str(row[0])] = row[2]
    for locus in locus_tag_list:
        for taxon in ordered_taxons:
            if taxon not in taxon2locus2identity_closest[locus]:
                taxon2locus2identity_closest[locus][taxon] = '-'



        return taxon2locus2identity_closest

def get_locus2taxon2n_paralogs(biodb, locus_tag_list):

    '''
    get presence/absence of pfam domain(s) in all organisms of database "biodb"
    return it as a dictionnary taxon[pfam_id] --> count
    :param biodb:
    :param id_list: lit of COgs/interpro/pfam identifiers
    :param type: COG/Pfam/interpro or any comparative table stored in "comparative_tables" mysql database
    :return:
    '''

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db =manipulate_biosqldb.load_db(biodb)

    #print 'get_locus2taxon2n_paralogs!!!!!!!!!!!!!!!!'

    if db_driver == 'mysql':
        sql = 'show columns from comparative_tables_orthology'
    elif db_driver == 'sqlite':
        sql = 'PRAGMA table_info(comparative_tables_orthology);'

    ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)][1:]

    taxon2locus2n_paralogs = {}
    for locus in locus_tag_list:
        taxon2locus2n_paralogs[locus] = {}


    for locus in locus_tag_list:

        orthogroup_sql = 'select orthogroup from chlamdb.orthology_detail where locus_tag="%s";' % (biodb,locus)

        orthogroup = server.adaptor.execute_and_fetchall(orthogroup_sql,)[0][0]

        sql = 'select taxon_id,count(*) from chlamdb.orthology_detail where orthogroup="%s" group by taxon_id;' % (orthogroup)
        taxon_id2count = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

        for taxon in ordered_taxons:
            if taxon not in taxon_id2count:
                taxon2locus2n_paralogs[locus][str(taxon)] = 0
            else:
                taxon2locus2n_paralogs[locus][str(taxon)] = taxon_id2count[str(taxon)]
    return taxon2locus2n_paralogs


def combined_profiles_heatmap(biodb,
                              column_labels,
                              taxon2group2count,
                              taxon2motif2count,
                              ec2orthogroups,
                              rotate=False,
                              column_scale=False):

    '''
    motives (or EC or cogs) are sometimes missed by the annotation. One possibility is to
    propagate the prediction of one orthogroup member to the whole orthogroup. This is what this vizualization do:
    create a profile with hilights:
    blue: any of the group with target profile is present
    red: target profile itself is present
    :param biodb:
    :param column_labels:
    :param taxon2group2count:
    :param taxon2motif2count:
    :return:
    '''

    #print ec2orthogroups
    from ete3 import TreeStyle
    from chlamdb.biosqldb import manipulate_biosqldb
    from matplotlib.colors import rgb2hex

    server, db = manipulate_biosqldb.load_db(biodb)

    sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

    tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]

    t1 = Tree(tree)

    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    t1.ladderize()

    taxon_id2organism_name = manipulate_biosqldb.taxon_id2genome_description(server, biodb)

    tss=TreeStyle()
    tss.show_branch_support = False
    tss.draw_guiding_lines = True
    tss.guiding_lines_color = "gray"
    tss.show_leaf_name = False
    if rotate:
        tss.rotation = 90

    if column_scale:
        import matplotlib.cm as cm
        from matplotlib.colors import rgb2hex
        import matplotlib as mpl
        column2scale = {}
        column2max = {}
        for column in column_labels:
            values = [float(i) for i in taxon2motif2count[column].values()]
            #print values, column
            norm = mpl.colors.Normalize(vmin=min(values), vmax=max(values))

            cmap = cm.OrRd
            m = cm.ScalarMappable(norm=norm, cmap=cmap)
            column2scale[column] = m
            if min(values) == max(values):
                column2max[column] = max(values)*2
            else:
                column2max[column] = max(values)

    head = True
    for lf in t1.iter_leaves():
        lf.branch_vertical_margin = 0
        first_column = True
        for col, value in enumerate(column_labels):
            if head:
                #'first row, print gene names'
                #print 'ok!'
                n = TextFace(' %s ' % str(value))
                #n.vt_align = 2
                #n.hz_align = 2
                n.rotation = 270
                #n.margin_top = 2
                #n.margin_right = 2
                #n.margin_left = 2
                n.margin_bottom = 6

                n.inner_background.color = "white"
                n.opacity = 1.
                if rotate:
                    n.rotation = 270
                tss.aligned_header.add_face(n, col)
            try:
                if rotate:
                    if int(taxon2motif2count[value][lf.name]) < 10:
                        n = TextFace('  %s  ' % str(taxon2motif2count[value][lf.name]))
                    elif int(taxon2motif2count[value][lf.name]) < 100:
                        n = TextFace(' %s ' % str(taxon2motif2count[value][lf.name]))
                    else:
                        n = TextFace('%s ' % str(taxon2motif2count[value][lf.name]))
                else:
                    n = TextFace(' %s ' % str(taxon2motif2count[value][lf.name]))
                if column_scale:
                    if float(taxon2motif2count[value][lf.name]) >= 0.6*column2max[value]:
                        n.fgcolor = 'white'
            except:
                 n = TextFace(' - ')
            n.margin_top = 2
            n.margin_right = 2
            n.margin_left = 2
            n.margin_bottom = 2
            # if motif + ==> red
            try:
                count = taxon2motif2count[value][lf.name]
            except:
                count = 0
            if count >0:
                # red
                #n.inner_background.color = "#FA5858"
                # blue
                if not column_scale:
                    n.inner_background.color = "#58ACFA"
                else:
                    n.inner_background.color = rgb2hex(column2scale[value].to_rgba(float(count)))

            else:
                orthologue = False
                try:
                    grp_lst = ec2orthogroups[value]
                    for orthogroup in grp_lst:
                        if taxon2group2count[orthogroup][lf.name] > 0:
                            orthologue = True
                except:
                    orthologue = False

                # if orthogroup + ==> blue
                if orthologue:
                    n.inner_background.color = "#9FF781"
                # if no orthologue ==> white
                else:
                    n.inner_background.color = 'white'
            if rotate:
               n.rotation= 270
            lf.add_face(n, col, position="aligned")
        n = TextFace(taxon_id2organism_name[lf.name], fgcolor = "black", fsize = 12, fstyle = 'italic')
        lf.add_face(n, 0)

        head=False
    for n in t1.traverse():
       from ete3 import NodeStyle
       nstyle = NodeStyle()
       if n.support > 1:
           limit = 100
       else:
           limit = 1

       if n.support < limit:
           nstyle["fgcolor"] = "black"
           nstyle["size"] = 4
           n.set_style(nstyle)
       else:
           nstyle["fgcolor"] = "red"
           nstyle["size"] = 0
           n.set_style(nstyle)


    return t1, tss


def reverse_colourmap(cmap, name = 'my_cmap_r'):
    import matplotlib as mpl
    """
    In:
    cmap, name
    Out:
    my_cmap_r

    Explanation:
    t[0] goes from 0 to 1
    row i:   x  y0  y1 -> t[0] t[1] t[2]
                   /
                  /
    row i+1: x  y0  y1 -> t[n] t[1] t[2]

    so the inverse should do the same:
    row i+1: x  y1  y0 -> 1-t[0] t[2] t[1]
                   /
                  /
    row i:   x  y1  y0 -> 1-t[n] t[2] t[1]
    """
    reverse = []
    k = []

    for key in cmap._segmentdata:
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:
            data.append((1-t[0],t[2],t[1]))
        reverse.append(sorted(data))

    LinearL = dict(zip(k,reverse))
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL)
    return my_cmap_r

def pathways_heatmap(biodb,
                      category2maps,
                      category2taxon2map2count):

    '''
    motives (or EC or cogs) are sometimes missed by the annotation. One possibility is to
    propagate the prediction of one orthogroup member to the whole orthogroup. This is what this vizualization do:
    create a profile with hilights:
    blue: any of the group with target profile is present
    red: target profile itself is present
    :param biodb:
    :param column_labels:
    :param taxon2group2count:
    :param taxon2motif2count:
    :return:
    '''

    from chlamdb.biosqldb import manipulate_biosqldb
    import matplotlib.cm as cm
    from matplotlib.colors import rgb2hex
    import matplotlib as mpl





    server, db = manipulate_biosqldb.load_db(biodb)

    sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

    tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]

    t1 = Tree(tree)

    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    t1.ladderize()




    taxon_id2organism_name = manipulate_biosqldb.taxon_id2genome_description(server, biodb)

    taxon2map2color = {}
    for category in category2maps:
        for map in category2maps[category]:
            count_list = []
            taxon2count = {}
            for taxon in taxon_id2organism_name.keys():
                taxon=int(taxon)
                try:
                    count_list.append(float(category2taxon2map2count[category][taxon][map[0]][0]))
                    taxon2count[taxon] = category2taxon2map2count[category][taxon][map[0]][0]
                except:
                    continue
            if len(count_list)>0:
                #print 'count list',map, count_list, category2taxon2map2count[category][taxon]
                #if min(count_list) == max(count_list):
                count_list.append(0)
                norm = mpl.colors.Normalize(vmin=min(count_list), vmax=max(count_list)) # map2count[map[0]][0]
                cmap_blue = cm.Blues
                m2 = cm.ScalarMappable(norm=norm, cmap=cmap_blue)
                taxon2map2color[map[0]] = {}
                norm = mpl.colors.Normalize(vmin=min(count_list), vmax=max(count_list)) # map2count[map[0]][0]
                norm2 = mpl.colors.Normalize(vmin=min(count_list), vmax=max(count_list))
                cmap_bw = cm.afmhot#YlOrBr
                #cmap_bw_r = reverse_colourmap(cmap_bw)
                m1 = cm.ScalarMappable(norm=norm2, cmap=cmap_bw)
                m2 = cm.ScalarMappable(norm=norm, cmap=cmap_blue)
                taxon2map2color[map[0]] = {}



                for taxon in taxon2count:
                    taxon2map2color[map[0]][taxon] = [rgb2hex(m2.to_rgba(float(taxon2count[taxon]))),
                                                      rgb2hex(m1.to_rgba(float(taxon2count[taxon])))]


    head = True
    for lf in t1.iter_leaves():
        lf.branch_vertical_margin = 0
        first_column = True
        col = 0

        colors = []

        for x, category in enumerate(category2taxon2map2count):

            for y, map in enumerate(category2maps[category]):
                col +=1
                #print "column", col
                taxon = lf.name
                #print category2taxon2map2count[category][int(taxon)]
                try:
                    map_data = category2taxon2map2count[category][int(taxon)][map[0]]
                except:
                    #print 'map', map[0]
                    if head:

                        #print 'ok!'
                        n = TextFace('%s (%s)' % (map[1], category))
                        n.vt_align = 2
                        n.hz_align = 2
                        n.rotation= 270
                        n.margin_top = 0
                        n.margin_right = 0
                        n.margin_left = 0
                        n.margin_bottom = 0
                        if x%2==0:
                            n.inner_background.color = "#F79F81"
                        else:
                            n.inner_background.color = "white"
                        n.opacity = 1.
                        lf.add_face(n, col, position="aligned")

                    n = TextFace('  ')
                    n.margin_top = 1
                    n.margin_right = 1
                    n.margin_left = 1
                    n.margin_bottom = 1
                    # if motif + ==> red
                    n.inner_background.color = "#E0F8E0"

                else:
                    #print "map_data", map_data
                    #print 'map', map
                    if head:

                        #print 'ok!'
                        n = TextFace('%s (%s)' % (map[1], category))
                        n.vt_align = 2
                        n.hz_align = 2
                        n.rotation= 270
                        n.margin_top = 0
                        n.margin_right = 0
                        n.margin_left = 0
                        n.margin_bottom = 0

                        if x%2==0:
                            n.inner_background.color = "#F79F81"
                        else:
                            n.inner_background.color = "white"
                        n.opacity = 1.
                        lf.add_face(n, col, position="aligned")

                    n = TextFace(' %s ' % str(map_data[0]))
                    n.margin_top = 0
                    n.margin_right = 0
                    n.margin_left = 0
                    n.margin_bottom = 0
                    # if motif + ==> red
                    if map_data[0] > 0:
                        #n.inner_background.color = "#FA5858"
                        try:
                            n.inner_background.color = taxon2map2color[map[0]][int(taxon)][0] #"#58ACFA"
                            n.fgcolor = taxon2map2color[map[0]][int(taxon)][1]
                        except:
                            n.inner_background.color = "#58ACFA"

                    else:
                        n.inner_background.color = "#9FF781"
                    if taxon_id2organism_name[lf.name]=="Rhabdochlamydia helveticae T3358":
                        n.inner_background.color = "red"
                        n.fgcolor='white'


                lf.add_face(n, col, position="aligned")
        lf.name = taxon_id2organism_name[lf.name]
        head=False

    return t1

def pathways_heatmapV2(biodb,
                      category2maps,
                      category2taxon2map2count):

    '''
    motives (or EC or cogs) are sometimes missed by the annotation. One possibility is to
    propagate the prediction of one orthogroup member to the whole orthogroup. This is what this vizualization do:
    create a profile with hilights:
    blue: any of the group with target profile is present
    red: target profile itself is present
    :param biodb:
    :param column_labels:
    :param taxon2group2count:
    :param taxon2motif2count:
    :return:
    '''

    from chlamdb.biosqldb import manipulate_biosqldb
    import matplotlib.cm as cm
    from matplotlib.colors import rgb2hex
    import matplotlib as mpl





    server, db = manipulate_biosqldb.load_db(biodb)

    sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

    tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]

    acc_list = [u'NZ_CP006571', u'NZ_AWUS01000000', u'NZ_APJW00000000', u'NC_002620', u'NZ_CCJF00000000', u'NZ_AYKJ01000000', u'NC_000117', u'LNES01000000', u'LJUH01000000', u'NC_004552', u'NC_003361', u'NC_007899', u'NC_015408', u'NC_000922', u'NC_015470', u'NZ_CCEJ000000000', u'CWGJ01000001', u'NZ_JSDQ00000000', u'NZ_BASK00000000', u'NZ_JRXI00000000', u'NZ_BAWW00000000', u'NZ_ACZE00000000', u'NC_015702', u'NZ_BBPT00000000', u'NZ_JSAN00000000', u'NC_005861', u'FCNU01000001', u'NZ_LN879502', u'NZ_BASL00000000', u'Rht', u'CCSC01000000', u'NC_015713', u'NC_014225']

    filter = '"' + '","'.join(acc_list) + '"'
    sql = 'select taxon_id from bioentry where biodatabase_id=102 and accession in (%s)' % filter

    taxon_list = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql,)]



    t1 = Tree(tree)

    #for leaf in t1.iter_leaves():
    #    print type(leaf.name)

    t2 = t1.prune(taxon_list)


    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    t1.ladderize()




    taxon_id2organism_name = manipulate_biosqldb.taxon_id2genome_description(server, biodb)

    taxon2map2color = {}
    for category in category2maps:
        for map in category2maps[category]:
            count_list = []
            taxon2count = {}
            for taxon in taxon_id2organism_name.keys():
                if str(taxon) not in taxon_list:
                    continue
                taxon=int(taxon)
                try:
                    count_list.append(float(category2taxon2map2count[category][taxon][map[0]][0]))
                    taxon2count[taxon] = category2taxon2map2count[category][taxon][map[0]][0]
                except:
                    continue
            if len(count_list)>0:
                #print 'count list',map, count_list, category2taxon2map2count[category][taxon]
                #if min(count_list) == max(count_list):
                count_list.append(0)
                norm = mpl.colors.Normalize(vmin=min(count_list), vmax=max(count_list)) # map2count[map[0]][0]
                cmap_blue = cm.Blues
                m2 = cm.ScalarMappable(norm=norm, cmap=cmap_blue)
                taxon2map2color[map[0]] = {}
                norm = mpl.colors.Normalize(vmin=min(count_list), vmax=max(count_list)) # map2count[map[0]][0]
                norm2 = mpl.colors.Normalize(vmin=min(count_list), vmax=max(count_list))
                cmap_bw = cm.afmhot#YlOrBr
                #cmap_bw_r = reverse_colourmap(cmap_bw)
                m1 = cm.ScalarMappable(norm=norm2, cmap=cmap_bw)
                m2 = cm.ScalarMappable(norm=norm, cmap=cmap_blue)
                taxon2map2color[map[0]] = {}



                for taxon in taxon2count:
                    taxon2map2color[map[0]][taxon] = [rgb2hex(m2.to_rgba(float(taxon2count[taxon]))),
                                                      rgb2hex(m1.to_rgba(float(taxon2count[taxon])))]

    longest = 0
    longest_name = ''
    for x, category in enumerate(category2taxon2map2count):

        for y, map in enumerate(category2maps[category]):
            name = '%s (%s)' % (map[0], map[1].split('=>')[-1].split('[')[0])
            if len(name) > longest:
                longest = len(name)
                longest_name = name

    head = True
    for lf in t1.iter_leaves():
        lf.branch_vertical_margin = 0
        first_column = True
        col = 0

        colors = []

        for x, category in enumerate(category2taxon2map2count):

            for y, map in enumerate(category2maps[category]):
                col +=1
                #print "column", col
                taxon = lf.name
                #print category2taxon2map2count[category][int(taxon)]
                try:
                    map_data = category2taxon2map2count[category][taxon][map[0]]
                except:
                    #print 'map', map[0]
                    if head:

                        #print 'ok!'
                        name = '%s (%s)' % (map[0], map[1].split('=>')[-1].split('[')[0])
                        diff = longest - len(name) +2
                        name = name + ('_'*diff)
                        #print name, len(name)
                        n = TextFace(name)
                        #n.vt_align = 1
                        #n.hz_align = 1
                        n.rotation= 270
                        n.margin_top = 0
                        n.margin_right = 0
                        n.margin_left = 0
                        n.margin_bottom = 0
                        if x%2==0:
                            n.inner_background.color = "white"
                        else:
                            n.inner_background.color = "white"
                        n.opacity = 1.
                        lf.add_face(n, col, position="aligned")

                    n = TextFace(' 0 ')
                    n.margin_top = 1
                    n.margin_right = 1
                    n.margin_left = 1
                    n.margin_bottom = 1
                    # if motif + ==> red
                    n.fgcolor = "#E0F8F7"
                    n.inner_background.color = "#E0F8F7" #"#E0F8E0"

                else:
                    # first row has no enzyme for the considered module
                    #print "map_data", map_data
                    #print 'map', map
                    if head:

                        name = '%s (%s)' % (map[0], map[1].split('=>')[-1].split('[')[0])
                        diff = longest - len(name)
                        name = name + ('_'*diff)
                        #print name, len(name)
                        n = TextFace(name)
                        #n.vt_align = 1
                        #n.hz_align = 1
                        n.rotation= 270
                        n.margin_top = 0
                        n.margin_right = 0
                        n.margin_left = 0
                        n.margin_bottom = 0

                        if x%2==0:
                            n.inner_background.color = "white"
                        else:
                            n.inner_background.color = "white"
                        n.opacity = 1.
                        lf.add_face(n, col, position="aligned")

                    n = TextFace(' %s ' % str(map_data[0]))
                    n.rotation = 270
                    n.margin_top = 0
                    n.margin_right = 0
                    n.margin_left = 0
                    n.margin_bottom = 0
                    # if motif + ==> red
                    if map_data[0] > 0:
                        #n.inner_background.color = "#FA5858"
                        try:
                            n.inner_background.color = taxon2map2color[map[0]][taxon][0] #"#58ACFA"
                            n.fgcolor = taxon2map2color[map[0]][taxon][1]
                        except:
                            n.inner_background.color = "#58ACFA"

                    else:
                        n.inner_background.color = "#9FF781"
                    if taxon_id2organism_name[lf.name]=="helveticae T3358":
                        n.inner_background.color = "red"
                        n.fgcolor='white'


                lf.add_face(n, col, position="aligned")
        lf.name = taxon_id2organism_name[lf.name]
        head=False

    return t1




def multiple_profiles_heatmap(biodb,
                              column_labels,
                              group2taxon2count,
                              reference_taxon=False,
                              reference_column = False,
                              taxon2group2value=False,
                              highlight_first_column=False,
                              identity_scale=False,
                              column_scale=False,
                              show_labels=True,
                              tree=False,
                              as_float=False,
                              rotate=False,
                              sqlite3=False,
                              leaf_to_name = None):

    '''

    ATTENTION: dico inverse: il faut group2taxon2values

    :param biodb:
    :param column_labels:
    :param taxon2group2count:
    :param reference_taxon:
    :param reference_column:
    :param taxon2group2value: color or not background based on
    :param highlight_first_column:
    :return:
    '''
    if isinstance(reference_taxon, dict):
        import copy
        reference_taxon_dico = copy.copy(reference_taxon)

        reference_taxon = False
    else:
        reference_taxon_dico = False


    from chlamdb.biosqldb import manipulate_biosqldb

    if identity_scale:
        import matplotlib.cm as cm
        from matplotlib.colors import rgb2hex
        import matplotlib as mpl

        norm = mpl.colors.Normalize(vmin=0, vmax=100)
        cmap = cm.OrRd
        m = cm.ScalarMappable(norm=norm, cmap=cmap)

    elif column_scale:
        import matplotlib.cm as cm
        from matplotlib.colors import rgb2hex
        import matplotlib as mpl
        column2scale = {}
        column2max = {}
        for column in column_labels:
            values = [float(i) for i in group2taxon2count[column].values()]
            norm = mpl.colors.Normalize(vmin=min(values), vmax=max(values))
            cmap = cm.OrRd
            m = cm.ScalarMappable(norm=norm, cmap=cmap)
            column2scale[column] = m
            if min(values) == max(values):
                column2max[column] = max(values)*2
            else:
                column2max[column] = max(values)

    if not tree:
        server, db = manipulate_biosqldb.load_db(biodb)
        sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

        tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]
    if not isinstance(tree, Tree):
        t1 = Tree(tree)

    else:
        t1 = tree
    tss = TreeStyle()
    tss.show_branch_support = False
    tss.draw_guiding_lines = True
    tss.guiding_lines_color = "gray"
    tss.show_leaf_name = False
    tss.draw_aligned_faces_as_table=True

    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    t1.ladderize()

    if leaf_to_name == None:
        server, db = manipulate_biosqldb.load_db(biodb)
        leaf_to_name = manipulate_biosqldb.taxon_id2genome_description(server, biodb,filter_names=True)

    head = True
    leaf_list = [i for i in t1.iter_leaves()]
    n_leaves = len(leaf_list)

    longest_label = max(len(i) for i in column_labels)
    for lf_count, lf in enumerate(t1.iter_leaves()):
        lf.branch_vertical_margin = 0

        first_column = True
        for col, value in enumerate(column_labels):
            if lf_count == 0:
                    # add labels
                    diff = (longest_label - len(value))
                    n = TextFace('%s%s' % (str(value), diff * ' '))
                    n.vt_align = 0
                    n.hz_align = 0
                    n.rotation= 270
                    if col == 0:
                        n.margin_left = 3
                    else:
                        n.margin_bottom = 2
                        n.margin_left = 2
                    #n.margin_top = 2
                    #n.margin_right = 2
                    #n.margin_left = 2
                    n.margin_bottom = 6

                    from ete3 import NodeStyle
                    n.inner_background.color = "white"
                    n.opacity = 1.
                    tss.aligned_header.add_face(n, col)
                    #lf.add_face(n, col, position="aligned")

                    #nstyle = NodeStyle()
                    #nstyle["fgcolor"] = "red"
                    #nstyle["size"] = 15
                    #lf.set_style(nstyle)

            if first_column and not reference_column and highlight_first_column:
                
                # highlight of the first column only (red)
                try:
                    #print(value, group2taxon2count[value][lf.name])
                    n = TextFace(' %s ' % str(group2taxon2count[value][lf.name]))

                except:
                    n = TextFace(' - ')
                    #n = RectFace(width=10, height=10, fgcolor="red", bgcolor="blue", label='-')
                    if value not in group2taxon2count:
                        group2taxon2count[value] = {}
                    group2taxon2count[value][lf.name] = 0
                ref_data = str(value)
                n.margin_top = 1
                n.margin_right = 1
                n.margin_left = 1
                n.margin_bottom = 1

                if group2taxon2count[value][lf.name] > 0 and group2taxon2count[value][lf.name] != '-':
                    n.inner_background.color = "#FA5858"
                else:
                    n.inner_background.color = 'white'
                first_column = False
            else:
                if not taxon2group2value:
                    if value not in group2taxon2count:
                        n = TextFace(' - ')
                    else:
                        if lf.name not in group2taxon2count[value]:
                            n = TextFace(' - ')
                        else:
                            if group2taxon2count[value][lf.name] == '-':
                                n = TextFace(' - ')
                            else:
                                # if identity scale: 2 digit format
                                local_label = str(group2taxon2count[value][lf.name])
                                if not identity_scale:# and not column_scale:
                                    #print 'not identity scale!'
                                    if as_float:
                                        local_label = "%.2f" % group2taxon2count[value][lf.name]
                                    else:
                                        try:
                                            #print('ok!')
                                            local_label = "%s" % int(group2taxon2count[value][lf.name])
                                        except:
                                            local_label = "%s" % group2taxon2count[value][lf.name]
                                else:
                                    if round(float(group2taxon2count[value][lf.name]), 2) < 100:
                                        if type(group2taxon2count[value][lf.name]) != int:
                                            local_label = "%.2f" % round(group2taxon2count[value][lf.name], 2)
                                        else:
                                            local_label = "%s" % group2taxon2count[value][lf.name]
                                    else:
                                        if type(group2taxon2count[value][lf.name]) != int:
                                            local_label = "%.1f" % round(group2taxon2count[value][lf.name], 2)
                                        else:
                                            local_label = "%s" % group2taxon2count[value][lf.name]
                                if show_labels:
                                    try:
                                        if round(group2taxon2count[value][lf.name], 2) < 100 and column_scale:
                                            if round(group2taxon2count[value][lf.name], 2) < 10:
                                                #print('less than 10: %s' % group2taxon2count[value][lf.name])
                                                n = TextFace(' %s ' % local_label)
                                            else:
                                                #print('Beteen 10 and 100: %s' % group2taxon2count[value][lf.name])
                                                n = TextFace('%s' % local_label)
                                        else:
                                            #print('more than 100: %s' % group2taxon2count[value][lf.name])
                                            n = TextFace('%s' % local_label)
                                        if column_scale:
                                            if float(group2taxon2count[value][lf.name]) >= 0.6 * column2max[value]:
                                                n.fgcolor = 'white'
                                    # labels are not floats
                                    except TypeError:
                                        #print("not float")
                                        n = TextFace(' %s ' % group2taxon2count[value][lf.name])
                                else:
                                    n = TextFace(' - ')

                    '''
                    except:

                        if not reference_taxon_dico:
                            n = TextFace(' - ')

                        # if doctionnary locus2taxon report which taxon is the reference
                        # color the cell in blue
                        else:
                            if str(reference_taxon_dico[value]) == lf.name and show_labels:
                                n = TextFace(' 100.0 ')
                                n.fgcolor = '#58ACFA'
                            else:
                                n = TextFace(' - ')
                    '''
                    if show_labels:
                        n.margin_top = 1
                        n.margin_right = 1
                        n.margin_left = 1
                        n.margin_bottom = 1
                    else:
                        n.margin_top = 0
                        n.margin_right = 0
                        n.margin_left = 0
                        n.margin_bottom = 5
                    try:
                        count = int(group2taxon2count[value][lf.name])
                    except:
                        count = 0
                    #print value, lf.name
                    #print taxon2group2count[value][lf.name]
                    if count > 0 and group2taxon2count[value][lf.name] != '-':
                        if not reference_column:
                            if lf.name != str(reference_taxon):
                                # clor given a continuous scale
                                if identity_scale:
                                    n.inner_background.color = rgb2hex(m.to_rgba(float(count)))
                                    if float(count) >90:
                                        n.fgcolor='white'
                                    if not show_labels:
                                        n.fgcolor = rgb2hex(m.to_rgba(float(count)))
                                elif column_scale:
                                    n.inner_background.color = rgb2hex(column2scale[value].to_rgba(float(count)))
                                else:
                                    n.inner_background.color = "#58ACFA"
                            else:
                                n.inner_background.color = '#58ACFA'
                                n.fgcolor = "white"

                        else:
                            if lf.name == str(reference_taxon) or col == reference_column:
                                n.inner_background.color = "#FA5858"
                                n.fgcolor = "white"
                            else:
                                n.inner_background.color = "#58ACFA"
                    else:

                        if not reference_taxon_dico:
                            n.inner_background.color = 'white'
                        # if doctionnary locus2taxon report which taxon is the reference
                        # color the cell in blue
                        else:
                            if str(reference_taxon_dico[value]) == lf.name:
                                n.inner_background.color = '#58ACFA'
                            else:
                                n.inner_background.color = 'white'
                else:

                    try:
                        n = TextFace(' %s ' % str(group2taxon2count[value][lf.name]))
                    except:
                        group2taxon2count[value][lf.name] = 0
                        n = TextFace(' - ')
                    n.margin_top = 1
                    n.margin_right = 1
                    n.margin_left = 1
                    n.margin_bottom = 1

                    if group2taxon2count[value][lf.name] > 0 and group2taxon2count[value][lf.name] != '-':
                        if not reference_column:
                            #print 'no ref column'
                            if not reference_taxon:
                                if lf.name != str(reference_taxon):
                                    #print 'group', value
                                    #print 'ref data', ref_data
                                    #print taxon2group2value[str(lf.name)][value]
                                    #print(taxon2group2value)
                                    #print("test", taxon2group2value[str(lf.name)][value])
                                    try:
                                        if ref_data not in taxon2group2value[str(lf.name)][value]:
                                            n.inner_background.color = "#9FF781" #58ACFA
                                        else:
                                            n.inner_background.color = "#F78181" ##F6D8CE
                                    except:
                                        n.inner_background.color = "#9FF781"
                                else:
                                    n.inner_background.color = "#FA5858"
                            else:
                                #print 'group', lf.name
                                if ref_data not in taxon2group2value[lf.name][value]:
                                    n.inner_background.color = "#FA5858"
                                else:
                                    n.inner_background.color = "#FA5858"
                        else:
                            #print 'ref column'
                            if lf.name == str(reference_taxon) or col == reference_column:
                                n.inner_background.color = "#FA5858"
                                n.fgcolor = "white"

                            else:
                                n.inner_background.color = "#58ACFA"

                    else:
                        n.inner_background.color = 'white'
            if rotate:
               n.rotation= 270
            #print ('margin left: %s' % n.margin_left)
            lf.add_face(n, col, position="aligned")

        #lf.name = taxon_id2organism_name[lf.name]
        try:
            n = TextFace(leaf_to_name[lf.name], fgcolor = "black", fsize = 12, fstyle = 'italic')
        except:
            n = TextFace(lf.name, fgcolor = "black", fsize = 12, fstyle = 'italic')
        lf.add_face(n, 0)
        head=False
    for n in t1.traverse():
       from ete3 import NodeStyle
       nstyle = NodeStyle()
       if n.support > 1:
           limit = 100
       else:
           limit = 1

       if n.support < limit:
           nstyle["fgcolor"] = "black"
           nstyle["size"] = 4
           n.set_style(nstyle)
       else:
           nstyle["fgcolor"] = "red"
           nstyle["size"] = 0
           n.set_style(nstyle)

    return t1, tss

def multiple_profiles_heatmap_nobiodb(column_labels,
                              group2taxon2count,
                              taxon_id2organism_name,
                              reference_taxon=False,
                              reference_column = False,
                              taxon2group2value=False,
                              highlight_first_column=False,
                              identity_scale=False,
                              column_scale=False,
                              show_labels=True,
                              tree=False,
                              as_float=False,
                              rotate=False):

    '''

    ATTENTION: dico inverse: il faut group2taxon2values

    :param biodb:
    :param column_labels:
    :param taxon2group2count:
    :param reference_taxon:
    :param reference_column:
    :param taxon2group2value: color or not background based on
    :param highlight_first_column:
    :return:
    '''

    if isinstance(reference_taxon, dict):
        import copy
        reference_taxon_dico = copy.copy(reference_taxon)

        reference_taxon = False
    else:
        reference_taxon_dico = False


    from chlamdb.biosqldb import manipulate_biosqldb

    if identity_scale:
        import matplotlib.cm as cm
        from matplotlib.colors import rgb2hex
        import matplotlib as mpl

        norm = mpl.colors.Normalize(vmin=0, vmax=100)
        cmap = cm.OrRd
        m = cm.ScalarMappable(norm=norm, cmap=cmap)

    elif column_scale:
        import matplotlib.cm as cm
        from matplotlib.colors import rgb2hex
        import matplotlib as mpl
        column2scale = {}
        for column in column_labels:
            values = [float(i) for i in group2taxon2count[column].values()]
            #print values, column
            norm = mpl.colors.Normalize(vmin=min(values), vmax=max(values))
            cmap = cm.OrRd
            m = cm.ScalarMappable(norm=norm, cmap=cmap)
            column2scale[column] = m

    if not isinstance(tree, Tree):
        t1 = Tree(tree)

    else:
        t1 = tree
    tss = TreeStyle()
    tss.show_branch_support = False
    tss.draw_guiding_lines = True
    tss.guiding_lines_color = "gray"
    tss.show_leaf_name = False

    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    t1.ladderize()



    head = True
    leaf_list = [i for i in t1.iter_leaves()]
    n_leaves = len(leaf_list)



    for lf_count, lf in enumerate(t1.iter_leaves()):
        lf.branch_vertical_margin = 0

        first_column = True
        for col, value in enumerate(column_labels):
            if lf_count == 0:
                    # add labels
                    n = TextFace(' %s ' % str(value))
                    n.vt_align = 2
                    n.hz_align = 2
                    n.rotation= 270
                    n.margin_top = 2
                    n.margin_right = 2
                    n.margin_left = 2
                    n.margin_bottom = 6

                    from ete3 import NodeStyle
                    n.inner_background.color = "white"
                    n.opacity = 1.
                    tss.aligned_header.add_face(n, col)
                    #lf.add_face(n, col, position="aligned")

                    #nstyle = NodeStyle()
                    #nstyle["fgcolor"] = "red"
                    #nstyle["size"] = 15
                    #lf.set_style(nstyle)

            if first_column and not reference_column and highlight_first_column:
                print("highlight of the first column only (red)")
                # highlight of the first column only (red)
                try:

                    n = TextFace(' %s ' % str(group2taxon2count[value][lf.name]))

                except:
                    n = TextFace(' - ')
                    #n = RectFace(width=10, height=10, fgcolor="red", bgcolor="blue", label='-')
                    group2taxon2count[value][lf.name] = 0
                ref_data = str(value)
                n.margin_top = 2
                n.margin_right = 2
                n.margin_left = 30
                n.margin_bottom = 2

                if group2taxon2count[value][lf.name] > 0 and group2taxon2count[value][lf.name] != '-':
                    n.inner_background.color = "#FA5858"
                else:
                    n.inner_background.color = 'white'
                first_column = False
            else:
                print("no not highlight of the first column only (red)")
                if not taxon2group2value:
                    if value not in group2taxon2count:

                        #print 'not-------------------', value

                        n = TextFace(' - ')
                    else:
                        if lf.name not in group2taxon2count[value]:
                            n = TextFace(' - ')
                        else:
                            if group2taxon2count[value][lf.name] == '-':
                                n = TextFace(' - ')
                            else:
                                # if identity scale: 2 digit format
                                local_label = str(group2taxon2count[value][lf.name])
                                if not identity_scale:# and not column_scale:
                                    if as_float:
                                        local_label = "%.2f" % group2taxon2count[value][lf.name]
                                    else:
                                        try:
                                            local_label = "%s" % int(group2taxon2count[value][lf.name])
                                        except:
                                            local_label = "%s" % group2taxon2count[value][lf.name]
                                else:
                                    if round(float(group2taxon2count[value][lf.name]), 2) < 100:
                                        if type(group2taxon2count[value][lf.name]) != int:
                                            local_label = "%.2f" % round(group2taxon2count[value][lf.name], 2)
                                        else:
                                            local_label = "%s" % group2taxon2count[value][lf.name]
                                    else:
                                        if type(group2taxon2count[value][lf.name]) != int:
                                            local_label = "%.1f" % round(group2taxon2count[value][lf.name], 2)
                                        else:
                                            local_label = "%s" % group2taxon2count[value][lf.name]
                                if show_labels:
                                    try:
                                        if round(group2taxon2count[value][lf.name], 2) < 100 and column_scale:
                                            if round(group2taxon2count[value][lf.name], 2) < 10:
                                                n = TextFace('  %s  ' % local_label)
                                            else:
                                                n = TextFace(' %s  ' % local_label)
                                        else:
                                            n = TextFace('%s ' % local_label)

                                    # labels are not floats
                                    except TypeError:
                                        n = TextFace(' %s ' % group2taxon2count[value][lf.name])
                                else:
                                    n = TextFace(' - ')

                    '''
                    except:

                        if not reference_taxon_dico:
                            n = TextFace(' - ')

                        # if doctionnary locus2taxon report which taxon is the reference
                        # color the cell in blue
                        else:
                            if str(reference_taxon_dico[value]) == lf.name and show_labels:
                                n = TextFace(' 100.0 ')
                                n.fgcolor = '#58ACFA'
                            else:
                                n = TextFace(' - ')
                    '''
                    if show_labels:
                        n.margin_top = 2
                        n.margin_right = 2
                        n.margin_left = 2
                        n.margin_bottom = 2
                    else:
                        n.margin_top = 0
                        n.margin_right = 0
                        n.margin_left = 0
                        n.margin_bottom = 5
                    try:
                        count = group2taxon2count[value][lf.name]
                    except:
                        count = 0
                    #print value, lf.name
                    #print taxon2group2count[value][lf.name]
                    if count > 0 and group2taxon2count[value][lf.name] != '-':
                        if not reference_column:
                            if lf.name != str(reference_taxon):
                                # clor given a continuous scale
                                if identity_scale:
                                    n.inner_background.color = rgb2hex(m.to_rgba(float(count)))
                                    if not show_labels:
                                        n.fgcolor = rgb2hex(m.to_rgba(float(count)))
                                elif column_scale:
                                    n.inner_background.color = rgb2hex(column2scale[value].to_rgba(float(count)))
                                else:
                                    n.inner_background.color = "#58ACFA"
                            else:
                                n.inner_background.color = '#58ACFA'
                                n.fgcolor = "white"

                        else:
                            if lf.name == str(reference_taxon) or col == reference_column:
                                n.inner_background.color = "#FA5858"
                                n.fgcolor = "white"
                            else:
                                n.inner_background.color = "#58ACFA"
                    else:

                        if not reference_taxon_dico:
                            n.inner_background.color = 'white'
                        # if doctionnary locus2taxon report which taxon is the reference
                        # color the cell in blue
                        else:
                            if str(reference_taxon_dico[value]) == lf.name:
                                n.inner_background.color = '#58ACFA'
                            else:
                                n.inner_background.color = 'white'
                else:
                    # taxon2group2value set
                    # use it to color background
                    try:
                        n = TextFace(' %s ' % str(group2taxon2count[value][lf.name]))
                    except:
                        group2taxon2count[value][lf.name] = 0
                        n = TextFace(' - ')
                    n.margin_top = 2
                    n.margin_right = 2
                    n.margin_left = 2
                    n.margin_bottom = 2
                    if group2taxon2count[value][lf.name] > 0 and group2taxon2count[value][lf.name] != '-':
                        if not reference_column:
                            #print 'no ref column'
                            if not reference_taxon:
                                if lf.name != str(reference_taxon):
                                    #print 'group', value
                                    #print 'ref data', ref_data
                                    #print taxon2group2value[str(lf.name)][value]
                                    try:
                                        if ref_data not in taxon2group2value[lf.name][value]:
                                            n.inner_background.color = "#9FF781" #58ACFA
                                        else:
                                            n.inner_background.color = "#F78181" ##F6D8CE
                                    except:
                                        n.inner_background.color = "#9FF781"
                                else:
                                    n.inner_background.color = "#FA5858"
                            else:
                                #print 'group', lf.name
                                if ref_data not in taxon2group2value[lf.name][value]:
                                    n.inner_background.color = "#FA5858"
                                else:
                                    n.inner_background.color = "#FA5858"
                        else:
                            #print 'ref column'
                            if lf.name == str(reference_taxon) or col == reference_column:
                                n.inner_background.color = "#FA5858"
                                n.fgcolor = "white"

                            else:
                                n.inner_background.color = "#58ACFA"

                    else:
                        n.inner_background.color = 'white'
            if rotate:
               n.rotation= 270
            lf.add_face(n, col, position="aligned")

        #lf.name = taxon_id2organism_name[lf.name]
        n = TextFace(taxon_id2organism_name[lf.name], fgcolor = "black", fsize = 12, fstyle = 'italic')
        lf.add_face(n, 0)
        head=False

    return t1, tss

def multiple_profiles_heatmap_nobiodb(column_labels,
                              group2taxon2count,
                              taxon_id2organism_name,
                              reference_taxon=False,
                              reference_column = False,
                              taxon2group2value=False,
                              highlight_first_column=False,
                              identity_scale=False,
                              column_scale=False,
                              show_labels=True,
                              tree=False,
                              as_float=False,
                              rotate=False):

    '''

    ATTENTION: dico inverse: il faut group2taxon2values

    :param biodb:
    :param column_labels:
    :param taxon2group2count:
    :param reference_taxon:
    :param reference_column:
    :param taxon2group2value: color or not background based on
    :param highlight_first_column:
    :return:
    '''

    if isinstance(reference_taxon, dict):
        import copy
        reference_taxon_dico = copy.copy(reference_taxon)

        reference_taxon = False
    else:
        reference_taxon_dico = False


    from chlamdb.biosqldb import manipulate_biosqldb

    if identity_scale:
        import matplotlib.cm as cm
        from matplotlib.colors import rgb2hex
        import matplotlib as mpl

        norm = mpl.colors.Normalize(vmin=0, vmax=100)
        cmap = cm.OrRd
        m = cm.ScalarMappable(norm=norm, cmap=cmap)

    elif column_scale:
        import matplotlib.cm as cm
        from matplotlib.colors import rgb2hex
        import matplotlib as mpl
        column2scale = {}
        for column in column_labels:
            values = [float(i) for i in group2taxon2count[column].values()]
            #print values, column
            norm = mpl.colors.Normalize(vmin=min(values), vmax=max(values))
            cmap = cm.OrRd
            m = cm.ScalarMappable(norm=norm, cmap=cmap)
            column2scale[column] = m

    if not isinstance(tree, Tree):
        t1 = Tree(tree)

    else:
        t1 = tree
    tss = TreeStyle()
    tss.show_branch_support = False
    tss.draw_guiding_lines = True
    tss.guiding_lines_color = "gray"
    tss.show_leaf_name = False

    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    t1.ladderize()



    head = True
    leaf_list = [i for i in t1.iter_leaves()]
    n_leaves = len(leaf_list)

    for lf_count, lf in enumerate(t1.iter_leaves()):
        lf.branch_vertical_margin = 0

        first_column = True
        for col, value in enumerate(column_labels):
            if lf_count == 0:
                    # add labels
                    n = TextFace(' %s ' % str(value))
                    n.vt_align = 2
                    n.hz_align = 2
                    n.rotation= 270
                    n.margin_top = 2
                    n.margin_right = 2
                    n.margin_left = 2
                    n.margin_bottom = 6

                    from ete3 import NodeStyle
                    n.inner_background.color = "white"
                    n.opacity = 1.
                    tss.aligned_header.add_face(n, col)
                    #lf.add_face(n, col, position="aligned")

                    #nstyle = NodeStyle()
                    #nstyle["fgcolor"] = "red"
                    #nstyle["size"] = 15
                    #lf.set_style(nstyle)

            if first_column and not reference_column and highlight_first_column:
                # highlight of the first column only (red)
                try:

                    n = TextFace(' %s ' % str(group2taxon2count[value][lf.name]))

                except:
                    n = TextFace(' - ')
                    #n = RectFace(width=10, height=10, fgcolor="red", bgcolor="blue", label='-')
                    group2taxon2count[value][lf.name] = 0
                ref_data = str(value)
                n.margin_top = 2
                n.margin_right = 2
                n.margin_left = 30
                n.margin_bottom = 2

                if group2taxon2count[value][lf.name] > 0 and group2taxon2count[value][lf.name] != '-':
                    n.inner_background.color = "#FA5858"
                else:
                    n.inner_background.color = 'white'
                first_column = False
            else:
                if not taxon2group2value:
                    if value not in group2taxon2count:

                        #print 'not-------------------', value

                        n = TextFace(' - ')
                    else:
                        if lf.name not in group2taxon2count[value]:
                            n = TextFace(' - ')
                        else:
                            if group2taxon2count[value][lf.name] == '-':
                                n = TextFace(' - ')
                            else:
                                # if identity scale: 2 digit format
                                local_label = str(group2taxon2count[value][lf.name])
                                if not identity_scale:# and not column_scale:
                                    if as_float:
                                        local_label = "%.2f" % group2taxon2count[value][lf.name]
                                    else:
                                        try:
                                            local_label = "%s" % int(group2taxon2count[value][lf.name])
                                        except:
                                            local_label = "%s" % group2taxon2count[value][lf.name]
                                else:
                                    if round(float(group2taxon2count[value][lf.name]), 2) < 100:
                                        if type(group2taxon2count[value][lf.name]) != int:
                                            local_label = "%.2f" % round(group2taxon2count[value][lf.name], 2)
                                        else:
                                            local_label = "%s" % group2taxon2count[value][lf.name]
                                    else:
                                        if type(group2taxon2count[value][lf.name]) != int:
                                            local_label = "%.1f" % round(group2taxon2count[value][lf.name], 2)
                                        else:
                                            local_label = "%s" % group2taxon2count[value][lf.name]
                                if show_labels:
                                    try:
                                        if round(group2taxon2count[value][lf.name], 2) < 100 and column_scale:
                                            if round(group2taxon2count[value][lf.name], 2) < 10:
                                                n = TextFace('  %s  ' % local_label)
                                            else:
                                                n = TextFace(' %s  ' % local_label)
                                        else:
                                            n = TextFace('%s ' % local_label)

                                    # labels are not floats
                                    except TypeError:
                                        n = TextFace(' %s ' % group2taxon2count[value][lf.name])
                                else:
                                    n = TextFace(' - ')

                    '''
                    except:

                        if not reference_taxon_dico:
                            n = TextFace(' - ')

                        # if doctionnary locus2taxon report which taxon is the reference
                        # color the cell in blue
                        else:
                            if str(reference_taxon_dico[value]) == lf.name and show_labels:
                                n = TextFace(' 100.0 ')
                                n.fgcolor = '#58ACFA'
                            else:
                                n = TextFace(' - ')
                    '''
                    if show_labels:
                        n.margin_top = 2
                        n.margin_right = 2
                        n.margin_left = 2
                        n.margin_bottom = 2
                    else:
                        n.margin_top = 0
                        n.margin_right = 0
                        n.margin_left = 0
                        n.margin_bottom = 5
                    try:
                        count = group2taxon2count[value][lf.name]
                    except:
                        count = 0
                    #print value, lf.name
                    #print taxon2group2count[value][lf.name]
                    if count > 0 and group2taxon2count[value][lf.name] != '-':
                        if not reference_column:
                            if lf.name != str(reference_taxon):
                                # clor given a continuous scale
                                if identity_scale:
                                    n.inner_background.color = rgb2hex(m.to_rgba(float(count)))
                                    if not show_labels:
                                        n.fgcolor = rgb2hex(m.to_rgba(float(count)))
                                elif column_scale:
                                    n.inner_background.color = rgb2hex(column2scale[value].to_rgba(float(count)))
                                else:
                                    n.inner_background.color = "#58ACFA"
                            else:
                                n.inner_background.color = '#58ACFA'
                                n.fgcolor = "white"

                        else:
                            if lf.name == str(reference_taxon) or col == reference_column:
                                n.inner_background.color = "#FA5858"
                                n.fgcolor = "white"
                            else:
                                n.inner_background.color = "#58ACFA"
                    else:

                        if not reference_taxon_dico:
                            n.inner_background.color = 'white'
                        # if doctionnary locus2taxon report which taxon is the reference
                        # color the cell in blue
                        else:
                            if str(reference_taxon_dico[value]) == lf.name:
                                n.inner_background.color = '#58ACFA'
                            else:
                                n.inner_background.color = 'white'
                else:
                    try:
                        n = TextFace(' %s ' % str(group2taxon2count[value][lf.name]))
                    except:
                        group2taxon2count[value][lf.name] = 0
                        n = TextFace(' - ')
                    n.margin_top = 2
                    n.margin_right = 2
                    n.margin_left = 2
                    n.margin_bottom = 2
                    if group2taxon2count[value][lf.name] > 0 and group2taxon2count[value][lf.name] != '-':
                        if not reference_column:
                            #print 'no ref column'
                            if not reference_taxon:
                                if lf.name != str(reference_taxon):
                                    #print 'group', value
                                    #print 'ref data', ref_data
                                    #print taxon2group2value[str(lf.name)][value]
                                    try:
                                        if ref_data not in taxon2group2value[int(lf.name)][value]:
                                            n.inner_background.color = "#9FF781" #58ACFA
                                        else:
                                            n.inner_background.color = "#F78181" ##F6D8CE
                                    except:
                                        n.inner_background.color = "#9FF781"
                                else:
                                    n.inner_background.color = "#FA5858"
                            else:
                                #print 'group', lf.name
                                if ref_data not in taxon2group2value[lf.name][value]:
                                    n.inner_background.color = "#FA5858"
                                else:
                                    n.inner_background.color = "#FA5858"
                        else:
                            #print 'ref column'
                            if lf.name == str(reference_taxon) or col == reference_column:
                                n.inner_background.color = "#FA5858"
                                n.fgcolor = "white"

                            else:
                                n.inner_background.color = "#58ACFA"

                    else:
                        n.inner_background.color = 'white'
            if rotate:
               n.rotation= 270
            lf.add_face(n, col, position="aligned")

        #lf.name = taxon_id2organism_name[lf.name]
        n = TextFace(taxon_id2organism_name[lf.name], fgcolor = "black", fsize = 12, fstyle = 'italic')
        lf.add_face(n, 0)
        head=False

    return t1, tss


def multiple_orthogroup_heatmap(biodb, reference_orthogroup, max_distance=2.2):

    from chlamdb.biosqldb import manipulate_biosqldb
    import pandas
    import matplotlib.cm as cm
    from matplotlib.colors import rgb2hex
    import matplotlib as mpl

    server, db = manipulate_biosqldb.load_db(biodb)

    #queries = ['selv']

    sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

    tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]
    #print tree
    t1 = Tree(tree)
    #t.populate(8)
    # Calculate the midpoint node
    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    t1.ladderize()

    taxon_id2organism_name = manipulate_biosqldb.taxon_id2genome_description(server, biodb)


    sql = 'select * from phylogenetic_profiles_euclidian_distance_%s' \
          ' where group_1="%s" or group_2="%s" and euclidian_dist <=%s;' % (biodb,
                                                                          reference_orthogroup,
                                                                          reference_orthogroup,
                                                                            max_distance)

    data = list(server.adaptor.execute_and_fetchall(sql,))

    data_frame = pandas.DataFrame(data)
    sorted_data_frame = data_frame.sort(2)


    ordered_orthogroups = [reference_orthogroup] + list(sorted_data_frame[sorted_data_frame.columns[0]])

    l = sorted(set(sorted_data_frame[sorted_data_frame.columns[2]]))

    colmap = dict(zip(l,range(len(l))[::-1]))
    #print dict(colmap)
    norm = mpl.colors.Normalize(vmin=-2, vmax=len(l))

    cmap = cm.OrRd
    cmap_blue = cm.Blues
    #m2 = cm.ScalarMappable(norm=norm, cmap=cmap)
    m2 = cm.ScalarMappable(norm=norm, cmap=cmap_blue)

    orthogroup2distance = {}
    distances = []

    for one_pair in sorted_data_frame.itertuples(index=False):
        distances.append(one_pair[2])
        if one_pair[0] == one_pair[1]:
            orthogroup2distance[one_pair[0]] = one_pair[2]
        elif one_pair[0] == reference_orthogroup:
            orthogroup2distance[one_pair[1]] = one_pair[2]
        elif one_pair[1] == reference_orthogroup:
            orthogroup2distance[one_pair[0]] = one_pair[2]
        else:
            raise 'Error: unexpected combination of groups'
    ordered_distances = sorted(distances)

    if db_driver == 'mysql':
        sql = 'show columns from comparative_tables_orthology'
    elif db_driver == 'sqlite':
        sql = 'PRAGMA table_info(comparative_tables_orthology);'
        
    ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)][1:]

    #print 'taxons!', ordered_taxons

    ortho_sql = '"' + '","'.join(orthogroup2distance.keys()) + '"' + ',"%s"' % reference_orthogroup

    sql = 'select * from comparative_tables_orthology where orthogroup in (%s)' % (biodb, ortho_sql)

    profile_tuples = list(server.adaptor.execute_and_fetchall(sql,))
    #print "profile_tuples", profile_tuples
    #profile_data = pandas.DataFrame(profile_tuples)
    #print "profile_data", profile_data

    taxon2group2n_homologs = {}
    #all_orthogroups = list(profile_data[profile_data.columns[0]])
    #import sys
    #sys.exit()
    #print "ordered_taxons", len(ordered_taxons), ordered_taxons
    for i, tuple in enumerate(profile_tuples):
        # get position of the group based on score
        # get colum of taxon i
        taxon2group2n_homologs[tuple[0]] = {}
        for i, taxon in enumerate(ordered_taxons):
            taxon2group2n_homologs[tuple[0]][str(taxon)] = tuple[i+1]
    #print taxon2group2n_homologs
    # and set it as tree outgroup
    head = True
    for lf in t1.iter_leaves():
        #lf.add_face(AttrFace("name", fsize=20), 0, position="branch-right")
        lf.branch_vertical_margin = 0
        #data = [random.randint(0,2) for x in xrange(3)]

        for col, value in enumerate(ordered_orthogroups):
            #print 'value', value
            if head:

                    #'first row, print gene names'
                    #print 'ok!'
                    n = TextFace(' %s ' % str(value))
                    n.rotation= 270
                    n.margin_top = 4
                    n.margin_right = 4
                    n.margin_left = 4
                    n.margin_bottom = 4
                    if value == reference_orthogroup:
                        n.inner_background.color = "red"
                    else:
                        n.inner_background.color = "white"
                    n.opacity = 1.
                    lf.add_face(n, col, position="aligned")



            #if float(identity_value) >70:
            #    if str(identity_value) == '100.00' or str(identity_value) == '100.0':
            #        identity_value = '100'
            #    else:
            #        identity_value = str(round(float(identity_value), 1))
            n = TextFace(' %s ' % str(taxon2group2n_homologs[value][lf.name]))
            n.margin_top = 4
            n.margin_right = 4
            n.margin_left = 4
            n.margin_bottom = 4
            if taxon2group2n_homologs[value][lf.name] >0:
                if value == reference_orthogroup:
                    n.inner_background.color = "red"
                else:

                    n.inner_background.color = rgb2hex(m2.to_rgba(float(colmap[orthogroup2distance[value]])))

            else:
                n.inner_background.color = 'white'


            #n.inner_background.color = rgb2hex(m.to_rgba(float(identity_value)))
            lf.add_face(n, col, position="aligned")
        lf.name = taxon_id2organism_name[lf.name]
        head=False


    # , tree_style=ts
    t1.render("test.png", dpi=800, h=400)
    #t1.write(format=0, outfile="new_tree.nw")

def get_pfam_data(orthogroup, biodb, aa_alignment=False):

    '''

    :param orthogroup: target orthogroup id
    :param biodb: biodb to which the orthogroup belong
    :param aa_alignment: get aa align as well (optional)
    :return: a dictionnary with locus tags as keys, with all data necessary to draw a PFAM phylogeny figure with the
             function "draw_pfam_tree"
    '''

    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    # inefficient query
    '''
    sql = 'select A.locus_tag, B.start, B.stop, A.organism, A.sequence_length, B.signature_accession, B.signature_description, A.taxon_id ' \
          ' from (select taxon_id, orthogroup,locus_tag, protein_id, organism, length(translation) as sequence_length from orthology_detail ' \
          ' where orthogroup="%s" ) A ' \
          ' left join (select * from interpro where orthogroup="%s" and analysis="Pfam") B ' \
          ' on A.locus_tag=B.locus_tag;' % (biodb, orthogroup, biodb, orthogroup)
    '''
    # correct with filter on join
    sql = 'select t1.locus_tag,t2.start,t2.stop,t1.organism,length(t1.translation),t2.signature_accession,t2.signature_description,t1.taxon_id ' \
          ' from orthology_detail t1 ' \
          ' left join interpro t2 on t1.seqfeature_id=t2.seqfeature_id and t2.orthogroup="%s" and t2.analysis="Pfam" ' \
          ' where t1.orthogroup="%s";' % (orthogroup, orthogroup)

    #print(sql)
    data = server.adaptor.execute_and_fetchall(sql,)
    #print("OK")
    locus2aa_seq = {}
    # getting aa alignment
    if aa_alignment:
        from Bio import AlignIO
        alignment = AlignIO.read(aa_alignment, "fasta")
        for record in alignment:
            locus2aa_seq[record.id] = str(record.seq)

    locus2data = {}
    for one_locus in data:
        if one_locus[0] not in locus2data:
            if aa_alignment:
                locus2data[one_locus[0]] = [list(one_locus[1:len(one_locus)]) + [locus2aa_seq[one_locus[0]]]]
            else:
                locus2data[one_locus[0]] = [list(one_locus[1:len(one_locus)])]
        else:
            if aa_alignment:
                locus2data[one_locus[0]].append(list(one_locus[1:len(one_locus)])+ [locus2aa_seq[one_locus[0]]])
            else:
                locus2data[one_locus[0]].append(list(one_locus[1:len(one_locus)]))
    return locus2data


def get_TM_data(biodb,
                orthogroup=False,
                aa_alignment=False,
                signal_peptide=False):
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    if orthogroup:
        #print("orthogroup!")
        sql = 'select locus_tag, start, stop, organism, sequence_length, signature_accession, signature_description  ' \
          ' from interpro as t2 where orthogroup="%s" and analysis="Phobius" and signature_accession="TRANSMEMBRANE"' % (orthogroup)
        if db_driver == 'mysql':
            sql2 = 'select locus_tag, char_length(translation), organism from orthology_detail where orthogroup="%s";' % (orthogroup)
        if db_driver == 'sqlite':
            sql2 = 'select locus_tag, length(translation), organism from orthology_detail where orthogroup="%s";' % (orthogroup)
        sql_signalp = 'select locus_tag, t2.start, t2.stop, t6.description, sequence_length, signature_accession, signature_description' \
                       ' from annotation_seqfeature_id2locus t1' \
                       ' inner join interpro_interpro t2 on t1.seqfeature_id=t2.seqfeature_id' \
                       ' inner join interpro_signature t3 on t2.signature_id=t3.signature_id' \
                       ' inner join orthology_seqfeature_id2orthogroup t4 on t1.seqfeature_id=t4.seqfeature_id' \
                       ' inner join orthology_orthogroup t5 on t4.orthogroup_id=t5.orthogroup_id' \
                       ' inner join bioentry t6 on t1.bioentry_id=t6.bioentry_id' \
                       ' where t5.orthogroup_name="%s" and signature_description="Signal peptide region";' % (orthogroup)

        locus2seq_length = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
    else:
        sql = 'select locus_tag, start, stop, organism, sequence_length, signature_accession, signature_description  ' \
          ' from interpro as t2 where analysis="Phobius" and signature_accession="TRANSMEMBRANE"'
    #print(sql)
    data = server.adaptor.execute_and_fetchall(sql,)

    if signal_peptide:
        data2 = server.adaptor.execute_and_fetchall(sql_signalp,)
        data += data2

    locus2aa_seq = {}
    # getting aa alignment
    if aa_alignment:
        from Bio import AlignIO
        alignment = AlignIO.read(aa_alignment, "fasta")
        for record in alignment:
            locus2aa_seq[record.id] = str(record.seq)

    protein2data = {}
    for one_locus in data:
        if one_locus[0] not in protein2data:
            if aa_alignment:
                protein2data[one_locus[0]] = [list(one_locus[1:len(one_locus)])+[locus2aa_seq[one_locus[0]]]]
            else:
                protein2data[one_locus[0]] = [list(one_locus[1:len(one_locus)])]
        else:
            if aa_alignment:
                protein2data[one_locus[0]].append(list(one_locus[1:len(one_locus)])+[locus2aa_seq[one_locus[0]]])
            else:
                protein2data[one_locus[0]].append(one_locus[1:len(one_locus)])
    if orthogroup:
        for locus in locus2seq_length:
            if locus not in protein2data:
                protein2data[locus] = list(locus2seq_length[locus])
    return protein2data


seq = "LHGRISQQVEQSRSQVQAIGEKVSLAQAKIEKIKGSKKAIKVFSSAKYPAPERLQEYGSIFTDAQDPGLQRRPRHRIQSKQRPLDERALQEKLKDFPVCVSTKPEPEDDAEEGLGGLPSNISSVSSLLLFNTTENLYKKYVFLDPLAGAVTKTHVMLGAETEEKLFDAPLSISKREQLEQQVPENYFYVPDLGQVPEIDVPSYLPDLPGIANDLMYIADLGPGIAPSAPGTIPELPTFHTEVAEPLKVGELGSGMGAGPGTPAHTPSSLDTPHFVFQTYKMGAPPLPPSTAAPVGQGARQDDSSSSASPSVQGAPREVVDPSGGWATLLESIRQAGGIGKAKLRSMKERKLEKQQQKEQEQVRATSQGGHLMSDLFNKLVMRRKGISGKGPGAGDGPGGAFARVSDSIPPLPPPQQPQAEDEDDWES"
motifs = [
    # seq.start, seq.end, shape, width, height, fgcolor, bgcolor
    [0, len(seq), "-", 1, 1, "black", None, None],
    [10, 50, "[]", 40, 10, "black", "PaleGreen", "arial|8|red|TLC"],
    [150, 280, "()", 130, 10, "black", "PaleGreen", "arial|8|red|TLC"],
    [80, 120, "()", 40, 10, "black", "grey", "arial|3|black|TLC"],
    #[120, 190, "^", None, 14, "black", "yellow", None],
    [191, 200, "v", None, 12, "black", "rgradient:orange", None],
    [185, 190, "o", None, 12, "black", "brown", None],
    [198, 200, "<>", None, 15, "black", "rgradient:gold", None],
    [210, 420, "compactseq", 2, 10, None, None, None],
    [1, 40, "seq", 10, 10, None, None, None],
    [310, 320, "<>", None, 30, "black", "rgradient:black", None],
    [0, len(seq), "-", None, 10, "black", None, None],
    [11, 30, "()", None, 20, "blue", "blue", None],
    [300, 310, "()", None, 40, "green", "green", None],
    ]

def layout(node):
    if node.is_leaf():
        motifs_seqface = SeqMotifFace(seq, motifs, scale_factor=1)

        add_face_to_node(motifs_seqface, node, 1, position="aligned")


def organism2color(locus2data, taxon_id2family=False):


    if not taxon_id2family:
        organism_list = []
        for locus in locus2data:
            # case in which we only hace seq length
            if type(locus2data[locus][0]) != list:
                #print locus2data[locus]
                #print '2!!!!'
                if locus2data[locus][1] not in organism_list:
                    organism_list.append(locus2data[locus][1])
            else:
                #print 'not 2!!!!!'
                if locus2data[locus][0][2] not in organism_list:
                    organism_list.append(locus2data[locus][0][2])
        colors = _get_colors(len(organism_list))
        #print organism_list, colors
        return dict(zip(organism_list, colors))
    else:
        family_list = []
        for locus in locus2data:
            taxon_id = locus2data[locus][0][-1]
            if taxon_id2family[str(taxon_id)] not in family_list:
                family_list.append(taxon_id2family[str(taxon_id)])
        colors = _get_colors(len(family_list))

        return dict(zip(family_list, colors))

def interpro_tsv2pfam_data(interpro_tsv_file,
                           locus2organism):
    '''

    interpro_tsv_file: interpro output file for locus of interest. Consider first colum as locus_tags (or protein IDs)
    locus2organism_name_table: Dictionnary to convert locus_tags to organism name for leaf labels for the tree.
    :return: a dictionarry woth locus_tags as keys with all data necessary to draw a PFAM phylogeny figure with the
             function "draw_pfam_tree"

    0 domain start
    1 domain stop
    2 organism name
    3 sequence_length
    4 signature_accession
    5 signature_description
    6 taxon_id
    7 sequence (optional)

    '''

    locus2data = {}
    with open(interpro_tsv_file, 'r') as f:
        for row in f:
            data = row.rstrip().split('\t')
            if len(data) != 13:
                #print (len(data))
                #print (data)
                raise IOError('Unexpected number of columns in the tsv file (expected 13 of the standard 13 columns format from 2018)')
            else:
                locus_tag = data[0]
                domain_start = int(data[6])
                domain_stop = int(data[7])
                organism_name = locus2organism[locus_tag]
                seq_length = int(data[2])
                signature_accession = data[4]
                signature_description = data[5]
                taxon_id = 0
                if locus_tag not in locus2data:
                    locus2data[data[0]] = [[domain_start,
                                            domain_stop,
                                            organism_name,
                                            seq_length,
                                            signature_accession,
                                            signature_description,
                                            taxon_id]]
                else:
                    locus2data[data[0]].append([domain_start,
                                            domain_stop,
                                            organism_name,
                                            seq_length,
                                            signature_accession,
                                            signature_description,
                                            taxon_id])
    #print locus2data
    return locus2data


def draw_pfam_tree(tree_name, locus2data,
                   locus2protein_id=False,
                   taxon_id2family=False,
                   color_by_family=False):

    '''
    format: locus2data
    data: list of lists with:

    0 domain start
    1 domain stop
    2 organism name
    3 sequence_length
    4 signature_accession
    5 signature_description
    6 taxon_id
    7 sequence (optional)
    '''

    # from the original set
    from ete3 import Tree, TreeStyle, faces, AttrFace
    t = Tree(tree_name)
    ts = TreeStyle()
    ts.draw_guiding_lines = True
    ts.guiding_lines_color = "gray"   
    #ts.show_leaf_name = False
    ts.show_branch_support = True
    #ts.layout_fn = layout    
    #t.populate(8)
    # Calculate the midpoint node
    #print (t)
    R = t.get_midpoint_outgroup()
    #print (R)
    # and set it as tree outgroup
    t.set_outgroup(R)

    color_dico = organism2color(locus2data, taxon_id2family)

    leaf_number = 0
    for i, l in enumerate(t.iter_leaves()):
        # top leaf: add headers
        if i == 0:
            
            n = TextFace('Locus tag')
            n.margin_top = 1
            n.margin_right = 1
            n.margin_left = 20
            n.margin_bottom = 1
            n.inner_background.color = "white"
            n.opacity = 1.
            n.rotation = -25
            #lf.add_face(n, 7, position="aligned")
            ts.aligned_header.add_face(n, 0)
 
            n = TextFace('Pfam domain(s)')
            n.margin_top = 1
            n.margin_right = 1
            n.margin_left = 20
            n.margin_bottom = 1
            n.inner_background.color = "white"
            n.opacity = 1.
            n.rotation = -25
            #lf.add_face(n, 7, position="aligned")
            ts.aligned_header.add_face(n, 1)
            
        leaf_number+=1
        if locus2protein_id:
            #protein_id = locus2protein_id[str(l)[3:len(str(l))]]
            data = locus2data[l]
        else:
            data = locus2data[str(l)[3:len(str(l))]]
        # [[None, None, 'Chlamydia pecorum (Chlamydophila pecorum) L1', 54, None, None, 13]]
        # [[147, 258, 'Chlamydia pecorum E58', 303, 'PF13091', 'PLD-like domain', 4]]
        seq_motifs = []

        #print 'data!', data

        #l.img_style['hz_line_type'] = 0
        #l.img_style['size'] = 10

        for motif in data:
            # check if motif is not None
            if motif[0]:
                '''
                #print 'motif', motif[0], motif[1]
                seq = data[0][-1]

                count=0
                for n, letter in enumerate(seq):
                    if letter != '-':
                        count+=1
                    if count == motif[0]:
                        align_start = n+1
                    if count == motif[1]:
                        align_end = n+1
                '''



                seq_motifs.append([motif[0], motif[1], "[]", None, 10, "black", "PaleGreen", "arial|8|red|%s" % motif[4]])
        # check if alignment is available or not
        if len(data[0]) == 8:
            #print 'tata', len(data[0][-1])
            seqFace = SeqMotifFace(data[0][-1],
                                   motifs=seq_motifs,
                                   width=10,
                                   height=12,
                                   gap_format='-',
                                   seq_format='-',
                                   gapcolor='white') #seq, seq_motifs, scale_factor=1, intermotif_format=None) #intermotif_format="seq", seqtail_format=None, seqtype='aa')
                                   # seqtail_format="-",
        else:
            # no alignment available
            seqFace = SeqMotifFace(data[0][3]*'N',
                                   motifs=seq_motifs,
                                   width=10,
                                   height=12,
                                   gap_format='-',
                                   seq_format='-',
                                   gapcolor='white') #seq, seq_motifs, scale_factor=1, intermotif_format=None) #intermotif_format="seq", seqtail_format=None, seqtype='aa')
                                   #seqtail_format="-",

        seqFace.margin_bottom = 2
        seqFace.margin_top = 2
        seqFace.opacity = 1.0

        l.add_face(seqFace, column=1, position="aligned")
        locus = TextFace(str(l)[3:len(str(l))])
        l.name = data[0][2]

        if not taxon_id2family:
            if color_by_family:
                l.img_style['fgcolor'] = color_dico[data[0][2]]
                l.img_style['size'] = 6

        else:
            taxon_id = data[0][-1]
            family = taxon_id2family[str(taxon_id)]
            col = color_dico[family]

            ff = AttrFace("name", fsize=12)
            #ff.background.color = 'red'
            ff.fgcolor = col
            if color_by_family:
                #print 'color----------'
                l.add_face(ff, column=0)

        locus.margin_right = 10
        locus.margin_left = 22
        locus.margin_bottom = 0
        l.add_face(locus, column=0, position="aligned")


    return t, ts, leaf_number


def draw_TM_tree(tree_name, locus2data):
    # from the original set
    from ete3 import Tree, TreeStyle, faces, AttrFace
    
    t = Tree(tree_name)
    ts = TreeStyle()
    ts.draw_guiding_lines = True
    ts.guiding_lines_color = "gray"    
    #t.populate(8)
    # Calculate the midpoint node
    R = t.get_midpoint_outgroup()
    # and set it as tree outgroup
    t.set_outgroup(R)

    #print("draw_TM_tree--------------", locus2data)
    color_dico = organism2color(locus2data)

    for leaf_number, l in enumerate(t.iter_leaves()):


        if leaf_number == 0:
            
            n = TextFace('Locus tag')
            n.margin_top = 1
            n.margin_right = 1
            n.margin_left = 20
            n.margin_bottom = 1
            n.inner_background.color = "white"
            n.opacity = 1.
            n.rotation = -25
            ts.aligned_header.add_face(n, 0)
 
            n = TextFace('SP/TM domain(s)')
            n.margin_top = 1
            n.margin_right = 1
            n.margin_left = 20
            n.margin_bottom = 1
            n.inner_background.color = "white"
            n.opacity = 1.
            n.rotation = -25
            ts.aligned_header.add_face(n, 1)

        locus_name = str(l)[3:len(str(l))]
        locus_name = locus_name.split('|')[0]
        #print(locus_name, locus_name in locus2data)
        try:
            data = locus2data[locus_name]
            seq_motifs = []
            l.img_style['fgcolor'] = color_dico[data[0][2]]
            l.img_style['hz_line_type'] = 0
            l.img_style['size'] = 6
        except:
            seq_motifs = []
            l.img_style['fgcolor'] = color_dico[data[1]]
            l.img_style['hz_line_type'] = 0
            l.img_style['size'] = 6
        # case in which we have more than seq length
        #print("len data", len(data), data)
        if isinstance(data[0], list):
            for motif in data:
                #print(motif[-1])
                if motif[-2] != "SIGNAL_PEPTIDE":
                    seq_motifs.append([motif[0], motif[1], "()", None, 10, "black", "PaleGreen", "arial|8|red|"])
                else:
                    seq_motifs.append([motif[0], motif[1], "[]", None, 10, "black", "red", "arial|8|red|"])
            #print("motifs:", seq_motifs)
        if isinstance(data[0], int):
            seqFace = SeqMotifFace(data[0]*'N',
                                    motifs=[],
                                    width=10,
                                    height=12,
                                    gap_format='-',
                                    seq_format='-',
                                    gapcolor='white')
        else:
            if isinstance(data[0], int):
                seqFace = SeqMotifFace(data[0][-1],
                                       motifs=seq_motifs,
                                       width=10,
                                       height=12,
                                       gap_format='-',
                                       seq_format='-',
                                       gapcolor='white') #seq, seq_motifs, scale_factor=1, intermotif_format=None) #intermotif_format="seq", seqtail_format=None, seqtype='aa')


            else:
                seqFace = SeqMotifFace(data[0][3]*'N',
                                       motifs=seq_motifs,
                                       width=10,
                                       height=12,
                                       gap_format='-',
                                       seq_format='-',
                                       gapcolor='white')

        seqFace.margin_bottom = 2
        seqFace.margin_top = 2
        seqFace.opacity = 1.0

        l.add_face(seqFace, column=1, position="aligned")
        locus = TextFace(str(l)[3:len(str(l))])
        try:
            l.name = data[0][2]
        except:
            l.name = data[1]
        locus.margin_right = 10
        locus.margin_left = 15
        locus.margin_bottom = 0
        l.add_face(locus, column=0, position="aligned")

    
    #ts.layout_fn = layout
    return t, ts, leaf_number+1


def draw_basic_tree(newick_tree, 
                    leaf2description):
    # from the original set
    from ete3 import Tree, TreeStyle, faces, AttrFace
    
    t = Tree(newick_tree)
    ts = TreeStyle()
    ts.draw_guiding_lines = True
    ts.guiding_lines_color = "gray"    
    #t.populate(8)
    # Calculate the midpoint node
    R = t.get_midpoint_outgroup()
    # and set it as tree outgroup
    t.set_outgroup(R)

    for leaf_number, l in enumerate(t.iter_leaves()):


        if leaf_number == 0:
            
            n = TextFace('Locus tag')
            n.margin_top = 1
            n.margin_right = 1
            n.margin_left = 20
            n.margin_bottom = 1
            n.inner_background.color = "white"
            n.opacity = 1.
            n.rotation = -25
            ts.aligned_header.add_face(n, 0)
 
        locus_name = str(l)[3:len(str(l))]
        locus_name = locus_name.split('|')[0]
        locus = TextFace(str(l)[3:len(str(l))])
        locus.margin_right = 10
        locus.margin_left = 15
        locus.margin_bottom = 0
        l.add_face(locus, column=0, position="aligned")

        l.name = leaf2description[l.name]

    return t, ts, leaf_number+1

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t",'--tree',type=str,help="newick tree")



    args = parser.parse_args()


    #locus2pfam_data = get_pfam_data("group_20", 'chlamydia_03_15')

    #t, ts = draw_pfam_tree(args.tree, locus2pfam_data)
    '''
    locus2TM_data = get_TM_data('chlamydia_03_15')
    t, ts = draw_TM_tree(args.tree, locus2TM_data)
    t.render("motifs.svg", w=1200, dpi=800, tree_style=ts)
    #t.show(tree_style=ts)
    '''
    multiple_COGs_heatmap("chlamydia_12_15", ['COG0835','COG2201','COG0840','COG2201','COG2208'])
