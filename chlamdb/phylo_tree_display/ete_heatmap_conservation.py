#!/usr/bin/env python

import sys
from random import sample
from random import randint
from ete3 import Tree, SeqMotifFace, TreeStyle, add_face_to_node, TextFace, StackedBarFace


import numpy as np
import colorsys
import matplotlib.colors as pcolors


def _get_colors(num_colors):
    colors=[]
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i/360.
        lightness = (50 + np.random.rand() * 10)/100.
        saturation = (90 + np.random.rand() * 10)/100.
        rgb = colorsys.hls_to_rgb(hue, lightness, saturation)
        colors.append(pcolors.rgb2hex(rgb))
    return colors




def organism2color(locus2data):

    organism_list = []
    for locus in locus2data:
        if locus2data[locus][0][2] not in organism_list:
            organism_list.append(locus2data[locus][0][2])
    colors = _get_colors(len(organism_list))

    return dict(zip(organism_list, colors))


def plot_heat_tree_V1(taxid2n, tree_file, genes, taxid2st=False, leaf_label_conversion_dico=False):
    '''
    Plot heatmap next to a tree. The order of the heatmap **MUST** be the same,
    as order of the leafs on the tree. The tree must be in the Newick format. If
    *output_file* is specified, then heat-tree will be rendered as a PNG,
    otherwise interactive browser will pop-up with your heat-tree.

    TODO ajouter en option la possibilite d'ajouter en option la valeur dans la cellule

    Parameters
    ----------


    tree_file: str
        Path to the tree file in Newick format. The leaf node labels should
        be the same as as row names in the heatmap file. E.g. row1, row2.

    output_file: str, optional
        If specified the heat-tree will be rendered in that file as a PNG image,
        otherwise interactive browser will pop-up. **N.B.** program will wait
        for you to exit the browser before continuing.
    '''

    t1 = Tree(tree_file)
    tss = TreeStyle()
    #t.populate(8)
    # Calculate the midpoint node
    R = t1.get_midpoint_outgroup()
    # and set it as tree outgroup
    t1.set_outgroup(R)    # To operate with numbers efficiently

    import matplotlib.cm as cm
    from matplotlib.colors import rgb2hex
    import matplotlib as mpl
    norm = mpl.colors.Normalize(vmin=0.8, vmax=1) # map2count[map[0]][0]
    cmap_blue = cm.Blues
    m2 = cm.ScalarMappable(norm=norm, cmap=cmap_blue)

    leaf_number = 0
    for lf in t1.iter_leaves():
        leaf_number+=1
        lf.branch_vertical_margin = 0

        try:
            data = taxid2n[str(lf.name)]
        except:
            data=[0]

        try:
            st = taxid2st[lf.name]
        except:
            st = False
            '''
            if "taxon2accession_list" not in locals():
                from chlamdb.biosqldb import manipulate_biosqldb
                server, db = manipulate_biosqldb.load_db("k_cosson_05_16")
                sql = 'select taxon_id, accession from bioentry where biodatabase_id=104'
                data_tax = server.adaptor.execute_and_fetchall(sql,)
                taxon2accession_list = {}
                for i in data_tax:
                    if i[0] not in taxon2accession_list:
                        taxon2accession_list[i[0]] = [i[1]]
                    else:
                       taxon2accession_list[i[0]].append(i[1])
            else:
                for taxon in taxon2accession_list:
                    if lf.name in taxon2accession_list[taxon]:
                        for accession in taxon2accession_list[taxon]:
                            print lf.name, accession
                            try:
                                st = taxid2st[accession]
                                data = taxid2n[accession]
                                print 'st ok!!', st
                                break
                            except:
                                continue
             '''


        if accession2description:
            try:
                lf.name = accession2description[lf.name]
            except:
                pass
        if st:
            lf.name = lf.name + ' (' + st + ')'
        else:
            pass
        for col, value in enumerate(data):

            if leaf_number==1:
                n = TextFace('%s' % (genes[col]), fsize=6)
                n.vt_align = 2
                n.hz_align = 2
                n.rotation= 270
                n.margin_top = 0
                n.margin_right = 0
                n.margin_left = 4
                n.margin_bottom = 0
                n.inner_background.color = "white"
                n.opacity = 1.
                tss.aligned_header.add_face(n, col)
                #lf.add_face(n, col, position="aligned")

            if value > 0:
                n = TextFace('  ')
                n.margin_top = 0
                n.margin_right = 0
                n.margin_left = 0
                n.margin_bottom = 0
                n.inner_background.color = rgb2hex(m2.to_rgba(float(value))) #'#140718' #"#81BEF7"
                n.opacity = 1.
                lf.add_face(n, col, position="aligned")

            else:
                n = TextFace('  ')
                n.margin_top = 0
                n.margin_right = 0
                n.margin_left = 0
                n.margin_bottom = 0
                n.inner_background.color = "white"
                n.opacity = 1.
                lf.add_face(n, col, position="aligned")

    return t1, leaf_number, tss

def plot_heat_tree_V0(heatmap_file, tree_file, output_file=None):
    '''
    Plot heatmap next to a tree. The order of the heatmap **MUST** be the same,
    as order of the leafs on the tree. The tree must be in the Newick format. If
    *output_file* is specified, then heat-tree will be rendered as a PNG,
    otherwise interactive browser will pop-up with your heat-tree.

    Parameters
    ----------
    heatmap_file: str
        Path to the heatmap file. The first row must have '#Names' as first
        element of the header.
            e.g. #Names, A, B, C, D
                row1, 2, 4, 0, 4
                row2, 4, 6, 2, -1

    tree_file: str
        Path to the tree file in Newick format. The leaf node labels should
        be the same as as row names in the heatmap file. E.g. row1, row2.

    output_file: str, optional
        If specified the heat-tree will be rendered in that file as a PNG image,
        otherwise interactive browser will pop-up. **N.B.** program will wait
        for you to exit the browser before continuing.
    '''
    import numpy

    from ete3.treeview.faces import add_face_to_node
    from ete3 import ClusterTree, TreeStyle, AttrFace, ProfileFace


    # To operate with numbers efficiently

    # Loads tree and array
    t = ClusterTree(tree_file, heatmap_file)
    t.ladderize()
    R = t.get_midpoint_outgroup()
    t.set_outgroup(R)
    # nodes are linked to the array table
    array =  t.arraytable

    # Calculates some stats on the matrix. Needed to establish the color
    # gradients.
    matrix_dist = [i for r in xrange(len(array.matrix))\
                   for i in array.matrix[r] if numpy.isfinite(i)]
    matrix_max = numpy.max(matrix_dist)
    matrix_min = numpy.min(matrix_dist)
    matrix_avg = matrix_min+((matrix_max-matrix_min)/2)

    # Creates a profile face that will represent node's profile as a
    # heatmap
    profileFace  = ProfileFace(1., 0., 0.5, 1000, 14, "heatmap",colorscheme=1)

    nameFace = AttrFace("name", fsize=8)
    # Creates my own layout function that uses previous faces
    def mylayout(node):
        # If node is a leaf
        if node.is_leaf():
            # And a line profile
            add_face_to_node(profileFace, node, 0, aligned=True)
            node.img_style["size"]=0
            add_face_to_node(nameFace, node, 1, aligned=True)

    # Use my layout to visualize the tree
    ts = TreeStyle()
    ts.layout_fn = mylayout
    t.render("test.svg",tree_style=ts)



    '''
    import random


    for lf in t1.iter_leaves():
        # Each leaf node must have a profile and a deviation vector, which will be based on your source matrix of values
        data = [random.randint(0,2) for x in xrange(2)]
        print data
        lf.add_features(profile = data)
        lf.img_style["size"]=0
        # if no std deviation for the vector values, just use 0
        lf.add_features(deviation = [0 for x in xrange(10)])

        # Add a ProfileFace to each leaf node (you probably want aligned position). Choose among the following types.

        ProFace =  ProfileFace(max_v=1, min_v=0.0, center_v=0.5, width=40, height=20, style='heatmap', colorscheme=0)
        ProFace.margin_bottom = 2
        ProFace.margin_top = 2
        ProFace.opacity = 0.5
        lf.add_face(ProFace, column=0, position="aligned")
        lf.add_face(TextFace("hola "), column=1)
        #lf.add_face(ProfileFace(max_v=1, min_v=0.0, center_v=0.5, width=200, height=40, style='heatmap', colorscheme=1), column=1, position="aligned")
        #lf.add_face(ProfileFace(max_v=1, min_v=0.0, center_v=0.5, width=200, height=40, style='heatmap', colorscheme=2), column=2, position="aligned")
        #lf.add_face(ProfileFace(max_v=1, min_v=0.0, center_v=0.5, width=200, height=40, style='lines', colorscheme=2), column=1, position="aligned")
        #lf.add_face(ProfileFace(max_v=1, min_v=0.0, center_v=0.5, width=200, height=40, style='bars', colorscheme=2), column=2, position="aligned")
        #lf.add_face(ProfileFace(max_v=1, min_v=0.0, center_v=0.5, width=200, height=40, style='cbars', colorscheme=2), column=3, position="aligned")
    '''







def plot_heat_tree(biodb, taxid2n, tree_file):
    '''
    Plot heatmap next to a tree. The order of the heatmap **MUST** be the same,
    as order of the leafs on the tree. The tree must be in the Newick format. If
    *output_file* is specified, then heat-tree will be rendered as a PNG,
    otherwise interactive browser will pop-up with your heat-tree.

    Parameters
    ----------
    heatmap_file: str
        Path to the heatmap file. The first row must have '#Names' as first
        element of the header.
            e.g. #Names, A, B, C, D
                row1, 2, 4, 0, 4
                row2, 4, 6, 2, -1

    tree_file: str
        Path to the tree file in Newick format. The leaf node labels should
        be the same as as row names in the heatmap file. E.g. row1, row2.

    output_file: str, optional
        If specified the heat-tree will be rendered in that file as a PNG image,
        otherwise interactive browser will pop-up. **N.B.** program will wait
        for you to exit the browser before continuing.
    '''

    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    taxid2organism = manipulate_biosqldb.taxon_id2genome_description(server, biodb, True)

    t1 = Tree(tree_file)

    # Calculate the midpoint node
    R = t1.get_midpoint_outgroup()
    # and set it as tree outgroup
    t1.set_outgroup(R)

    leaf_number = 0
    for lf in t1.iter_leaves():
        leaf_number+=1
        lf.branch_vertical_margin = 0
        try:
            data = [taxid2n[str(lf.name)]]
        except:
            data=[0]
        #print 'taxon', int(lf.name)
        lf.name = taxid2organism[int(lf.name)]
        for col, value in enumerate(data):
            if value > 0:
                n = TextFace(' %s ' % str(value))
                n.margin_top = 2
                n.margin_right = 2
                n.margin_left = 2
                n.margin_bottom = 2
                n.inner_background.color = "#81BEF7"
                n.opacity = 1.
                lf.add_face(n, col, position="aligned")

            else:
                n = TextFace(' %s ' % str(value))
                n.margin_top = 2
                n.margin_right = 2
                n.margin_left = 2
                n.margin_bottom = 2
                n.inner_background.color = "white"
                n.opacity = 1.
                lf.add_face(n, col, position="aligned")

    return t1, leaf_number


def plot_heatmap_tree_locus(biodb,
                            tree_file,
                            taxid2count,
                            taxid2identity=False,
                            taxid2locus=False,
                            reference_taxon=False,
                            n_paralogs_barplot=False):

    '''

    plot tree and associated heatmap with count of homolgs
    optional:
        - add identity of closest homolog
        - add locus tag of closest homolog

    '''


    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    taxid2organism = manipulate_biosqldb.taxon_id2genome_description(server, biodb, True)

    t1 = Tree(tree_file)
    ts = TreeStyle()
    ts.draw_guiding_lines = True
    ts.guiding_lines_color = "gray"
    # Calculate the midpoint node
    R = t1.get_midpoint_outgroup()
    # and set it as tree outgroup
    t1.set_outgroup(R)

    leaf_number = 0

    for lf in t1.iter_leaves():

        if str(lf.name) not in taxid2count:
            taxid2count[str(lf.name)] = 0

    max_count = max([taxid2count[str(lf.name)] for lf in t1.iter_leaves()])

    for i, lf in enumerate(t1.iter_leaves()):
        
        # top leaf, add header
        if i == 0:
            
            n = TextFace('Number of homologs')
            n.margin_top = 1
            n.margin_right = 1
            n.margin_left = 20
            n.margin_bottom = 1
            n.inner_background.color = "white"
            n.opacity = 1.
            n.rotation = -25
            #lf.add_face(n, 7, position="aligned")
            ts.aligned_header.add_face(n, 1)
 
            if taxid2identity:
                n = TextFace('Protein identity')
                n.margin_top = 1
                n.margin_right = 1
                n.margin_left = 20
                n.margin_bottom = 1
                n.inner_background.color = "white"
                n.opacity = 1.
                n.rotation = -25
                #lf.add_face(n, 7, position="aligned")
                ts.aligned_header.add_face(n, 2)
            if taxid2locus:
                n = TextFace('Locus tag')
                n.margin_top = 1
                n.margin_right = 1
                n.margin_left = 20
                n.margin_bottom = 1
                n.inner_background.color = "white"
                n.opacity = 1.
                n.rotation = -25
                #lf.add_face(n, 7, position="aligned")
                ts.aligned_header.add_face(n, 3)
        
        leaf_number+=1

        lf.branch_vertical_margin = 0

        data = [taxid2count[str(lf.name)]]

        # possibility to add one or more columns
        for col, value in enumerate(data):
            col_index = col
            if value > 0:
                n = TextFace(' %s ' % str(value))
                n.margin_top = 2

                n.margin_right = 2
                if col == 0:
                    n.margin_left = 20
                else:
                    n.margin_left = 2
                n.margin_bottom = 2
                n.inner_background.color = "white" # #81BEF7
                n.opacity = 1.
                lf.add_face(n, col, position="aligned")

            else:
                n = TextFace(' %s ' % str(value))
                n.margin_top = 2
                n.margin_right = 2
                if col == 0:
                    n.margin_left = 20
                else:
                    n.margin_left = 2
                n.margin_bottom = 2
                n.inner_background.color = "white"
                n.opacity = 1.
                lf.add_face(n, col, position="aligned")
        # optionally indicate number of paralogs as a barplot
        if n_paralogs_barplot:
            col_index += 1
            percent = (float(value)/max_count)*100
            n = StackedBarFace([percent, 100-percent], width=150, height=18, colors=['#6699ff', 'white'], line_color='white')
            n.rotation= 0
            n.inner_border.color = "white"
            n.inner_border.width = 0
            n.margin_right = 15
            n.margin_left = 0
            lf.add_face(n, col+1, position="aligned")

        # optionally add additionnal column with identity
        if taxid2identity:
            import matplotlib.cm as cm
            from matplotlib.colors import rgb2hex
            import matplotlib as mpl

            norm = mpl.colors.Normalize(vmin=0, vmax=100)
            cmap = cm.OrRd
            m = cm.ScalarMappable(norm=norm, cmap=cmap)

            try:
                if round(taxid2identity[str(lf.name)], 2) != 100:
                    value = "%.2f" % round(taxid2identity[str(lf.name)], 2)
                else:
                    value = "%.1f" % round(taxid2identity[str(lf.name)], 2)
            except:
                value = '-'
            if str(lf.name) == str(reference_taxon):
                value = '         '
            n = TextFace(' %s ' % value)
            n.margin_top = 2
            n.margin_right = 2
            n.margin_left = 20
            n.margin_bottom = 2
            if not value.isspace() and value is not '-':
                n.inner_background.color = rgb2hex(m.to_rgba(float(value)))
                if float(value) > 82:
                    n.fgcolor = 'white'
            n.opacity = 1.
            if str(lf.name) == str(reference_taxon):
                n.inner_background.color = '#800000'

            lf.add_face(n, col_index+1, position="aligned")
        # optionaly add column with locus name
        if taxid2locus:
            try:
                value = str(taxid2locus[str(lf.name)])
            except:
                value = '-'
            n = TextFace(' %s ' % value)
            n.margin_top = 2
            n.margin_right = 2
            n.margin_left = 2
            n.margin_bottom = 2
            if str(lf.name) != str(reference_taxon):
                n.inner_background.color = "white"
            else:
                n.fgcolor = '#ff0000'
                n.inner_background.color = "white"
            n.opacity = 1.
            lf.add_face(n, col_index+2, position="aligned")
        lf.name = taxid2organism[str(lf.name)]

    return t1, leaf_number, ts


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    parser = argparse.ArgumentParser()
    parser.add_argument("-t",'--tree',type=str,help="newick tree")
    parser.add_argument("-m",'--matrix',type=str,help="matrix (tab file)")
    parser.add_argument("-s",'--mlst',type=str,help="mlst file")
    parser.add_argument("-a",'--accession2description',default=False, type=str,help="tab file with accessions and corresponding descriptions for leaf labels")
    parser.add_argument("-g",'--gbk_files', default=False, nargs='+', help="genbank files to get accession2description conversion")

    args = parser.parse_args()

    '''
    server, db = manipulate_biosqldb.load_db('chlamydia_03_15')
    biodb = 'chlamydia_03_15'
    orthogroup = 'group_15'

    sql_pfam = 'select taxon_id, count(*) from interpro_chlamydia_03_15 where analysis="Pfam" and signature_accession="PF00137" group by taxon_id;'

    sql_grp = 'select taxon_id,count(*) from  orthology_detail where orthogroup="%s" group by organism;' % (biodb, orthogroup)

    sql_cog = 'select t1.taxon_id, count(*) from (select * from COG.locus_tag2gi_hit_chlamydia_03_15 where COG_id="COG2911") A inner join biosqldb.orthology_detail_chlamydia_03_15 as t1 on A.locus_tag=t1.locus_tag group by t1.taxon_id;'

    taxid2n = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_cog,))

    plot_heat_tree('chlamydia_03_15', taxid2n, args.matrix, args.tree)


    '''

    if args.accession2description:
        accession2description = {}
        with open(args.accession2description, 'r') as f:
            for row in f:
                data = row.rstrip().split('\t')
                accession2description[data[0]] = data[1]
    elif args.gbk_files:
        import gbk2accessiontodefinition
        accession2description = gbk2accessiontodefinition.get_coressp(args.gbk_files)
    else:
        accession2description = False


    taxid2st = {}
    with open(args.mlst) as f:
        for n, row in enumerate(f):
            if n==0:
                continue
            else:
                data = row.rstrip().split('\t')
                taxid2st[data[0]] = data[2]
    taxid2n = {}
    with open(args.matrix) as f:
        for n, row in enumerate(f):
            if n==0:
                genes = row.rstrip().split('\t')[1:]
            else:
                data = row.rstrip().split('\t')
                taxid2n[data[0]] = [float(i) for i in data[1:]]

    t, n, style = plot_heat_tree_V1(taxid2n, args.tree, genes,taxid2st, accession2description)
    t.render("test.svg", tree_style=style)
