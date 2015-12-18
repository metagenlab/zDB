#!/usr/bin/env python

import sys
from random import sample
from random import randint
from ete2 import Tree, SeqMotifFace, TreeStyle, add_face_to_node, TextFace


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

        print rgb
        colors.append(pcolors.rgb2hex(rgb))
    return colors




def organism2color(locus2data):

    organism_list = []
    for locus in locus2data:
        if locus2data[locus][0][2] not in organism_list:
            organism_list.append(locus2data[locus][0][2])
    colors = _get_colors(len(organism_list))

    return dict(zip(organism_list, colors))


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
    import numpy

    from ete2.treeview.faces import add_face_to_node
    from ete2 import ClusterTree, TreeStyle, AttrFace, ProfileFace



    import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)


    taxid2organism = manipulate_biosqldb.taxon_id2genome_description(server, biodb, True)

    print taxid2organism

    t1 = Tree(tree_file)
    #t.populate(8)
    # Calculate the midpoint node
    R = t1.get_midpoint_outgroup()
    # and set it as tree outgroup
    t1.set_outgroup(R)    # To operate with numbers efficiently

    # Loads tree and array


    '''
    # Calculates some stats on the matrix. Needed to establish the color
    # gradients.
    matrix_dist = [i for r in xrange(len(array.matrix))\
                   for i in array.matrix[r] if numpy.isfinite(i)]
    matrix_max = numpy.max(matrix_dist)
    matrix_min = numpy.min(matrix_dist)
    matrix_avg = matrix_min+((matrix_max-matrix_min)/2)

    # Creates a profile face that will represent node's profile as a
    # heatmap
    profileFace  = ProfileFace(1., 0., 0.5, 100, 14, "heatmap",colorscheme=2)

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
    t.show(tree_style=ts)
    '''
    import random

    '''
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
    leaf_number = 0
    for lf in t1.iter_leaves():
        leaf_number+=1
        print lf
        #lf.add_face(AttrFace("name", fsize=20), 0, position="branch-right")
        lf.branch_vertical_margin = 0
        #data = [random.randint(0,2) for x in xrange(3)]
        try:
            data = [taxid2n[str(lf.name)]]
        except:
            data=[0]
        print 'taxon', int(lf.name)
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


    #a = t1.render(output_file, dpi=800, h=i*15)
    return t1, leaf_number


if __name__ == '__main__':
    import argparse
    import manipulate_biosqldb
    parser = argparse.ArgumentParser()
    parser.add_argument("-t",'--tree',type=str,help="newick tree")
    parser.add_argument("-m",'--matrix',type=str,help="matrix (tab file)")



    args = parser.parse_args()
    server, db = manipulate_biosqldb.load_db('chlamydia_03_15')
    biodb = 'chlamydia_03_15'
    orthogroup = 'group_15'

    sql_pfam = 'select taxon_id, count(*) from interpro_chlamydia_03_15 where analysis="Pfam" and signature_accession="PF00137" group by taxon_id;'

    sql_grp = 'select taxon_id,count(*) from  orthology_detail_%s where orthogroup="%s" group by organism;' % (biodb, orthogroup)

    sql_cog = 'select t1.taxon_id, count(*) from (select * from COG.locus_tag2gi_hit_chlamydia_03_15 where COG_id="COG2911") A inner join biosqldb.orthology_detail_chlamydia_03_15 as t1 on A.locus_tag=t1.locus_tag group by t1.taxon_id;'

    taxid2n = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_cog,))

    plot_heat_tree('chlamydia_03_15', taxid2n, args.matrix, args.tree)

