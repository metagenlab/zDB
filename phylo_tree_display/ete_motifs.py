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

        colors.append(pcolors.rgb2hex(rgb))
    return colors



def get_pfam_data(orthogroup, biodb):
    import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select A.protein_id, B.start, B.stop, A.organism, A.sequence_length, B.signature_accession, B.signature_description, A.taxon_id ' \
          ' from (select taxon_id, orthogroup, protein_id, organism, length(translation) as sequence_length from orthology_detail_%s ' \
          ' where orthogroup="%s" ) A ' \
          ' left join (select * from interpro_%s where orthogroup="%s" and analysis="Pfam") B ' \
          ' on A.protein_id=B.accession;' % (biodb, orthogroup, biodb, orthogroup)
    '''
    sql = 'select accession, start, stop, organism, sequence_length, signature_accession, signature_description  ' \
          ' from interpro_%s as t2 where orthogroup="%s" and analysis="Pfam"' % (biodb, orthogroup)
    '''
    data = server.adaptor.execute_and_fetchall(sql,)

    protein2data = {}
    for one_locus in data:
        if one_locus[0] not in protein2data:
            protein2data[one_locus[0]] = [one_locus[1:len(one_locus)]]
        else:
            protein2data[one_locus[0]].append(one_locus[1:len(one_locus)])
    return protein2data

def get_TM_data(biodb, orthogroup=False):
    import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    if orthogroup:
        sql = 'select locus_tag, start, stop, organism, sequence_length, signature_accession, signature_description  ' \
          ' from interpro_%s as t2 where orthogroup="%s" and analysis="Phobius" and signature_accession="TRANSMEMBRANE"' % (biodb, orthogroup)
    else:
        sql = 'select locus_tag, start, stop, organism, sequence_length, signature_accession, signature_description  ' \
          ' from interpro_%s as t2 where analysis="Phobius" and signature_accession="TRANSMEMBRANE"' % (biodb)

    data = server.adaptor.execute_and_fetchall(sql,)

    protein2data = {}
    for one_locus in data:
        if one_locus[0] not in protein2data:
            protein2data[one_locus[0]] = [one_locus[1:len(one_locus)]]
        else:
            protein2data[one_locus[0]].append(one_locus[1:len(one_locus)])
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
            if locus2data[locus][0][2] not in organism_list:
                organism_list.append(locus2data[locus][0][2])
        colors = _get_colors(len(organism_list))

        return dict(zip(organism_list, colors))
    else:
        family_list = []
        for locus in locus2data:
            taxon_id = locus2data[locus][0][-1]
            if taxon_id2family[str(taxon_id)] not in family_list:
                family_list.append(taxon_id2family[str(taxon_id)])
        colors = _get_colors(len(family_list))

        return dict(zip(family_list, colors))


def draw_pfam_tree(tree_name, locus2data, locus2protein_id = False, taxon_id2family=False):
    # Create a random tree and add to each leaf a random set of motifs
    # from the original set
    t = Tree(tree_name)
    #t.populate(8)
    # Calculate the midpoint node
    R = t.get_midpoint_outgroup()
    # and set it as tree outgroup
    t.set_outgroup(R)

    color_dico = organism2color(locus2data, taxon_id2family)

    leaf_number = 0
    for l in t.iter_leaves():
        leaf_number+=1
        if locus2protein_id:
            protein_id = locus2protein_id[str(l)[3:len(str(l))]]
            data = locus2data[protein_id]
        else:
            data = locus2data[str(l)[3:len(str(l))]]

        seq_motifs = []
        if not taxon_id2family:
            l.img_style['fgcolor'] = color_dico[data[0][2]]
        else:
            taxon_id = data[0][-1]
            family = taxon_id2family[str(taxon_id)]
            l.img_style['fgcolor'] = color_dico[family]
        l.img_style['hz_line_type'] = 0
        l.img_style['size'] = 10

        for motif in data:
            if motif[0]:

                seq_motifs.append([motif[0], motif[1], "[]", None, 10, "black", "PaleGreen", "arial|8|red|%s" % motif[4]])

        seqFace = SeqMotifFace(data[0][3]*'N',motifs=seq_motifs, width=10,height=12, intermotif_format='-', seqtail_format="-", seq_format='-') #seq, seq_motifs, scale_factor=1, intermotif_format=None) #intermotif_format="seq", seqtail_format=None, seqtype='aa')
        seqFace.margin_bottom = 2
        seqFace.margin_top = 2
        seqFace.opacity = 1.0

        l.add_face(seqFace, column=1, position="aligned")
        locus = TextFace(str(l)[3:len(str(l))])
        l.name = data[0][2]
        locus.margin_right = 10
        locus.margin_left = 10
        locus.margin_bottom = 0
        l.add_face(locus, column=0, position="aligned")



    ts = TreeStyle()
    #ts.layout_fn = layout
    return t, ts, leaf_number

def draw_TM_tree(tree_name, locus2data):
    # Create a random tree and add to each leaf a random set of motifs
    # from the original set
    t = Tree(tree_name)
    #t.populate(8)
    # Calculate the midpoint node
    R = t.get_midpoint_outgroup()
    # and set it as tree outgroup
    t.set_outgroup(R)


    color_dico = organism2color(locus2data)

    for l in t.iter_leaves():
        locus_name = str(l)[3:len(str(l))]
        locus_name = locus_name.split('|')[0]
        data = locus2data[locus_name]

        seq_motifs = []
        l.img_style['fgcolor'] = color_dico[data[0][2]]
        l.img_style['hz_line_type'] = 0
        l.img_style['size'] = 10


        for motif in data:

            seq_motifs.append([motif[0], motif[1], "[]", None, 10, "black", "PaleGreen", "arial|8|red|TM"])


        seqFace = SeqMotifFace(data[0][3]*'N',motifs=seq_motifs, width=10,height=12, intermotif_format='-', seqtail_format="-", seq_format='-') #seq, seq_motifs, scale_factor=1, intermotif_format=None) #intermotif_format="seq", seqtail_format=None, seqtype='aa')
        seqFace.margin_bottom = 2
        seqFace.margin_top = 2
        seqFace.opacity = 1.0

        l.add_face(seqFace, column=1, position="aligned")
        locus = TextFace(str(l)[3:len(str(l))])
        l.name = data[0][2]
        print dir(l), '########################'
        locus.margin_right = 10
        locus.margin_left = 10
        locus.margin_bottom = 0
        l.add_face(locus, column=0, position="aligned")

    ts = TreeStyle()
    #ts.layout_fn = layout
    return t, ts



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t",'--tree',type=str,help="newick tree")



    args = parser.parse_args()


    #locus2pfam_data = get_pfam_data("group_20", 'chlamydia_03_15')

    #t, ts = draw_pfam_tree(args.tree, locus2pfam_data)

    locus2TM_data = get_TM_data('chlamydia_03_15')
    t, ts = draw_TM_tree(args.tree, locus2TM_data)
    t.render("motifs.svg", w=1200, dpi=800, tree_style=ts)
    #t.show(tree_style=ts)


