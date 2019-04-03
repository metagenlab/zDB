#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

def plot_phylo(nw_tree,
               out_name,
               parenthesis_classif=True,
               show_support=False,
               radial_mode=False,
               root=False):



    from ete3 import Tree, AttrFace,TreeStyle,NodeStyle, TextFace
    import orthogroup2phylogeny_best_refseq_uniprot_hity

    ete2_tree = Tree(nw_tree, format=0)
    if root:
        R = ete2_tree.get_midpoint_outgroup()
        # and set it as tree outgroup
        ete2_tree.set_outgroup(R)
    ete2_tree.set_outgroup('Bacillus subtilis')
    ete2_tree.ladderize()

    if parenthesis_classif:
        print ('parenthesis_classif!')
        name2classif = {}
        for lf in ete2_tree.iter_leaves():
            print (lf)
            try:
                classif = lf.name.split('_')[-2][0:-1]
                print ('classif', classif)
                #lf.name = lf.name.split('(')[0]
                name2classif[lf.name] = classif
            except:
                pass
        classif_list = list(set(name2classif.values()))
        classif2col = dict(zip(classif_list, orthogroup2phylogeny_best_refseq_uniprot_hity.get_spaced_colors(len(classif_list))))


    for lf in ete2_tree.iter_leaves():

        #try:
        if parenthesis_classif:
            try:
                col = classif2col[name2classif[lf.name]]
            except:
                col = 'black'
        else:
            col = 'black'
            #print col
            #lf.name = '%s|%s-%s' % (lf.name, accession2name_and_phylum[lf.name][0],accession2name_and_phylum[lf.name][1])

        if radial_mode:
            ff = AttrFace("name", fsize=12, fstyle = 'italic')
        else:
            ff = AttrFace("name", fsize=12, fstyle = 'italic')
        #ff.background.color = 'red'
        ff.fgcolor = col

        lf.add_face(ff, column=0)

        if not show_support:
            print('support')
            for n in ete2_tree.traverse():
               print (n.support)
               nstyle = NodeStyle()
               if float(n.support) < 1:
                   nstyle["fgcolor"] = "red"
                   nstyle["size"] = 4
                   n.set_style(nstyle)
               else:
                   nstyle["fgcolor"] = "red"
                   nstyle["size"] = 0
                   n.set_style(nstyle)
        else:
            for n in ete2_tree.traverse():
               nstyle = NodeStyle()
               nstyle["fgcolor"] = "red"
               nstyle["size"] = 0
               n.set_style(nstyle)
         
        #nameFace = AttrFace(lf.name, fsize=30, fgcolor=phylum2col[accession2name_and_phylum[lf.name][1]])
        #faces.add_face_to_node(nameFace, lf, 0, position="branch-right")
        #
        #nameFace.border.width = 1
        '''
        except:
            col = 'red'
            print col
            lf.name = '%s| %s' % (lf.name, locus2organism[lf.name])

            ff = AttrFace("name", fsize=12)
            #ff.background.color = 'red'
            ff.fgcolor = col

            lf.add_face(ff, column=0)
        '''
        #n = TextFace(lf.name, fgcolor = "black", fsize = 12, fstyle = 'italic')
        #lf.add_face(n, 0)
    '''
    for n in ete2_tree.traverse():
       nstyle = NodeStyle()
       if n.support < 90:
           nstyle["fgcolor"] = "black"
           nstyle["size"] = 4
           n.set_style(nstyle)
       else:
           nstyle["fgcolor"] = "red"
           nstyle["size"] = 0
           n.set_style(nstyle)
    '''
    ts = TreeStyle()
    ts.show_leaf_name = False
    #ts.scale=2000
    #ts.scale=20000
    ts.show_branch_support = show_support

    if radial_mode:
        ts.mode = "c"
        ts.arc_start = -90
        ts.arc_span = 360
    ts.tree_width = 370
    ts.complete_branch_lines_when_necessary = True
    ete2_tree.render(out_name, tree_style=ts, w=900)

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_tree', type=str, help="input newick tree")
    parser.add_argument("-s", '--support', action="store_true", help="show support")
    parser.add_argument("-rm", '--root_midpoint', action="store_true", help="Root the tree at midpoint")
    parser.add_argument("-r", '--radial', action="store_true", help="radial display")
    parser.add_argument("-c", '--color', action="store_true", help="color based on taxon: taxon2color table")
    parser.add_argument("-c2", '--color2', action="store_true", help="color based on parenthesis content")

    args = parser.parse_args()

    plot_phylo(args.input_tree, 'test.svg',
               parenthesis_classif=args.color2,
               show_support=args.support,
               radial_mode=args.radial,
               root=args.root_midpoint)
