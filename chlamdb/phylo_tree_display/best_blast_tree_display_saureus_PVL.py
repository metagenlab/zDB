#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

import sys
from random import sample
from random import randint
from ete2 import Tree, SeqMotifFace, TreeStyle, add_face_to_node, TextFace, faces


import numpy as np
import colorsys
import matplotlib.colors as pcolors
import matplotlib.cm as cm
from matplotlib.colors import rgb2hex
import matplotlib as mpl

norm = mpl.colors.Normalize(vmin=70, vmax=100)
cmap = cm.OrRd
cmap_blue = cm.Blues
m = cm.ScalarMappable(norm=norm, cmap=cmap)
m2 = cm.ScalarMappable(norm=norm, cmap=cmap_blue)
def parse_mlst_results(mlst_file):
    '''
    Parse mlst result from script wrote by torsen seeman
    return accession2ST_type
    :return:
    '''

    accession2st_type = {}
    with open(mlst_file, 'r') as f:
        for line in f:
            line = line.rstrip().split('\t')
            accession2st_type[line[0]] = line[2]
    return accession2st_type





def plot_blast_result(tree_file, blast_result_file_list, id2description, id2mlst):
    '''
    Projet Staph aureus PVL avec Laure Jaton
    Script pour afficher une phylog�nie et la conservation de facteurs de firulence c�te � c�te
    N�cessite r�sultats MLST, ensemble des r�sultats tblastn (facteurs de virulence vs chromosomes),
    ainsi qu'une correspondance entre les accession des g�nomes et les noms qu'on veut afficher dans la phylog�nie. Icemn
    pour les identifiants molis des patients, on les remplace par CHUV n.
    :param tree_file: phylog�nie au format newick avec identifiants correspondants � tous les dico utilis�s
    :param blast_result_file_list: r�sultats tblastn virulence factors vs chromosome (seulement best blast)
    :param id2description: identifiants g�nome utiis� dans l'arbre et description correspondante (i.e S aureus Newman)
    :param id2mlst: identitifiants arbre 2 S. aureus ST type
    :return:
    '''

    blast2data = {}
    queries = []
    for one_blast_file in blast_result_file_list:
        with open(one_blast_file, 'r') as f:
            for line in f:
                line = line.split('\t')
                if line[1] not in blast2data:
                    blast2data[line[1]] = {}
                    blast2data[line[1]][line[0]] = [float(line[2]), int(line[8]), int(line[9])]
                else:
                     blast2data[line[1]][line[0]] = [float(line[2]), int(line[8]), int(line[9])]
                if line[0] not in queries:
                    queries.append(line[0])
    print blast2data
    print queries


    for one_blast in blast2data.keys():
        for ref_gene in blast2data[one_blast].keys():

            for query_gene in blast2data[one_blast].keys():
                overlap = False
                if ref_gene == query_gene:
                    continue
                if one_blast == 'NC_002745' and ref_gene == 'selm':
                    print 'target:', ref_gene, blast2data[one_blast][ref_gene]
                    print query_gene, blast2data[one_blast][query_gene]
                # check if position is overlapping
                try:
                    sorted_coordinates = sorted(blast2data[one_blast][ref_gene][1:3])
                    if blast2data[one_blast][query_gene][1] <= sorted_coordinates[1] and blast2data[one_blast][query_gene][1]>= sorted_coordinates[0]:
                        print 'Overlaping locations!'
                        print one_blast, ref_gene, blast2data[one_blast][ref_gene]
                        print one_blast, query_gene, blast2data[one_blast][query_gene]
                        overlap =True
                    sorted_coordinates = sorted(blast2data[one_blast][query_gene][1:3])
                    if blast2data[one_blast][ref_gene][1] <= sorted_coordinates[1] and blast2data[one_blast][ref_gene][1]>= sorted_coordinates[0]:
                        print 'Overlapping locations!'
                        print one_blast, ref_gene, blast2data[one_blast][ref_gene]
                        print one_blast, query_gene, blast2data[one_blast][query_gene]
                        overlap =True
                    if overlap:
                        if blast2data[one_blast][ref_gene][0] > blast2data[one_blast][query_gene][0]:
                            del blast2data[one_blast][query_gene]
                            print 'removing', query_gene
                        else:
                            del blast2data[one_blast][ref_gene]
                            print 'removing', ref_gene
                            break
                except KeyError:
                    print 'colocation already resolved:', query_gene, ref_gene


    '''
    rm_genes = ['selv','spsmA1','psmB1','psmB2','ses','set','sel','selX','sek','sel2','LukF', 'LukM', 'hly', 'hld'
        , 'hlgA', 'hlgB', 'hlgC', 'selq1', 'sec3', 'sek2', 'seq2', 'lukD', 'lukE', 'eta', 'etb', 'sec', 'tst']
    #rm_genes = ['icaR','icaA','icaB','icaC','icaD', 'sdrF', 'sdrH']

    for gene in rm_genes:
        queries.pop(queries.index(gene))
    '''
    #queries = ['selv']
    t1 = Tree(tree_file)
    #t.populate(8)
    # Calculate the midpoint node
    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    t1.ladderize()
    #t2=t1

    #for lf in t2.iter_leaves():
    #    try:
    #        lf.name = ' %s (%s)' % (id2description[lf.name], id2mlst[lf.name])
    #    except:
    #        #lf.name = ' %s (%s)' % (lf.name, lf.name)
    #
    #        a = TextFace(' %s (%s)' % (lf.name, id2mlst[lf.name]))
    #        a.fgcolor = "red"

    #        lf.name = a
    #t2.render("test.svg", dpi=800, h=400)
    #import sys
    #sys.exit()
        
    # and set it as tree outgroup
    for lf in t1.iter_leaves():
        #lf.add_face(AttrFace("name", fsize=20), 0, position="branch-right")
        lf.branch_vertical_margin = 0
        #data = [random.randint(0,2) for x in xrange(3)]

        for col, value in enumerate(queries):

            if lf.name == "1505183575":
                    'first row, print gene names'
                    #print 'ok!'
                    n = TextFace(' %s ' % str(value))
                    n.margin_top = 4
                    n.margin_right = 4
                    n.margin_left = 4
                    n.margin_bottom = 4
                    n.inner_background.color = "white"
                    n.opacity = 1.
                    lf.add_face(n, col, position="aligned")

            try:

                identity_value = blast2data[lf.name][value][0]

                if 'CHUV' in id2description[lf.name]:

                    if float(identity_value) >70:
                        if str(identity_value) == '100.00' or str(identity_value) == '100.0':
                            identity_value = '100'
                        else:
                            identity_value = str(round(float(identity_value), 1))
                        n = TextFace(' %s ' % str(identity_value))
                        n.margin_top = 4
                        n.margin_right = 4
                        n.margin_left = 4
                        n.margin_bottom = 4
                        n.inner_background.color = rgb2hex(m.to_rgba(float(identity_value)))
                        if float(identity_value) >92:
                            n.fgcolor = "white"
                        n.opacity = 1.
                        lf.add_face(n, col, position="aligned")
                    else:
                        identity_value = '-'
                        n = TextFace(' %s ' % str(identity_value))
                        n.margin_top = 2
                        n.margin_right = 2
                        n.margin_left = 2
                        n.margin_bottom = 2
                        n.inner_background.color = "white"
                        n.opacity = 1.
                        lf.add_face(n, col, position="aligned")
                else:
                    if float(identity_value) >70:
                        if str(identity_value) == '100.00' or str(identity_value) == '100.0':
                            identity_value = '100'
                        else:
                            identity_value = str(round(float(identity_value), 1))
                        n = TextFace(' %s ' % str(identity_value))
                        n.margin_top = 2
                        n.margin_right = 2
                        n.margin_left = 2
                        n.margin_bottom = 2
                        n.inner_background.color = rgb2hex(m2.to_rgba(float(identity_value)))

                        if float(identity_value) >92:
                            n.fgcolor = "white"

                        n.opacity = 1.
                        lf.add_face(n, col, position="aligned")
                    else:
                        identity_value = '-'
                        n = TextFace(' %s ' % str(identity_value))
                        n.margin_top = 2
                        n.margin_right = 2
                        n.margin_left = 2
                        n.margin_bottom = 2
                        n.inner_background.color = "white"
                        n.opacity = 1.
                        lf.add_face(n, col, position="aligned")
            except KeyError:

                identity_value = '-'
                n = TextFace(' %s ' % str(identity_value))
                n.margin_top = 2
                n.margin_right = 2
                n.margin_left = 2
                n.margin_bottom = 2
                n.inner_background.color = "white"
                n.opacity = 1.
                lf.add_face(n, col, position="aligned")

        try:
            lf.name = ' %s (%s)' % (id2description[lf.name], id2mlst[lf.name])
        except:
            #lf.name = ' %s (%s)' % (lf.name, lf.name)

            a = TextFace(' %s (%s)' % (lf.name, id2mlst[lf.name]))
            a.fgcolor = "red"

            lf.name = a
            #.add_face(a, 0, position="aligned")
    # add boostrap suppot
    #for n in t1.traverse():
    #    if n.is_leaf():
    #        continue
    #    n.add_face(TextFace(str(t1.support)), column=0, position = "branch-bottom")
    #ts = TreeStyle()
    #ts.show_branch_support = True

    # , tree_style=ts
    t1.render("saureus_tree.svg", dpi=800, h=400)
    t1.write(format=0, outfile="new_tree.nw")


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    import gbk2accessiontodefinition
    parser = argparse.ArgumentParser()
    parser.add_argument("-t",'--tree',type=str,help="newick tree")
    parser.add_argument("-b",'--blast_files',type=str,help="Blast tab files (no headers)", nargs="+")
    parser.add_argument("-g",'--input_gbk',type=str,help="gbff files", nargs="+")
    parser.add_argument("-m",'--molis_table',type=str,help="molis table")
    parser.add_argument("-s",'--st_file',type=str,help="mlst result file (t. seeman script)")

    args = parser.parse_args()
    id2description = gbk2accessiontodefinition.get_coressp(args.input_gbk, args.molis_table)
    accession2st_type = parse_mlst_results(args.st_file)
    plot_blast_result(args.tree, args.blast_files, id2description, accession2st_type)
