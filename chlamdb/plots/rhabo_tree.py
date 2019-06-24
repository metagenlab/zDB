#!/usr/bin/python
#-*- coding: utf-8 -*-
#import heatmap
import numpy as np
from chlamdb.biosqldb import manipulate_biosqldb
import parse_newick_tree
from Bio import Phylo

# heatmap Chlamydiales pan-genome






def detailed_leaf_name(tree_name, id2quantity, id2Ct, id2location, id2canton, id2best_hit, id2best_hit_id, id2Ct2):
    from ete2 import Tree, SeqMotifFace, TreeStyle, add_face_to_node, TextFace
    t = Tree(tree_name)

        #t.populate(8)
    # Calculate the midpoint node
    #R = t.get_midpoint_outgroup()
    # and set it as tree outgroup
    #t.set_outgroup(R)
    i = 0
    for l in t.iter_leaves():
        i+=1
        #l.img_style['fgcolor'] = "red"
        #l.img_style['hz_line_type'] = 0
        #l.img_style
        # ['size'] = 10

        pcr_id =  str(l.name)
        #print 'pcr %s' % pcr_id
        #print id2quantity[pcr_id]
        try:
            if (float(id2quantity[pcr_id]) > 5000):
                print "%s\t%s" % (pcr_id, id2quantity[pcr_id])
                location = id2location[pcr_id]
                canton = id2canton[pcr_id]
                quantity = id2quantity[pcr_id]
                best_hit = id2best_hit[pcr_id]
                best_hit_id = id2best_hit_id[pcr_id]

                l.name = "%s %s (%s) %s %s %s" % (pcr_id, best_hit, best_hit_id, quantity, canton, location)
            else:
                l.name = i
        except:
            print '######### problem with: %s' % pcr_id
            l.name = i
    t.write(format=1, outfile="new_tree3.nw")

def parse_tique_tables(table_result, table_pos):
    with open(table_result, 'r') as f:
        id2quantity = {}
        id2Ct = {}
        for line in f:
            data = line.rstrip().split("\t")
            id2quantity[data[1].strip()] = data[8]
            id2Ct[data[1].strip()] = data[5]

    with open(table_pos, 'r') as g:
        id2location = {}
        id2canton = {}
        id2best_hit = {}
        id2best_hit_id = {}
        id2Ct2 = {}
        for line in g:
            data = line.rstrip().split("\t")
            id2location[data[0].strip()] = data[3]
            id2canton[data[0].strip()] = data[2]
            try:
                id2best_hit[data[0].strip()] = data[5]
                id2best_hit_id[data[0].strip()] = data[6]
            except:
                id2best_hit[data[0].strip()] = '-'
                id2best_hit_id[data[0].strip()] = '-'
            id2Ct2[data[0].strip()] = data[4]
    print len(id2location), id2location.keys()
    return (id2quantity, id2Ct, id2location, id2canton, id2best_hit, id2best_hit_id, id2Ct2)

########### MAIN ####################
#####################################

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--tree', type=str, help="newick tree")
    parser.add_argument("-p", '--pcr_result', type=str, help="PCR result table")
    parser.add_argument("-r", '--pcr_pos', type=str, help="PCR pos table")

    args = parser.parse_args()

    id2quantity, id2Ct, id2location, id2canton, id2best_hit, id2best_hit_id, id2Ct2 = parse_tique_tables(args.pcr_result, args.pcr_pos)
    #y§parse_tique_tables(args.pcr_result, args.pcr_pos, id2quantity, id2Ct, id2location, id2canton, id2best_hit, id2best_hit_id, id2Ct2)
    detailed_leaf_name(args.tree, id2quantity, id2Ct, id2location, id2canton, id2best_hit, id2best_hit_id, id2Ct2)