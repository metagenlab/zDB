#!/usr/bin/env python

import itertools


mat = '/home/tpillone/work/projets/2018_10_MYCOST/2018_07_4_strains/typing/freebayes_joint_genotyping/cgMLST/bwa/distances_in_snp.tsv'
# /home/tpillone/work/projets/2018_10_MYCOST/2018_07_4_strains/typing/freebayes_joint_genotyping/cgMLST/bwa
import numpy
import networkx
import pandas
import json
import matplotlib.pyplot as plt

m = pandas.read_csv(mat, delimiter='\t', header=0, index_col=0)


nodes = list(m.index)


missing = []


A = m.values

G = networkx.from_pandas_adjacency(m)#networkx.from_numpy_matrix(A)


def find_clusters(G, node_list):
    comb = itertools.combinations(node_list, 2)
    groups = []
    for i in comb:
        if i[0] == i[1]:
            continue
        try:
            print (G[i[0]][i[1]])
        except:
            if len(groups)==0:
                no_match=True
            else:
                no_match = False
            for n, group in enumerate(groups):
                if i[0] in group and i[1] in group:
                    continue
                elif i[0] in group and i[1] not in group:
                    groups[n].append(i[1])
                elif i[1] in group and i[0] not in group:
                    groups[n].append(i[1])
                else:
                    no_match=True
            if no_match:
                groups.append([i[0], i[1]])
    return groups

def merge_group_nodes(G, groups, node_list):
    for group in groups:
        median_dico = {}
        for node in node_list:
            if node in group:
                continue
            data = []
            for member in group:
                data.append(G[member][node]['weight'])
            m = numpy.median(data)
            median_dico[node] = m
        mapping = {group[0]: '\n'.join(group)}
        G = networkx.relabel_nodes(G, mapping)
        for i in group[1:len(group)]:
            G.remove_node(i)
        for i in node_list:
            if i in group:
                continue
            G['\n'.join(group)][i]['weight'] = median_dico[i]
    return G

groups = find_clusters(G, nodes)
print("groups", groups)
G = merge_group_nodes(G, groups, nodes)

T=networkx.minimum_spanning_tree(G)
#print(sorted(T.edges(data=True)))

#networkx.draw(T)
#plt.show()$

# this function is used to convert networkx to Cytoscape.js JSON format
# returns string of JSON
def convert2cytoscapeJSON(G):
    # load all nodes into nodes array
    final = {}
    final["nodes"] = []
    final["edges"] = []
    for node in G.nodes():
        nx = {}
        nx["data"] = {}
        nx["data"]["id"] = node
        nx["data"]["label"] = node
        final["nodes"].append(nx.copy())
    #load all edges to edges array
    for edge in G.edges(data=True):
        print(edge)
        nx = {}
        nx["data"]={}
        nx["data"]["id"]=edge[0]+edge[1]
        nx["data"]["strength"] = edge[2]["weight"]
        nx["data"]["source"]=edge[0]
        nx["data"]["target"]=edge[1]
        final["edges"].append(nx)
    return json.dumps(final)

print(convert2cytoscapeJSON(T))