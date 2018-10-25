#!/usr/bin/env python



mat = '/home/tpillone/work/projets/dev/dev_report/typing/gatk_gvcfs/cgMLST/bwa/distances_in_snp.tsv'

import numpy
import networkx
import pandas
import json
import matplotlib.pyplot as plt

m = pandas.read_csv(mat, delimiter='\t', header=0, index_col=0)
A = m.values
print(m)
G = networkx.from_pandas_adjacency(m)#networkx.from_numpy_matrix(A)
print(G.edges(data=True))

T=networkx.minimum_spanning_tree(G)
print(sorted(T.edges(data=True)))

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