#!/usr/bin/env python
# python 2.7.5 requires biopython




def contruct_graph():
    '''
    import manipulate_biosqldb
    import igraph

    server, db = manipulate_biosqldb.load_db('chlamydia_04_16')

    sql = 'select * from comparative_tables.phylogenetic_profiles_euclidian_distance_chlamydia_04_16 where group_1!=group_2 and euclidian_dist=0'

    data = server.adaptor.execute_and_fetchall(sql,)

    all_verticles = []
    for i in data:
        all_verticles.append(i[0])
        all_verticles.append(i[1])
    all_verticles = list(set(all_verticles))

    print len(all_verticles)

    g = igraph.Graph()

    g.add_vertices(len(all_verticles))

    g.vs["label"] = all_verticles

    all_indexes = []
    for i in data:
        index_1 = all_verticles.index(i[0])
        index_2 = all_verticles.index(i[1])
        all_indexes.append((index_1, index_2))
    g.add_edges(all_indexes)
    igraph.summary(g)
    communities = g.community_edge_betweenness(directed=False)
    clusters = communities.as_clustering()
    print len(clusters)
    print clusters[0]

    layout = dg[0].layout("kk")
    igraph.plot(dg[0], 'test.pdf' ,layout = layout)
    '''

    import networkx as nx
    import manipulate_biosqldb
    import igraph
    import matplotlib.pyplot as plt

    server, db = manipulate_biosqldb.load_db('chlamydia_04_16')

    sql = 'select * from comparative_tables.phylogenetic_profiles_euclidian_distance_chlamydia_04_16 where group_1!=group_2 and euclidian_dist<=2'

    data = server.adaptor.execute_and_fetchall(sql,)



    all_verticles = []
    for i in data:
        all_verticles.append(i[0])
        all_verticles.append(i[1])
    all_verticles = list(set(all_verticles))

    #print len(all_verticles)

    G = nx.Graph()
    G.add_nodes_from(all_verticles)

    label = {}
    for n, orthogroup in enumerate(all_verticles):
        label[n] = orthology_chlamydia_04_16orthogroup
    print label[10], label[2388]

    for i in data:
        index_1 = i[0]#all_verticles.index(i[0])
        index_2 = i[1]#all_verticles.index(i[1])
        weight = i[2]
        G.add_edge(index_1, index_2, weight=weight)
        #all_indexes.append((index_1, index_2))
    #G.add_edges_from(all_indexes)

    graph_list = [i for i in nx.connected_component_subgraphs(G)]

    print len(graph_list)
    '''
    for n, graph in enumerate(graph_list):
        print n, len(graph)
        nx.draw(graph_list[n],node_color='red')
        plt.savefig("path_%s.png" % n)
    '''

    for i in range(0, 40):
        print 'fig %s' % i
        print len(graph_list[i])
        gg = graph_list[i]
        #remove = [node for node,degree in gg.degree().items() if degree < 2]
        #gg.remove_nodes_from(remove)
        fig = plt.figure()
        nx.draw(gg, node_color='red', label=label, with_labels = True)
        fig.savefig("path_%s.png" % i)
        fig.clear()




if __name__ == '__main__':
    ###Argument handling.
    import argparse
    arg_parser = argparse.ArgumentParser(description='');
    #arg_parser.add_argument("coords_input", help="Directory to show-coords tab-delimited input file.");
    arg_parser.add_argument("-g", "--genbank", help="genbank")

    args = arg_parser.parse_args()

    contruct_graph()