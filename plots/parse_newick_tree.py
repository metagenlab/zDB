#!/usr/bin/env python


def get_accssion2ST(infile):
    f = open(infile, 'r')
    accession2st={}
    for line in f:
        row = line.strip().split('\t')
        st = row[0]
        for i in range(1, len(row)):
            accession2st[row[i].split('.')[0]]=st

    return accession2st

def convert_terminal_node_names(tree_newick_name, dictionnary, tree_format = 'phyloxml', write_corresp_table=True):
    #from Bio import Phylo
    from ete2 import Tree

    if write_corresp_table:
        f =  open('corresp_table.tab', 'w')

    print 'format', tree_format
    print 'tree', tree_newick_name
    try:
        tree = Tree(tree_newick_name, tree_format)
        leaf_names = {}
        for leaf in tree:
            #n+=1
            #print species.name
            #print dictionnary[str(species.name)]

            #accession = manipulate_biosqldb.taxon_id2accessions(server, str(species.name), "saureus_01_15")[0]
            try:
                print "%s --> %s" % (leaf.name, dictionnary[str(leaf.name)])
                if write_corresp_table:
                    f.write("%s\t%s\n" % (leaf.name, dictionnary[str(leaf.name)]))
                if dictionnary[str(leaf.name)] not in leaf_names:
                    leaf_names[dictionnary[str(leaf.name)]] = 1
                    leaf.name = dictionnary[str(leaf.name)] #+ ' (%s)' % my_accession2st[accession]
                else:
                    leaf_names[dictionnary[str(leaf.name)]] += 1
                    leaf.name = dictionnary[str(leaf.name)] + '%s' % leaf_names[dictionnary[str(leaf.name)]]#+ ' (%s)' % my_accession2st[accession]
            except:
                print "%s --> ?" % (leaf.name)
                if write_corresp_table:
                    f.write("%s\t?\n" % (leaf.name))
                pass
    except:
        from Bio import Phylo
        tree = Phylo.read(tree_newick_name, tree_format)

        leaves = tree.get_terminals()

        leaf_names = {}
        print dictionnary
        for leaf in leaves:

            try:
                print "%s --> %s" % (leaf.name, dictionnary[str(leaf.name)])
                if write_corresp_table:
                    f.write("%s\t%s\n" % (leaf.name, dictionnary[str(leaf.name)]))
                print "%s --> %s" % (leaf.name, dictionnary[str(leaf.name)])
                if dictionnary[str(leaf.name)] not in leaf_names:
                    leaf_names[dictionnary[str(leaf.name)]] = 1
                    leaf.name = dictionnary[str(leaf.name)] #+ ' (%s)' % my_accession2st[accession]
                else:
                    leaf_names[dictionnary[str(leaf.name)]] += 1
                    leaf.name = dictionnary[str(leaf.name)] + '%s' % leaf_names[dictionnary[str(leaf.name)]]#+ ' (%s)' % my_accession2st[accession]
            except:
                print "%s --> ?" % (leaf.name)
                if write_corresp_table:
                    f.write("%s\t?\n" % (leaf.name))
                print "%s --> ?" % (leaf.name)
                pass
            #species.name = dictionnary[str(species.name)] + ' (%s)' % 8
    return tree


color_map = {"C. abortus" : "blue",
            "C. avium" : "blue",
            "C. avium p." : "blue",
            "C. caviae" : "blue",
            "C. caviae p." : "blue",
            "C. felis" : "blue",
            "C. felis p." : "blue",
            "C. gallinacea" : "blue",
            "C. ibidis" : "blue",
            "C. muridarum" : "blue",
            "C. muridarum p." : "blue",
            "C. pecorum" : "blue",
            "C. pneumoniae" : "blue",
            "C. psittaci" : "blue",
            "C. psittaci p." : "blue",
            "Cr. sequanensis" : "red",
            "Cr. sequanensis p." : "red",
            "C. suis" : "blue",
            "C. trachomatis" : "blue",
            "E. lausannensis" : "red",
            "E. lausannensis p." : "red",
            "N. hartmannellae" : "green",
            "P. acanthamoebae" : "green",
            "Pr. amoebophila" : "green",
            "S. negevensis" : "orange",
            "S. negevensis p." : "orange",
            "W. chondrophila" : "magenta",
            "W. chondrophila p." : "magenta"
}





def combine_terminal_node_names(tree_newick_name, dictionnary):
    map_file = ""
    from Bio import Phylo
    import re
    trees = Phylo.parse(tree_newick_name, 'newick')
    trees = [tree for tree in trees]
    for tree in trees:
        #print tree
        #Phylo.draw_ascii(tree)
        tree.ladderize()   # Flip branches so deeper clades are displayed at top
        leaves = tree.get_terminals()
        for species in leaves:
            #print dictionnary[str(species.name)]
            species.name = re.sub("ElaC", "Estrella_lausannensis", species.name)
            species.name = re.sub("cibi", "cibid", species.name)
            species.name = re.sub("cibi_", "cibid_", species.name)
            species.name = re.sub("Cavi", "cavium", species.name)
            species.name = re.sub("Cavp", "caviump", species.name)
            if species.name == "Estrella_p":
                map_file += 'stroke:red Individual Estrella_p\n'
                pass

            elif species.name == "BAE81780.1":
               new_name = "CfeFp (%s)" "BAE81780.1"
               map_file += '\'fill:%s;stroke-width:3;stroke:%s\' Individual \'%s\'\n' % ("blue", "blue", new_name)
               species.name = new_name

            elif species.name == "BAE81781.1":
                new_name = "CfeFp (%s)" "BAE81781.1"
                map_file += '\'fill:%s;stroke-width:3;stroke:%s\' Individual \'%s\'\n' % ("blue", "blue", new_name)
                species.name = new_name
            elif species.name == "BAE81783.1":
                new_name = "CfeFp (%s)" "BAE81783.1"
                map_file += '\'fill:%s;stroke-width:3;stroke:%s\' Individual \'%s\'\n' % ("blue", "blue", new_name)
                species.name = new_name
            elif species.name == "BAE81782.1":
                new_name = "CfeFp (%s)" "BAE81782.1"
                map_file += '\'fill:%s;stroke-width:3;stroke:%s\' Individual "%s"\n' % ("blue", "blue", new_name)
                species.name = new_name
            elif species.name == "BAE81778.1":
                new_name = "CfeFp (%s)" "BAE81778.1"
                map_file += '\'fill:%s;stroke-width:3;stroke:%s\' Individual \'%s\'\n' % ("blue", "blue", new_name)
                species.name = new_name
            elif species.name == "BAE81779.1":
                new_name = "CfeFp (%s)" "BAE81779.1"
                map_file += '\'fill:%s;stroke-width:3;stroke:%s\' Individual \'%s\'\n' % ("blue", "blue", new_name)
                species.name = new_name
            elif species.name == "BAE81784.1":
                new_name = "CfeFp (%s)" "BAE81784.1"
                map_file += '\'fill:%s;stroke-width:3;stroke:%s\' Individual \'%s\'\n' % ("blue", "blue", new_name)
                species.name = new_name
            elif species.name == "BAE81785.1":
                new_name = "CfeFp (%s)" "BAE81785.1"
                map_file += '\'fill:%s;stroke-width:3;stroke:%s\' Individual \'%s\'\n' % ("blue", "blue",new_name)
                species.name = new_name


            else:
                new_name = dictionnary[str(species.name)] + "(%s)" % species.name
                map_file += '\'fill:%s;stroke-width:3;stroke:%s\' Individual \'%s\'\n' % (color_map[dictionnary[str(species.name)]], color_map[dictionnary[str(species.name)]], new_name)
                species.name = new_name
    return (trees, map_file)



