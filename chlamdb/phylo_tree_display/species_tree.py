

from chlamdb.biosqldb import manipulate_biosqldb

def species2taxon_id(biodb):

    server, db = manipulate_biosqldb.load_db(biodb)

    sql2 = 'select distinct family,taxon_id from taxid2species_%s t1 ' \
          ' inner join species_curated_taxonomy_%s t2 on t1.species_id=t2.species_id;' % (biodb, 
                                                                                          biodb)
    
    species_id2taxon_id = {}
    data = server.adaptor.execute_and_fetchall(sql2,)
    for row in data:
        if row[0] not in species_id2taxon_id:
            species_id2taxon_id[row[0]] = [str(row[1])]
        else:
            species_id2taxon_id[row[0]].append(str(row[1]))
    
    return species_id2taxon_id

def get_species_tree(biodb):
    
    from ete3 import Tree,TreeStyle
    
    server, db = manipulate_biosqldb.load_db(biodb)
    
    
    sql_tree = 'select tree from reference_phylogeny t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
               ' where t2.name="%s";' % biodb
               
    server, db = manipulate_biosqldb.load_db(biodb)
    complete_tree = Tree(server.adaptor.execute_and_fetchall(sql_tree,)[0][0])
    R = complete_tree.get_midpoint_outgroup()
    complete_tree.set_outgroup(R)

    sql = 'select distinct taxon_id,species from taxid2species_%s t1 ' \
          ' inner join species_curated_taxonomy_%s t2 on t1.species_id=t2.species_id;' % (biodb, 
                                                                                          biodb)
          
    taxon_id2species_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    
    
    # changing taxon id to species id
    for leaf in complete_tree.iter_leaves():
        #print '%s --> %s' % (leaf.name, str(taxon_id2species_id[str(leaf.name)]))
        leaf.name = "%s" % str(taxon_id2species_id[str(leaf.name)])

    # attributing unique id to each node
    # if all node descendant have the same name, use that name as node name
    n = 0
    for node in complete_tree.traverse():
        if node.name=='':
            desc_list = list(set([i.name for i in node.iter_descendants()]))
            try:
                desc_list.remove('')
            except ValueError:
                pass
            if len(desc_list) != 1:
                node.name = '%sbb' % n
            else:
                node.name = desc_list[0]
            n+=1
 
    # Collapsing nodes while traversing
    # http://etetoolkit.org/docs/latest/tutorial/tutorial_trees.html#collapsing-nodes-while-traversing-custom-is-leaf-definition
    node2labels = complete_tree.get_cached_content(store_attr="name")
    
    def collapsed_leaf(node):
        if len(node2labels[node]) == 1:
            return True
        else:
            return False

    species_tree = Tree(complete_tree.write(is_leaf_fn=collapsed_leaf))
    
    return complete_tree, species_tree