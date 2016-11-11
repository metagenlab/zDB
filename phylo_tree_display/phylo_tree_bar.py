#!/usr/bin/env python

from ete2 import Tree, SeqMotifFace, TreeStyle, add_face_to_node, TextFace, BarChartFace, StackedBarFace, NodeStyle

def plot_heat_tree(tree_file, biodb="chlamydia_04_16"):
    import manipulate_biosqldb
    import matplotlib.cm as cm
    from matplotlib.colors import rgb2hex
    import matplotlib as mpl




    server, db = manipulate_biosqldb.load_db(biodb)

    sql_biodatabase_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb
    db_id = server.adaptor.execute_and_fetchall(sql_biodatabase_id,)[0][0]

    t1 = Tree(tree_file)

    sql1 = 'select taxon_id, description from bioentry where biodatabase_id=%s and description not like "%%%%plasmid%%%%"' % db_id
    sql2 = 'select t2.taxon_id, t1.GC from genomes_info_%s as t1 inner join bioentry as t2 ' \
           ' on t1.accession=t2.accession where t2.biodatabase_id=%s and t1.description not like "%%%%plasmid%%%%";' % (biodb, db_id)
    sql3 = 'select t2.taxon_id, t1.genome_size from genomes_info_%s as t1 ' \
           ' inner join bioentry as t2 on t1.accession=t2.accession ' \
           ' where t2.biodatabase_id=%s and t1.description not like "%%%%plasmid%%%%";' % (biodb, db_id)
    sql4 = 'select taxon_id,sum(stop-start) as coding from orthology_detail_%s group by taxon_id;' % biodb

    taxon2description = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql1,))
    taxon2gc = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
    taxon2genome_size = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3,))
    taxon2coding_size = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql4,))

    taxon2coding_density = {}
    for taxon in taxon2coding_size:
        taxon2coding_density[taxon] = float(taxon2coding_size[taxon])/float(taxon2genome_size[taxon])

    print sql2
    print taxon2genome_size

    my_taxons = [lf.name for lf in t1.iter_leaves()]

    genome_sizes = [float(taxon2genome_size[i]) for i in my_taxons]
    gc_list = [float(taxon2gc[i]) for i in my_taxons]
    fraction_list = [float(taxon2coding_density[i]) for i in my_taxons]
    # Calculate the midpoint node
    R = t1.get_midpoint_outgroup()
    # and set it as tree outgroup
    t1.set_outgroup(R)

    value=1

    max_genome_size = max(genome_sizes)#3424182#
    max_gc = max(gc_list) #48.23

    cmap = cm.YlGnBu#YlOrRd#OrRd
    print 'maxmin', max(genome_sizes)
    norm = mpl.colors.Normalize(vmin=min(genome_sizes)-100000, vmax=max(genome_sizes))
    m1 = cm.ScalarMappable(norm=norm, cmap=cmap)
    norm = mpl.colors.Normalize(vmin=min(gc_list), vmax=max(gc_list))
    m2 = cm.ScalarMappable(norm=norm, cmap=cmap)
    norm = mpl.colors.Normalize(vmin=min(fraction_list)*100, vmax=max(fraction_list)*100)
    m3 = cm.ScalarMappable(norm=norm, cmap=cmap)



    #print "max_genome_size", max_genome_size
    for i, lf in enumerate(t1.iter_leaves()):
        #if taxon2description[lf.name] == 'Pirellula staleyi DSM 6068':
        #    lf.name = 'Pirellula staleyi DSM 6068'
        #    continue
        if i==0:
            n = TextFace('Size (Mbp)')
            n.margin_top = 1
            n.margin_right = 1
            n.margin_left = 20
            n.margin_bottom = 1
            n.inner_background.color = "white"
            n.opacity = 1.
            lf.add_face(n, 3, position="aligned")
            n = TextFace('GC (%)')
            n.margin_top = 1
            n.margin_right = 1
            n.margin_left = 20
            n.margin_bottom = 1
            n.inner_background.color = "white"
            n.opacity = 1.
            lf.add_face(n, 5, position="aligned")
            n = TextFace('')
            lf.add_face(n, 2, position="aligned")
            lf.add_face(n, 4, position="aligned")
            n = TextFace('Coding density (%)')
            n.margin_top = 1
            n.margin_right = 1
            n.margin_left = 20
            n.margin_bottom = 1
            n.inner_background.color = "white"
            n.opacity = 1.
            lf.add_face(n, 7, position="aligned")
            n = TextFace('')
            lf.add_face(n, 6, position="aligned")

        value+=1
        n = TextFace('  %s ' % str(round(taxon2genome_size[lf.name]/float(1000000),2)))
        n.margin_top = 1
        n.margin_right = 1
        n.margin_left = 20
        n.margin_bottom = 1
        n.inner_background.color = "white"
        n.opacity = 1.

        lf.add_face(n, 2, position="aligned")
        print taxon2genome_size[lf.name]
        fraction_biggest = (float(taxon2genome_size[lf.name])/max_genome_size)*100
        fraction_rest = 100-fraction_biggest
        print fraction_biggest, fraction_rest
        if taxon2description[lf.name] == 'Rhabdochlamydia helveticae T3358':
            col = 'black'
        else:
            col = rgb2hex(m1.to_rgba(float(taxon2genome_size[lf.name])))  # 'grey'

        b = StackedBarFace([fraction_biggest, fraction_rest], width=100, height=10,colors=[col, 'white'])
        b.rotation= 0
        b.inner_border.color = "grey"
        b.inner_border.width = 0
        b.margin_right = 15
        b.margin_left = 0
        lf.add_face(b, 3, position="aligned")

        fraction_biggest = (float(taxon2gc[lf.name])/max_gc)*100
        fraction_rest = 100-fraction_biggest
        if taxon2description[lf.name] == 'Rhabdochlamydia helveticae T3358':
            col = 'black'
        else:
            col = rgb2hex(m2.to_rgba(float(taxon2gc[lf.name])))
        b = StackedBarFace([fraction_biggest, fraction_rest], width=100, height=10,colors=[col, 'white'])
        b.rotation= 0
        b.inner_border.color = "grey"
        b.inner_border.width = 0
        b.margin_left = 0
        b.margin_right = 15


        lf.add_face(b, 5, position="aligned")
        n = TextFace('  %s ' % str(round(float(taxon2gc[lf.name]),2)))
        n.margin_top = 1
        n.margin_right = 0
        n.margin_left = 0
        n.margin_bottom = 1
        n.inner_background.color = "white"
        n.opacity = 1.
        lf.add_face(n, 4, position="aligned")


        if taxon2description[lf.name] == 'Rhabdochlamydia helveticae T3358':
            col = 'black'
        else:
            col = rgb2hex(m3.to_rgba(float(taxon2coding_density[lf.name])*100))
        n = TextFace('  %s ' % str(round(float(taxon2coding_density[lf.name])*100, 1)))
        n.margin_top = 1
        n.margin_right = 0
        n.margin_left = 0
        n.margin_right = 0
        n.margin_bottom = 1
        n.inner_background.color = "white"
        n.opacity = 1.
        lf.add_face(n, 6, position="aligned")
        fraction = round(float(taxon2coding_density[lf.name]), 2)*100
        fraction_rest = 100 - fraction


        b = StackedBarFace([fraction, fraction_rest], width=100, height=10,colors=[col, 'white'])# 1-round(float(taxon2coding_density[lf.name]), 2)
        b.rotation = 0
        b.margin_right = 1
        b.inner_border.color = "grey"
        b.inner_border.width = 0
        b.margin_left = 5
        lf.add_face(b, 7, position="aligned")

        lf.name = taxon2description[lf.name]

    for n in t1.traverse():
       nstyle = NodeStyle()
       if n.support < 1:
           nstyle["fgcolor"] = "red"
           nstyle["size"] = 6
           n.set_style(nstyle)
       else:
           nstyle["fgcolor"] = "red"
           nstyle["size"] = 0
           n.set_style(nstyle)

    return t1




if __name__ == '__main__':
    import argparse
    import manipulate_biosqldb
    parser = argparse.ArgumentParser()
    parser.add_argument("-t",'--tree',type=str,help="newick tree")
    parser.add_argument("-m",'--matrix',type=str,help="matrix (tab file)")
    parser.add_argument("-s",'--mlst',type=str,help="mlst file")
    parser.add_argument("-d",'--biodb',type=str,help="biodatabase name")
    parser.add_argument("-a",'--accession2description',default=False, type=str,help="tab file with accessions and corresponding descriptions for leaf labels")

    args = parser.parse_args()

    if not args.tree:
        import manipulate_biosqldb
        sql_tree = 'select tree from reference_phylogeny t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                   ' where t2.name="%s";' % args.biodb
        server, db = manipulate_biosqldb.load_db(args.biodb)
        tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]
        print tree
    else:
        tree = args.tree

    t = plot_heat_tree(tree, args.biodb)
    ts = TreeStyle()
    ts.show_branch_support = False
    t.render("test2.svg", tree_style=ts)
