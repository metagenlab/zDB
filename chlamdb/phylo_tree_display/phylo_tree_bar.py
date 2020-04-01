#!/usr/bin/env python

from ete3 import Tree, SeqMotifFace, TreeStyle, add_face_to_node, TextFace, BarChartFace, StackedBarFace, NodeStyle, faces


def plot_tree_text_metadata(tree_file,
                            header2taxon2text,
                            ordered_header_list, 
                            biodb):

    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)

    t1 = Tree(tree_file)

    taxon2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb, filter_names=True)

    # Calculate the midpoint node
    R = t1.get_midpoint_outgroup()
    # and set it as tree outgroup
    t1.set_outgroup(R)
    tss = TreeStyle()
    tss.draw_guiding_lines = True
    tss.guiding_lines_color = "gray"
    tss.show_leaf_name = False

    for i, leaf in enumerate(t1.iter_leaves()):

        # first leaf, add headers
        if i == 0:
            for column, header in enumerate(ordered_header_list):

                n = TextFace('%s' % (header))
                n.margin_top = 0
                n.margin_right = 1
                n.margin_left = 20
                n.margin_bottom = 1
                n.rotation = 270
                n.hz_align = 2
                n.vt_align = 2
                n.inner_background.color = "white"
                n.opacity = 1.
                tss.aligned_header.add_face(n, column)
        for column, header in enumerate(ordered_header_list):
            text = header2taxon2text[header][int(leaf.name)]
            n = TextFace('%s' % text)
            n.margin_top = 1
            n.margin_right = 1
            n.margin_left = 5
            n.margin_bottom = 1
            n.inner_background.color = "white"
            n.opacity = 1.
            #n.rotation = 270
            leaf.add_face(n, column+1, position="aligned")
        # rename leaf (taxon_id => description)        
        n = TextFace(taxon2description[leaf.name], fgcolor = "black", fsize = 12, fstyle = 'italic')
        leaf.add_face(n, 0)


    for n in t1.traverse():
        # rename leaf

        nstyle = NodeStyle()

        if n.support < 1:
            nstyle["fgcolor"] = "black"
            nstyle["size"] = 6
            n.set_style(nstyle)
        else:
            nstyle["fgcolor"] = "red"
            nstyle["size"] = 0
            n.set_style(nstyle)



    return t1, tss


def plot_tree_stacked_barplot(tree_file,
                             taxon2value_list_barplot=False,
                             header_list=False, # header stackedbarplots
                             taxon2set2value_heatmap=False,
                             taxon2label=False,
                             header_list2=False, # header counts columns
                             biodb=False,
                             column_scale=True,
                             general_max=False,
                             header_list3 =False,
                             set2taxon2value_list_simple_barplot=False,
                             set2taxon2value_list_simple_barplot_counts=True,
                             rotate=False,
                             taxon2description=False):

    '''

    taxon2value_list_barplot list of lists:
    [[bar1_part1, bar1_part2,...],[bar2_part1, bar2_part2]]
    valeures de chaque liste transformes en pourcentages

    :param tree_file:
    :param taxon2value_list:
    :param biodb:
    :param exclude_outgroup:
    :param bw_scale:
    :return:
    '''



    if biodb:
        from chlamdb.biosqldb import manipulate_biosqldb
        server, db = manipulate_biosqldb.load_db(biodb)

        taxon2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb, filter_names=True)


    t1 = Tree(tree_file)

    # Calculate the midpoint node
    R = t1.get_midpoint_outgroup()
    # and set it as tree outgroup
    t1.set_outgroup(R)

    colors2 = ["red","#FFFF00","#58FA58","#819FF7","#F781F3", "#2E2E2E","#F7F8E0", 'black']
    colors = ["#7fc97f","#386cb0","#fdc086","#ffffb3","#fdb462", "#f0027f","#F7F8E0", 'black'] # fdc086ff 386cb0ff f0027fff

    tss = TreeStyle()
    tss.draw_guiding_lines = True
    tss.guiding_lines_color = "gray"
    tss.show_leaf_name = False
    if column_scale and header_list2:
        import matplotlib.cm as cm
        from matplotlib.colors import rgb2hex
        import matplotlib as mpl
        column2scale = {}
        col_n = 0
        for column in header_list2:
            values = taxon2set2value_heatmap[column].values()
            #print values
            if min(values) == max(values):
                min_val = 0
                max_val = 1.5*max(values)
            else:
                min_val = min(values)
                max_val = max(values)
            #print 'min-max', min_val, max_val
            norm = mpl.colors.Normalize(vmin=min_val,vmax=max_val) # *1.1
            if col_n < 4:
                cmap = cm.OrRd#
            else:
                cmap = cm.YlGnBu#PuBu#OrRd

            m = cm.ScalarMappable(norm=norm, cmap=cmap)

            column2scale[column] = [m, float(max_val)] # *0.7
            col_n+=1

    for i, lf in enumerate(t1.iter_leaves()):


        #if taxon2description[lf.name] == 'Pirellula staleyi DSM 6068':
        #    lf.name = 'Pirellula staleyi DSM 6068'
        #    continue
        if i==0:

            if taxon2label:
                n = TextFace('  ')
                n.margin_top = 1
                n.margin_right = 1
                n.margin_left = 20
                n.margin_bottom = 1
                n.hz_align = 2
                n.vt_align = 2
                n.rotation = 270
                n.inner_background.color = "white"
                n.opacity = 1.

                tss.aligned_header.add_face(n, 0)
                col_add=1
            else:
                col_add=1
            if header_list:
                for col, header in enumerate(header_list):

                    n = TextFace('%s' % (header))
                    n.margin_top = 0
                    n.margin_right = 1
                    n.margin_left = 20
                    n.margin_bottom = 1
                    n.rotation = 270
                    n.hz_align = 2
                    n.vt_align = 2
                    n.inner_background.color = "white"
                    n.opacity = 1.
                    tss.aligned_header.add_face(n, col+col_add)
                col_add += col+1

            if header_list3:
                #print 'header_list 3!'
                col_tmp=0
                for header in header_list3:
                    n = TextFace('%s' % (header))
                    n.margin_top = 1
                    n.margin_right = 1
                    n.margin_left = 20
                    n.margin_bottom = 1
                    n.rotation = 270
                    n.hz_align = 2
                    n.vt_align = 2
                    n.inner_background.color = "white"
                    n.opacity = 1.



                    if set2taxon2value_list_simple_barplot_counts:
                        if col_tmp==0:
                            col_tmp+=1
                        tss.aligned_header.add_face(n, col_tmp+1+col_add)
                        n = TextFace('       ')
                        tss.aligned_header.add_face(n, col_tmp+col_add)
                        col_tmp+=2
                    else:
                        tss.aligned_header.add_face(n, col_tmp+col_add)
                        col_tmp+=1
                if set2taxon2value_list_simple_barplot_counts:
                    col_add += col_tmp
                else:
                    col_add += col_tmp

            if header_list2:
                for col, header in enumerate(header_list2):
                    n = TextFace('%s' % (header))
                    n.margin_top = 1
                    n.margin_right = 1
                    n.margin_left = 20
                    n.margin_bottom = 1
                    n.rotation = 270
                    n.hz_align = 2
                    n.vt_align = 2
                    n.inner_background.color = "white"
                    n.opacity = 1.
                    tss.aligned_header.add_face(n, col+col_add)
                col_add += col+1




        if taxon2label:
            try:
                n = TextFace('%s' % taxon2label[lf.name])
            except:
                try:
                    n = TextFace('%s' % taxon2label[int(lf.name)])
                except:
                    n = TextFace('-')
            n.margin_top = 1
            n.margin_right = 1
            n.margin_left = 20
            n.margin_bottom = 1
            n.inner_background.color = "white"
            n.opacity = 1.
            if rotate:
               n.rotation= 270
            lf.add_face(n, 1, position="aligned")
            col_add=2
        else:
            col_add=2

        if taxon2value_list_barplot:

            try:
                val_list_of_lists = taxon2value_list_barplot[lf.name]
            except:
                val_list_of_lists = taxon2value_list_barplot[int(lf.name)]

            #col_count = 0
            for col, value_list in enumerate(val_list_of_lists):

                total = float(sum(value_list))
                percentages = [(i/total)*100 for i in value_list]
                if col % 3 == 0:
                    col_list = colors2
                else:
                    col_list = colors
                b = StackedBarFace(percentages,
                                    width=150, height=18, colors=col_list[0:len(percentages)])
                b.rotation= 0
                b.inner_border.color = "white"
                b.inner_border.width = 0
                b.margin_right = 5
                b.margin_left = 5
                if rotate:
                    b.rotation = 270
                lf.add_face(b, col+col_add, position="aligned")
                #col_count+=1

            col_add+=col+1

        if set2taxon2value_list_simple_barplot:
            col_list = ['#fc8d59', '#91bfdb', '#99d594', '#c51b7d', '#f1a340', '#999999']
            color_i=0
            col=0
            for one_set in header_list3:
                if color_i>5:
                    color_i=0
                color = col_list[color_i]
                color_i+=1
                # values for all taxons
                values_lists = [float(i) for i in set2taxon2value_list_simple_barplot[one_set].values()]
                #print values_lists
                #print one_set
                value = set2taxon2value_list_simple_barplot[one_set][lf.name]

                if set2taxon2value_list_simple_barplot_counts:
                    if isinstance(value, float):
                        a = TextFace(" %s " % str(round(value,2)))
                    else:
                        a = TextFace(" %s " % str(value))
                    a.margin_top = 1
                    a.margin_right = 2
                    a.margin_left = 5
                    a.margin_bottom = 1
                    if rotate:
                        a.rotation = 270
                    lf.add_face(a, col+col_add, position="aligned")


                #print 'value and max', value, max(values_lists)
                fraction_biggest = (float(value)/max(values_lists))*100
                fraction_rest = 100-fraction_biggest

                #print 'fractions', fraction_biggest, fraction_rest
                b = StackedBarFace([fraction_biggest, fraction_rest], width=100, height=15,colors=[color, 'white'])
                b.rotation= 0
                b.inner_border.color = "grey"
                b.inner_border.width = 0
                b.margin_right = 15
                b.margin_left = 0
                if rotate:
                    b.rotation = 270
                if set2taxon2value_list_simple_barplot_counts:
                    if col == 0:
                        col+=1
                    lf.add_face(b, col+1+col_add, position="aligned")
                    col+=2
                else:
                    lf.add_face(b, col+col_add, position="aligned")
                    col+=1
            if set2taxon2value_list_simple_barplot_counts:
                col_add+=col

            else:
                col_add+=col

        if taxon2set2value_heatmap:
            i = 0
            #if not taxon2label:
            #    col_add-=1
            for col2, head in enumerate(header_list2):

                col_name = header_list2[i]
                try:
                    value = taxon2set2value_heatmap[col_name][str(lf.name)]
                except:
                    try:
                        value = taxon2set2value_heatmap[col_name][round(float(lf.name),2)]
                    except:
                        value = 0
                if header_list2[i] == 'duplicates':
                    print ('dupli', lf.name, value)
                #print 'val----------------', value
                if int(value) > 0:
                    if int(value) >=10 and int(value) < 100:
                        n = TextFace('%4i' % value)
                    elif int(value)>=100:
                        n = TextFace('%3i' % value)
                    else:

                        n = TextFace('%5i' % value)

                    n.margin_top = 1
                    n.margin_right = 2
                    n.margin_left = 5
                    n.margin_bottom = 1
                    n.hz_align = 1
                    n.vt_align = 1
                    if rotate:
                        n.rotation = 270
                    n.inner_background.color = rgb2hex(column2scale[col_name][0].to_rgba(float(value)))#"orange"
                    #print 'xaxaxaxaxa', value,
                    if float(value) > column2scale[col_name][1]:

                        n.fgcolor = 'white'
                    n.opacity = 1.
                    n.hz_align = 1
                    n.vt_align = 1
                    lf.add_face(n, col2+col_add, position="aligned")
                    i+=1
                else:
                    n = TextFace('')
                    n.margin_top = 1
                    n.margin_right = 1
                    n.margin_left = 5
                    n.margin_bottom = 1
                    n.inner_background.color = "white"
                    n.opacity = 1.
                    if rotate:
                        n.rotation = 270
                    lf.add_face(n, col2+col_add, position="aligned")
                    i+=1

        #lf.name = taxon2description[lf.name]
        n = TextFace(taxon2description[lf.name], fgcolor = "black", fsize = 12, fstyle = 'italic')
        lf.add_face(n, 0)

    for n in t1.traverse():
       nstyle = NodeStyle()

       if n.support < 1:
           nstyle["fgcolor"] = "black"
           nstyle["size"] = 6
           n.set_style(nstyle)
       else:
           nstyle["fgcolor"] = "red"
           nstyle["size"] = 0
           n.set_style(nstyle)

    return t1, tss


def plot_tree_barplot(tree_file,
                      taxon2value_list_barplot,
                      header_list,
                      taxon2set2value_heatmap=False,
                      header_list2=False,
                      presence_only=True,
                      biodb="chlamydia_04_16",
                      column_scale=True,
                      general_max=False,
                      barplot2percentage=False):

    '''

    display one or more barplot

    :param tree_file:
    :param taxon2value_list:
    :param biodb:
    :param exclude_outgroup:
    :param bw_scale:
    :param barplot2percentage: list of bool to indicates if the number are percentages and the range should be set to 0-100

    :return:
    '''

    from chlamdb.biosqldb import manipulate_biosqldb
    import matplotlib.cm as cm
    from matplotlib.colors import rgb2hex
    import matplotlib as mpl

    server, db = manipulate_biosqldb.load_db(biodb)

    taxon2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb, filter_names=True)

    #print isinstance(tree_file, Tree)
    #print type(tree_file)

    if isinstance(tree_file, Tree):
       t1 = tree_file
    else:
       t1 = Tree(tree_file)

    # Calculate the midpoint node
    R = t1.get_midpoint_outgroup()
    # and set it as tree outgroup
    t1.set_outgroup(R)


    tss = TreeStyle()
    value=1
    tss.draw_guiding_lines = True
    tss.guiding_lines_color = "gray"
    tss.show_leaf_name = False


    if column_scale and header_list2:
        import matplotlib.cm as cm
        from matplotlib.colors import rgb2hex
        import matplotlib as mpl
        column2scale = {}
        for column in header_list2:
            values = taxon2set2value_heatmap[column].values()

            norm = mpl.colors.Normalize(vmin=min(values), vmax=max(values))
            cmap = cm.OrRd
            m = cm.ScalarMappable(norm=norm, cmap=cmap)
            column2scale[column] = m

    cmap = cm.YlGnBu#YlOrRd#OrRd

    values_lists = taxon2value_list_barplot.values()

    scale_list = []
    max_value_list = []

    for n, header in enumerate(header_list):
        #print 'scale', n, header
        data = [float(i[n]) for i in values_lists]

        if barplot2percentage is False:
            max_value = max(data)#3424182#
            min_value = min(data) #48.23
        else:
            if barplot2percentage[n] is True:
                max_value = 100
                min_value = 0
            else:
                max_value = max(data)#3424182#
                min_value = min(data) #48.23
        norm = mpl.colors.Normalize(vmin=min_value, vmax=max_value)
        m1 = cm.ScalarMappable(norm=norm, cmap=cmap)
        scale_list.append(m1)
        if not general_max:
            max_value_list.append(float(max_value))
        else:
            max_value_list.append(general_max)

    for i, lf in enumerate(t1.iter_leaves()):

        #if taxon2description[lf.name] == 'Pirellula staleyi DSM 6068':
        #    lf.name = 'Pirellula staleyi DSM 6068'
        #    continue
        if i==0:

            col_add=0
            for col, header in enumerate(header_list):

                #lf.add_face(n, column, position="aligned")
                n = TextFace(' ')
                n.margin_top = 1
                n.margin_right = 2
                n.margin_left = 2
                n.margin_bottom = 1
                n.rotation = 90
                n.inner_background.color = "white"
                n.opacity = 1.
                n.hz_align = 2
                n.vt_align = 2

                tss.aligned_header.add_face(n, col_add)


                n = TextFace('%s' % header)
                n.margin_top = 1
                n.margin_right = 2
                n.margin_left = 2
                n.margin_bottom = 80
                n.rotation = 270
                n.inner_background.color = "white"
                n.opacity = 1.
                n.hz_align = 2
                n.vt_align = 2
                tss.aligned_header.add_face(n, col_add+1)
                col_add+=2

            if header_list2:
                for col, header in enumerate(header_list2):
                    n = TextFace('%s' % header)
                    n.margin_top = 1
                    n.margin_right = 20
                    n.margin_left = 2
                    n.margin_bottom = 1
                    n.rotation = 270
                    n.hz_align = 2
                    n.vt_align = 2
                    n.inner_background.color = "white"
                    n.opacity = 1.
                    tss.aligned_header.add_face(n, col+col_add)

        try:
            val_list = taxon2value_list_barplot[lf.name]
        except:
            try:
                val_list = taxon2value_list_barplot[int(lf.name)]
            except:
                val_list = [0]
        col_add=0
        for col, value in enumerate(val_list):

            # show value itself
            try:
                n = TextFace('  %s  ' % str(value))
            except:
                n = TextFace('  %s  ' % str(value))
            n.margin_top = 1
            n.margin_right = 10
            n.margin_left = 2
            n.margin_bottom = 1
            n.inner_background.color = "white"
            n.opacity = 1.

            lf.add_face(n, col_add, position="aligned")
            # show bar

            color = rgb2hex(scale_list[col].to_rgba(float(value)))
            try:
                percentage = (value/max_value_list[col])*100
                #percentage = value
            except:
                percentage = 0
            maximum_bar = ((max_value_list[col]-value)/max_value_list[col])*100
            #maximum_bar = 100-percentage
            b = StackedBarFace([percentage,
                                maximum_bar],
                                width=100, height=10, colors=[color, "white"])
            b.rotation= 0
            b.inner_border.color = "grey"
            b.inner_border.width = 0
            b.margin_right = 15
            b.margin_left = 0
            lf.add_face(b, col_add+1, position="aligned")
            col_add+=2


        if taxon2set2value_heatmap:
            shift = col+col_add+1

            i = 0
            for col, col_name in enumerate(header_list2):
                try:
                    value = taxon2set2value_heatmap[col_name][lf.name]
                except:
                    try:
                        value = taxon2set2value_heatmap[col_name][int(lf.name)]
                    except:
                        value = 0

                if int(value) >0:
                    if int(value)>9:
                        n = TextFace(' %i ' % int(value))
                    else:
                        n = TextFace(' %i   ' % int(value))
                    n.margin_top = 1
                    n.margin_right = 1
                    n.margin_left = 20
                    n.margin_bottom = 1
                    n.fgcolor = "white"
                    n.inner_background.color = rgb2hex(column2scale[col_name].to_rgba(float(value)))#"orange"
                    n.opacity = 1.
                    lf.add_face(n, col+col_add, position="aligned")
                    i+=1
                else:
                    n = TextFace('  ') #% str(value))
                    n.margin_top = 1
                    n.margin_right = 1
                    n.margin_left = 20
                    n.margin_bottom = 1
                    n.inner_background.color = "white"
                    n.opacity = 1.

                    lf.add_face(n, col+col_add, position="aligned")



        #lf.name = taxon2description[lf.name]
        try:
            n = TextFace(taxon2description[lf.name], fgcolor = "black", fsize = 12, fstyle = 'italic')
        except:
            n = TextFace(lf.name, fgcolor = "black", fsize = 12, fstyle = 'italic')
        lf.add_face(n, 0)

    for n in t1.traverse():
       nstyle = NodeStyle()
       if n.support < 1:
           nstyle["fgcolor"] = "black"
           nstyle["size"] = 6
           n.set_style(nstyle)
       else:
           nstyle["fgcolor"] = "red"
           nstyle["size"] = 0
           n.set_style(nstyle)
    #print t1
    return t1, tss

def plot_heat_tree(tree_file, biodb="chlamydia_04_16", exclude_outgroup=False, bw_scale=True):
    from chlamdb.biosqldb import manipulate_biosqldb
    import matplotlib.cm as cm
    from matplotlib.colors import rgb2hex
    import matplotlib as mpl

    server, db = manipulate_biosqldb.load_db(biodb)

    sql_biodatabase_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb
    db_id = server.adaptor.execute_and_fetchall(sql_biodatabase_id,)[0][0]
    if type(tree_file) == str:
        t1 = Tree(tree_file)
        try:
            R = t1.get_midpoint_outgroup()
            #print 'root', R
            # and set it as tree outgroup
            t1.set_outgroup(R)
        except:
            pass
    elif isinstance(tree_file, Tree):
        t1 = tree_file
    else:
        IOError('Unkown tree format')
    tss = TreeStyle()
    tss.draw_guiding_lines = True
    tss.guiding_lines_color = "gray"
    tss.show_leaf_name = False

    #print "tree", t1

    sql1 = 'select taxon_id, description from bioentry where biodatabase_id=%s and description not like "%%%%plasmid%%%%"' % db_id
    sql2 = 'select t2.taxon_id, t1.GC from genomes_info as t1 inner join bioentry as t2 ' \
           ' on t1.accession=t2.accession where t2.biodatabase_id=%s and t1.description not like "%%%%plasmid%%%%";' % (biodb, db_id)
    sql3 = 'select t2.taxon_id, t1.genome_size from genomes_info as t1 ' \
           ' inner join bioentry as t2 on t1.accession=t2.accession ' \
           ' where t2.biodatabase_id=%s and t1.description not like "%%%%plasmid%%%%";' % (biodb, db_id)
    sql4 = 'select t2.taxon_id,percent_non_coding from genomes_info as t1 ' \
           ' inner join bioentry as t2 on t1.accession=t2.accession ' \
           ' where t2.biodatabase_id=%s and t1.description not like "%%%%plasmid%%%%";' % (biodb, db_id)
           
    sql_checkm_completeness = 'select taxon_id, completeness from custom_tables_checkm;' % biodb
    sql_checkm_contamination = 'select taxon_id,contamination from custom_tables_checkm;' % biodb

    try:
        taxon_id2completeness = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_checkm_completeness))
        taxon_id2contamination = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_checkm_contamination))
    except:
        taxon_id2completeness = False
    #taxon2description = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql1,))

    taxon2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb, filter_names=True)

    taxon2gc = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
    taxon2genome_size = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3,))
    taxon2coding_density = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql4,))


    my_taxons = [lf.name for lf in t1.iter_leaves()]

    # Calculate the midpoint node

    if exclude_outgroup:
        excluded = str(list(t1.iter_leaves())[0].name)
        my_taxons.pop(my_taxons.index(excluded))


    genome_sizes = [float(taxon2genome_size[i]) for i in my_taxons]
    gc_list = [float(taxon2gc[i]) for i in my_taxons]
    fraction_list = [float(taxon2coding_density[i]) for i in my_taxons]


    value=1

    max_genome_size = max(genome_sizes)#3424182#
    max_gc = max(gc_list) #48.23

    cmap = cm.YlGnBu#YlOrRd#OrRd

    norm = mpl.colors.Normalize(vmin=min(genome_sizes)-100000, vmax=max(genome_sizes))
    m1 = cm.ScalarMappable(norm=norm, cmap=cmap)
    norm = mpl.colors.Normalize(vmin=min(gc_list), vmax=max(gc_list))
    m2 = cm.ScalarMappable(norm=norm, cmap=cmap)
    norm = mpl.colors.Normalize(vmin=min(fraction_list), vmax=max(fraction_list))
    m3 = cm.ScalarMappable(norm=norm, cmap=cmap)


    for i, lf in enumerate(t1.iter_leaves()):
        #if taxon2description[lf.name] == 'Pirellula staleyi DSM 6068':
        #    lf.name = 'Pirellula staleyi DSM 6068'
        #    continue
        if i==0:
            n = TextFace('Size (Mbp)')
            n.rotation = -25
            n.margin_top = 1
            n.margin_right = 1
            n.margin_left = 20
            n.margin_bottom = 1
            n.inner_background.color = "white"
            n.opacity = 1.
            #lf.add_face(n, 3, position="aligned")
            tss.aligned_header.add_face(n, 3)
            n = TextFace('GC (%)')
            n.rotation = -25
            n.margin_top = 1
            n.margin_right = 1
            n.margin_left = 20
            n.margin_bottom = 1
            n.inner_background.color = "white"
            n.opacity = 1.
            #lf.add_face(n, 5, position="aligned")
            tss.aligned_header.add_face(n, 5)
            n = TextFace('')
            #lf.add_face(n, 2, position="aligned")
            tss.aligned_header.add_face(n, 2)
            #lf.add_face(n, 4, position="aligned")
            tss.aligned_header.add_face(n, 4)
            n = TextFace('Non coding (%)')
            n.margin_top = 1
            n.margin_right = 1
            n.margin_left = 20
            n.margin_bottom = 1
            n.inner_background.color = "white"
            n.opacity = 1.
            n.rotation = -25
            #lf.add_face(n, 7, position="aligned")
            tss.aligned_header.add_face(n, 7)
            n = TextFace('')
            #lf.add_face(n, 6, position="aligned")
            tss.aligned_header.add_face(n, 6)
            
            if taxon_id2completeness:
                n = TextFace('Completeness (%)')
                n.margin_top = 1
                n.margin_right = 1
                n.margin_left = 20
                n.margin_bottom = 1
                n.inner_background.color = "white"
                n.opacity = 1.
                n.rotation = -25
                #lf.add_face(n, 7, position="aligned")
                tss.aligned_header.add_face(n, 9)
                n = TextFace('')
                #lf.add_face(n, 6, position="aligned")
                tss.aligned_header.add_face(n, 8) 

                n = TextFace('Contamination (%)')
                n.margin_top = 1
                n.margin_right = 1
                n.margin_left = 20
                n.margin_bottom = 1
                n.inner_background.color = "white"
                n.opacity = 1.
                n.rotation = -25
                #lf.add_face(n, 7, position="aligned")
                tss.aligned_header.add_face(n, 11)
                n = TextFace('')
                #lf.add_face(n, 6, position="aligned")
                tss.aligned_header.add_face(n, 10) 


        value+=1

        #print '------ %s' % lf.name
        if exclude_outgroup and i == 0:
            lf.name = taxon2description[lf.name]
            #print '#######################'
            continue

        n = TextFace('  %s ' % str(round(taxon2genome_size[lf.name]/float(1000000),2)))
        n.margin_top = 1
        n.margin_right = 1
        n.margin_left = 0
        n.margin_bottom = 1
        n.fsize = 7
        n.inner_background.color = "white"
        n.opacity = 1.

        lf.add_face(n, 2, position="aligned")
        #if max_genome_size > 3424182:
        #    max_genome_size = 3424182
        fraction_biggest = (float(taxon2genome_size[lf.name])/max_genome_size)*100
        fraction_rest = 100-fraction_biggest
        if taxon2description[lf.name] == 'Rhabdochlamydia helveticae T3358':
            col = '#fc8d59'
        else:
            if not bw_scale:
                col = rgb2hex(m1.to_rgba(float(taxon2genome_size[lf.name])))  # 'grey'
            else:
                col = '#fc8d59'

        b = StackedBarFace([fraction_biggest, fraction_rest], width=100, height=9,colors=[col, 'white'])
        b.rotation= 0
        b.inner_border.color = "black"
        b.inner_border.width = 0
        b.margin_right = 15
        b.margin_left = 0
        lf.add_face(b, 3, position="aligned")

        fraction_biggest = (float(taxon2gc[lf.name])/max_gc)*100
        fraction_rest = 100-fraction_biggest
        if taxon2description[lf.name] == 'Rhabdochlamydia helveticae T3358':
            col = '#91bfdb'
        else:
            if not bw_scale:
                col = rgb2hex(m2.to_rgba(float(taxon2gc[lf.name])))
            else:
                col = '#91bfdb'
        b = StackedBarFace([fraction_biggest, fraction_rest], width=100, height=9,colors=[col, 'white'])
        b.rotation= 0
        b.inner_border.color = "black"
        b.inner_border.width = 0
        b.margin_left = 0
        b.margin_right = 15


        lf.add_face(b, 5, position="aligned")
        n = TextFace('  %s ' % str(round(float(taxon2gc[lf.name]),2)))
        n.margin_top = 1
        n.margin_right = 0
        n.margin_left = 0
        n.margin_bottom = 1
        n.fsize = 7
        n.inner_background.color = "white"
        n.opacity = 1.
        lf.add_face(n, 4, position="aligned")


        if taxon2description[lf.name] == 'Rhabdochlamydia helveticae T3358':
            col = '#99d594'
        else:
            if not bw_scale:
                col = rgb2hex(m3.to_rgba(float(taxon2coding_density[lf.name])))
            else:
                col = '#99d594'
        n = TextFace('  %s ' % str(float(taxon2coding_density[lf.name])))
        n.margin_top = 1
        n.margin_right = 0
        n.margin_left = 0
        n.margin_right = 0
        n.margin_bottom = 1
        n.fsize = 7
        n.inner_background.color = "white"
        n.opacity = 1.
        lf.add_face(n, 6, position="aligned")
        fraction = (float(taxon2coding_density[lf.name])/max(taxon2coding_density.values()))*100
        fraction_rest = ((max(taxon2coding_density.values()) - taxon2coding_density[lf.name])/float(max(taxon2coding_density.values())))*100
        #print 'fraction, rest', fraction, fraction_rest
        b = StackedBarFace([fraction, fraction_rest], width=100, height=9,colors=[col, 'white'])# 1-round(float(taxon2coding_density[lf.name]), 2)
        b.rotation = 0
        b.margin_right = 1
        b.inner_border.color = "black"
        b.inner_border.width = 0
        b.margin_left = 5
        lf.add_face(b, 7, position="aligned")
        
        
        if taxon_id2completeness:
            n = TextFace('  %s ' % str(float(taxon_id2completeness[lf.name])))
            n.margin_top = 1
            n.margin_right = 0
            n.margin_left = 0
            n.margin_right = 0
            n.margin_bottom = 1
            n.fsize = 7
            n.inner_background.color = "white"
            n.opacity = 1.
            lf.add_face(n, 8, position="aligned")
            fraction = float(taxon_id2completeness[lf.name])
            fraction_rest = 100-fraction
            #print 'fraction, rest', fraction, fraction_rest
            b = StackedBarFace([fraction, fraction_rest], width=100, height=9,colors=["#d7191c", 'white'])# 1-round(float(taxon2coding_density[lf.name]), 2)
            b.rotation = 0
            b.margin_right = 1
            b.inner_border.color = "black"
            b.inner_border.width = 0
            b.margin_left = 5
            lf.add_face(b, 9, position="aligned")
            
            
            n = TextFace('  %s ' % str(float(taxon_id2contamination[lf.name])))
            n.margin_top = 1
            n.margin_right = 0
            n.margin_left = 0
            n.margin_right = 0
            n.margin_bottom = 1
            n.fsize = 7
            n.inner_background.color = "white"
            n.opacity = 1.
            lf.add_face(n, 10, position="aligned")
            fraction = float(taxon_id2contamination[lf.name])
            fraction_rest = 100-fraction
            #print 'fraction, rest', fraction, fraction_rest
            b = StackedBarFace([fraction, fraction_rest], width=100, height=9,colors=["black", 'white'])# 1-round(float(taxon2coding_density[lf.name]), 2)
            b.rotation = 0
            b.margin_right = 1
            b.inner_border.color = "black"
            b.inner_border.width = 0
            b.margin_left = 5
            lf.add_face(b, 11, position="aligned")
            
            
        
                #lf.name = taxon2description[lf.name]
        n = TextFace(taxon2description[lf.name], fgcolor = "black", fsize = 9, fstyle = 'italic')
        n.margin_right = 30
        lf.add_face(n, 0)

    for n in t1.traverse():
       nstyle = NodeStyle()
       if n.support < 1:
           nstyle["fgcolor"] = "black"
           nstyle["size"] = 6
           n.set_style(nstyle)
       else:
           nstyle["fgcolor"] = "red"
           nstyle["size"] = 0
           n.set_style(nstyle)

    return t1, tss

if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    parser = argparse.ArgumentParser()
    parser.add_argument("-t",'--tree',type=str,help="newick tree")
    parser.add_argument("-m",'--matrix',type=str,help="matrix (tab file)")
    parser.add_argument("-s",'--mlst',type=str,help="mlst file")
    parser.add_argument("-d",'--biodb',type=str,help="biodatabase name")
    parser.add_argument("-e",'--exclude_outgroup',action="store_true",help="exclude outroup from annotation")
    parser.add_argument("-w",'--bw_scale',action="store_true",help="use black and white color scale")
    parser.add_argument("-a",'--accession2description',default=False, type=str,help="tab file with accessions and corresponding descriptions for leaf labels")

    args = parser.parse_args()

    if not args.tree:
        from chlamdb.biosqldb import manipulate_biosqldb
        sql_tree = 'select tree from reference_phylogeny t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                   ' where t2.name="%s";' % args.biodb
        server, db = manipulate_biosqldb.load_db(args.biodb)
        tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]
        #print tree
    else:
        tree = args.tree

    t = plot_heat_tree(tree, args.biodb, args.exclude_outgroup, args.bw_scale)
    ts = TreeStyle()
    ts.show_branch_support = False
    t.render("test2.svg", tree_style=ts)
