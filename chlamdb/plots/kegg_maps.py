#!/usr/bin/python


def map2highlighted_map(map_id, ko_list, ko2freq, biodb, outpath = 'test.pdf', taxon_id=False, n_species=60):
    import re
    from chlamdb.biosqldb import shell_command
    from Bio.Graphics.KGML_vis import KGMLCanvas
    from Bio.Graphics import KGML_vis
    import urllib.request
    from Bio.KEGG.KGML.KGML_pathway import Pathway, Reaction, Relation
    import Bio.KEGG.KGML.KGML_pathway
    from Bio.KEGG.KGML import KGML_parser
    from Bio.Graphics.ColorSpiral import ColorSpiral
    import matplotlib.cm as cm
    from matplotlib.colors import rgb2hex
    import matplotlib as mpl

    values = [float(i) for i in ko2freq.values()]

    norm = mpl.colors.Normalize(vmin=0, vmax=n_species)
    cmap = cm.OrRd
    cmap2 = cm.Greens
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    m2 = cm.ScalarMappable(norm=norm, cmap=cmap2)

    url_template = 'http://rest.kegg.jp/get/%s/kgml' % re.sub('map', 'ko', map_id)
    print(url_template)
    f = urllib.request.urlopen(url_template)
    from Bio.Graphics import KGML_vis


    pathway = KGML_parser.read(f.read().decode('UTF-8'))

    kgml_map = KGMLCanvas(pathway, show_maps=True)

    # Let's use some arbitrary colours for the orthologs
    cs = ColorSpiral(a=2, b=0.2, v_init=0.85, v_final=0.5,
                     jitter=0.03)
    # Loop over the orthologs in the pathway, and change the
    # background colour
    orthologs = [e for e in pathway.orthologs]
    for o in orthologs:
        match = False
        if 'K00163' in o.name:
            print ('##################################')
        ko_temp_list = set([i.rstrip() for i in o.name.split('ko:') ])
        if len(ko_temp_list.intersection(set(ko2freq.keys())))>0:

            ko_keep = []
            for ko in ko_temp_list:
                if ko in ko2freq:
                    ko_keep.append(ko)
                if ko in ko_list:
                    match = True
            o.name = 'ko:'+' ko:'.join(ko_keep)
            total = sum([int(ko2freq[i]) for i in ko_temp_list.intersection(set(ko2freq.keys()))])

            for g in o.graphics:
                if match:
                    g.bgcolor = rgb2hex(m2.to_rgba(float(total)))
                else:
                    #print 'no match!!!!'
                    #print ko_temp_list
                    #print ko2freq.keys()
                    #print 'TOTAL:', total
                    g.bgcolor = rgb2hex(m.to_rgba(float(total)))
            o.name = "%s (%s)" % (o.name.split('ko:')[0], total)
        #else:
        #    for g in o.graphics:
        #        g.bgcolor = '#FFFFFF'

    # Default settings are for the KGML elements only


    # We need to use the image map, and turn off the KGML elements, to see
    # only the .png base map. We could have set these values on canvas
    # instantiation
    kgml_map.import_imagemap = True
    kgml_map.show_maps = True
    kgml_map.show_orthologs = True
    kgml_map.draw_relations = False
    kgml_map.show_compounds = False
    kgml_map.show_genes = False
    kgml_map.show_compounds = False
    kgml_map.show_genes = False
    kgml_map.draw(outpath)

    '''
    print 'DIRLISAT:', dir(pathway)
    maps = [m for m in pathway.maps]
    for map in maps:
        for g in map.graphics:
            print g.name
    '''

    #print re.sub('pdf', 'svg', outpath)
    shell_command.shell_command('inkscape %s --export-plain-svg=%s' % (outpath, re.sub('pdf', 'svg', outpath))) # 'pdf2svg %s %s all'
    t = edit_svg_map("%s" % re.sub('pdf', 'svg', outpath), ko2freq.keys(), biodb, map_id, taxon_id=taxon_id)
    #print "%s" % re.sub('pdf', 'svg', outpath)
    t.write("%s" % re.sub('pdf', 'svg', outpath))

#map2highlighted_map('ko00010',['K00001','K00162'])


def edit_svg_map(map_path, keep_ko_list, biodb_name, map_name,taxon_id=False):

    from chlamdb.biosqldb import manipulate_biosqldb
    import re

    server, db = manipulate_biosqldb.load_db(biodb_name)

    sql = 'select description,pathway_name from enzyme_kegg_pathway;'

    description2map = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    from xml.etree import ElementTree
    tree = ElementTree.parse(map_path)
    #print tree
    for element in tree.iter():
        if element.tag.split("}")[1] == 'text':
            #print element.tag
            #print element.attrib

            for child in element:
                #print child.tag
                #print child.attrib
                if child.text[0] != 'K':
                    #print child.text
                    try:
                        if not taxon_id:
                            add = 'window.open("/chlamdb/KEGG_mapp_ko/%s", "_top");' % ( description2map[child.text])
                        else:
                            add = 'window.open("/chlamdb/KEGG_mapp_ko_organism/%s/%s", "_top");' % (description2map[child.text], taxon_id)
                    except:
                        continue
                    mystyle = element.get("style")

                    add4 = "this.style.stroke = '#ff0000'; this.style['stroke-width'] = 1;"
                    add5 = "this.style.stroke = '#000000'; this.style['stroke-width'] = 0;"

                    element.set("onclick", add)
                    #element.set("target", add2)
                    element.set("onmouseover", add4)
                    element.set("onmouseout", add5)

                if child.text in keep_ko_list:
                    #print 'match-----------'
                    add = 'window.open("/chlamdb/fam/%s/ko", "_top");' % (child.text)
                    mystyle = element.get("style")

                    add4 = "this.style.stroke = '#ff0000'; this.style['stroke-width'] = 1;"
                    add5 = "this.style.stroke = '#000000'; this.style['stroke-width'] = 0;"

                    element.set("onclick", add)
                    #element.set("target", add2)
                    element.set("onmouseover", add4)
                    element.set("onmouseout", add5)
                if '...' in child.text:
                    #print 'match-----------'
                    add = 'window.open("/chlamdb/kegg_multi/%s/%s/", "_top");' % (map_name,
                                                                                     re.sub('\.\.\.', '', child.text))
                    mystyle = element.get("style")

                    add4 = "this.style.stroke = '#ff0000'; this.style['stroke-width'] = 1;"
                    add5 = "this.style.stroke = '#000000'; this.style['stroke-width'] = 0;"

                    element.set("onclick", add)
                    #element.set("target", add2)
                    element.set("onmouseover", add4)
                    element.set("onmouseout", add5)
    return tree
