#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# convert embl file to gbk
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: mars 1015
# ---------------------------------------------------------------------------

def edit_svg(svg_string, locus_list, biodb_name):
    print(len(locus_list))
    from xml.dom import minidom
    from xml.etree import ElementTree
    tree = ElementTree.ElementTree(ElementTree.fromstring(svg_string))
    from lxml import etree
    #tree = etree.parse(svg_handle)

    #locus_list=locus_list[::-1]

    i = 0
    #print locus_list
    for element in tree.iter():

        if element.tag.split("}")[1] == 'polygon':
            if len(element.get("points").split(",")) == 7 or len(element.get("points").split(",")) == 5:
                #print element.get("points")
                # window.location.href
                try:
                    # , "theFrame"
                    add = 'window.open("/chlamdb/locusx/%s/True", "_top");' % (locus_list[i])
                except IndexError:
                    add = 'window.open("/chlamdb/locusx/%s/True", "_top");' % (locus_list[0])
                #add = 'window.location.href="/chlamdb/locusx/%s/%s#target_div;' % (biodb_name, locus_list[i])
                #add2= "_blank"
                mystyle = element.get("style")
                add3 = mystyle + "cursor: pointer;"
                add4 = "this.style.stroke = '#ff0000'; this.style['stroke-width'] = 5;"
                add5 = "this.style.stroke = '#000000'; this.style['stroke-width'] = 0;"

                i+=1
                element.set("onclick", add)
                #element.set("target", add2)
                element.set("onmouseover", add4)
                element.set("onmouseout", add5)

                #nav = element.find('name')
                #element.addprevious(element_before)
                #element.addnext(element_after)
    # tree.write('output.svg')
    return tree

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_svg', type=str, help="input svg file")



    args = parser.parse_args()

    #f = open(args.input_svg, "rw")
    edit_svg(args.input_svg)
