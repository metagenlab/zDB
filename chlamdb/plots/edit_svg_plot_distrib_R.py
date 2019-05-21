#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# convert embl file to gbk
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: mars 1015
# ---------------------------------------------------------------------------

def edit_svg_file(svg_file):

    '''
    Test function to edit a R svg graph

    '''
    
    from xml.dom import minidom
    import xml.etree
    from xml.etree import ElementTree
    tree = ElementTree.parse(svg_file)
    from lxml import etree
    #tree = etree.parse(svg_handle)

    #locus_list=locus_list[::-1]

    i = 0

    for element in tree.iter():
        xmlRoot = tree.getroot()

        if element.tag.split("}")[1] == 'path':
            print element.attrib.keys()
            try:
                print 'len', len(element.attrib['d'])
                if len(element.attrib['d']) > 10000:
                    print element
                    element.attrib['class'] = 'line'
                    element.attrib['id'] = "line_%s" % i
                    i+=1
                    print i
                    print element.tag

            except:
                pass
            '''
            if len(element.get("points").split(","))==7:
                print element.get("points")
                # window.location.href
                add = 'window.open("/chlamdb/locusx/%s/%s", "theFrame");' % (biodb_name, locus_list[i])
                #add = 'window.location.href="/chlamdb/locusx/%s/%s#target_div;' % (biodb_name, locus_list[i])
                add2= "_blank"
                mystyle = element.get("style")
                add3 = mystyle + "cursor: pointer;"
                add4 = "this.style.stroke = '#ff0000'; this.style['stroke-width'] = 5;"
                add5 = "this.style.stroke = '#000000'; this.style['stroke-width'] = 0;"
                print add3
                i+=1
                element.set("onclick", add)
                element.set("target", add2)
                element.set("onmouseover", add4)
                element.set("onmouseout", add5)

                #nav = element.find('name')
                #element.addprevious(element_before)
                #element.addnext(element_after)
    # tree.write('output.svg')


    child = ElementTree.Element("text")
    print child

    child.attrib['x'] = "10"
    child.attrib['y'] = "10"
    child.attrib["class"] = 'legend'
    child.attrib["stroke"] = "red"
    fct = 'function(){' \
         ' var active   = redLine.active ? false : true ,' \
         '   newOpacity = active ? 0 : 1;' \
         '   d3.select("#line_1").style("opacity", newOpacity);' \
         '   line_1.active = active;})'
    child.set("onclick", fct)
    child.text = 'Red Line'
    root = tree.getroot()
    root.append(child)
     '''
    tree.write('output.svg')
    return tree


'''
// Add the red line title
svg.append("path")
	.attr("class", "line")
	.style("stroke", "red")
	.attr("id", "redLine")
	.attr("d", valueline2(data));


svg.append("text")
	.attr("x", 0)
	.attr("y", height + margin.top + 30)
	.attr("class", "legend")
	.style("fill", "red")
	.on("click", function(){
		// Determine if current line is visible
		var active   = redLine.active ? false : true ,
		  newOpacity = active ? 0 : 1;
		// Hide or show the elements
		d3.select("#redLine").style("opacity", newOpacity);
		d3.select("#redAxis").style("opacity", newOpacity);
		// Update whether or not the elements are active
		redLine.active = active;
	})
	.text("Red Line");
'''

if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_svg', type=str, help="input svg file")



    args = parser.parse_args()
    import StringIO
    import edit_svg

    #f = open(args.input_svg, "rw")
    edit_svg_file(args.input_svg)
