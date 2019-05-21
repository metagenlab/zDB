#!/usr/bin/env python

import urllib2


def pdb_entry2data(pdb_entry):

    import xmltodict
    from xml.etree import ElementTree

    link = "http://www.rcsb.org/pdb/rest/describeMol?structureId=%s" % (pdb_entry)
    req = urllib2.Request(link)
    try:
        page = urllib2.urlopen(req)
    except:
        print 'echec: %s' % pdb_entry

    tree = ElementTree.fromstring(page.read().decode('utf-8'))
    #print tree
    #root = tree.getroot()
    #print root
    for elem in tree.iter(tag='Taxonomy'):
        print "%s\t%s\t%s" % (pdb_entry,elem.attrib['name'], elem.attrib['id'])





if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i",'--pdb_id', default=False, type=str, help="pdb accession")

    args = parser.parse_args()
    pdb_entry2data(args.pdb_id)