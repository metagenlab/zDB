#!/usr/bin/env python

def process_tag(tag):
    return tag.split('}')[-1]

def accession2family(sequence, target_databases):
    import urllib2
    from bs4 import BeautifulSoup
    from xml.etree import cElementTree
    from xml.etree import ElementTree
    import re
    from urllib2 import URLError

    database_string = '&database=' .join(target_databases)

    link = "http://www.ebi.ac.uk/Tools/picr/rest/getUPIForSequence?sequence=%s&database=%s&includeattributes=true" % (sequence,
                                                                                                                      database_string)
    print link
    ns = 'http://model.picr.ebi.ac.uk'
    req = urllib2.Request(link)
    try:
        page = urllib2.urlopen(req)
        #data = page.read().decode('utf-8').split('\n')
        tree = ElementTree.parse(page)
    except:
        import time
        print 'connexion problem, trying again...'
        time.sleep(60)

    accession = ''
    version = ''
    databaseName = ''
    taxonId = ''

    #for element in parser:
    #    for event, element in parser:
    #        print event, process_tag(element.tag), element.text
    #for i in soup:
    #    print i
    print 'tree', tree.findall(ns+'identicalCrossReferences')
    for item in tree.findall('identicalCrossReferences'):
        print item

accession2family('MTDPLAPLASLPGVSDAAEAARDALSKVHRHRANLRGWPVTAAEAAVRAARSSSALDGGTMKLSADGAVEDPILAGALRVGQALDGDALAQLAAVWSRAPLQALARLHVLAATGMADEDTLGRPRPGVDTDRLELLAQLISGGTSVPAPILAAVTHGELLALNPFGSADGVVARAASRLVTVSRGLDPHGLGVPEVMWMRQSGRYRELSSAFAQGAPEAITAWIVFCCQALTAGAAEATSIADTAAG',
                 ['REFSEQ'])
