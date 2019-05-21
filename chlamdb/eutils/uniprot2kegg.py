#!/usr/bin/env python

from Bio import Entrez


Entrez.email = "trestan.pillonel@unil.ch"

def uniprot2kegg(uniprot_id_list):
    import urllib,urllib2

    url = 'http://www.uniprot.org/uploadlists/'

    params = {
    'from':'ACC',
    'to':'	KEGG_ID',
    'format':'tab',
    'query':'%s' % (' '.join(uniprot_id_list))
    }

    data = urllib.urlencode(params)
    request = urllib2.Request(url, data)
    contact = "trestan.pillonel@unil.ch" # Please set your email address here to help us debug in case of problems.
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib2.urlopen(request)
    page = [i for i in response.read(200000).split('\n')]
    accession2kegg_list = {}
    for n, line in enumerate(page):
        if n != 0:
            data = line.split('\t')
            if len(data)>1:
                if data[0] not in accession2kegg_list:
                    accession2kegg_list[data[0]] = [data[1]]
                else:
                    accession2kegg_list[data[0]].append(data[1])
    return accession2kegg_list
def kegg2ko_id(kegg_id):
    import urllib2
    url = 'http://rest.kegg.jp/link/ko/%s' % kegg_id
    req = urllib2.Request(url)
    page = urllib2.urlopen(req)
    data = page.read().decode('utf-8').split('\n')
    return data[0].split('\t')[1].split(':')[1]



if __name__ == '__main__':
    import argparse
    import sys
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",'--uniprot_accession', help="uniprot id", nargs="+")

    args = parser.parse_args()


    dico = uniprot2kegg(args.uniprot_accession)
    for accession in dico:
        print "%s\t%s\t%s" % (accession, dico[accession][0], kegg2ko_id(dico[accession][0]))