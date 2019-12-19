#!/usr/bin/env python

import urllib2

def string_id2pubmed_id_list(accession):

    link = 'http://string-db.org/api/tsv/abstractsList?identifiers=%s' % accession
    print link
    try:
        data = urllib2.urlopen(link).read().rstrip().decode('utf-8').split('\n')[1:]
    except urllib2.URLError:
        print ('echec', link)
        return False
    pid_list = [row.split(':')[1] for row in data]
    print ('list', pid_list)
    return pid_list

o = open("string_mapping_PMID.tab", "w")

string_mapping = "string_mapping.tab"

with open(string_mapping, 'r') as f:
    for n, row in enumerate(f):
        if n == 0:
            continue
        else:
            data = row.rstrip().split("	")
            pmid_list = string_id2pubmed_id_list(data[1])
            if pmid_list:
                for id in pmid_list:
                    o.write("%s\t%s\n" % (data[0], id))
            else:
                o.write("%s\tNone\n" % (data[0]))
