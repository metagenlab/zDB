#!/usr/bin/env python

from Bio import Entrez


Entrez.email = "trestan.pillonel@unil.ch"


def GOLD_id2data(gold_id):
    import urllib2
    from bs4 import BeautifulSoup
    import re
    from urllib2 import URLError

    link = "https://gold.jgi.doe.gov/projects?id=%s" % (gold_id)
    req = urllib2.Request(link)
    try:
        page = urllib2.urlopen(req)
        # data = page.read().decode('utf-8').split('\n')
        soup = BeautifulSoup(page, 'html.parser')
    except:
        import time
        print 'connexion problem, trying again...'
        time.sleep(60)
        return GOLD_id2data(gold_id)
    project_description = soup.find(
        "div", {"id": "id__22__type__String__column__"}).contents[0].rstrip()
    print soup.find("div", {"id": "id__25__type__List__column__"}).contents


if __name__ == '__main__':
    import argparse
    import sys
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--seq_proj_GOLD', default=False,
                        type=str, help="GOLD sequencing project ID")

    args = parser.parse_args()

    GOLD_id2data(args.seq_proj_GOLD)
