#!/usr/bin/env python

from Bio import Entrez


Entrez.email = "trestan.pillonel@unil.ch"

def accession2family(family_id):
    import urllib2
    from bs4 import BeautifulSoup
    import re
    from urllib2 import URLError

    link = "http://www.tcdb.org/search/result.php?tc=%s" % (family_id)
    req = urllib2.Request(link)
    try:
        page = urllib2.urlopen(req)
        #data = page.read().decode('utf-8').split('\n')
        soup = BeautifulSoup(page, 'html.parser')
    except:
        import time
        print 'connexion problem, trying again...'
        time.sleep(60)
    if len(family_id.split('.'))<=2:
        family_title = soup.find("div", {"class": "title"})
        return ' '.join(family_title.text.split(' ')[1:]).lstrip()
    elif len(family_id.split('.'))==4:
        family_title = soup.find_all("td", {"id": "right-border"} )
        for one_title in family_title:
            if family_id in one_title.text:
                return ' '.join(one_title.text.split(' ')[1:]).lstrip()
        return '-'
    elif len(family_id.split('.'))==5:
        match_table = soup.find("table")
        match_rows = match_table.findAll("tr")
        for one_row in match_rows:
            if family_id in one_row.text:
                return ' '.join(one_row.text.split(family_id)[1:]).lstrip()
        return '-'

    else:

        project_description = soup.find("div", {"class": "description"})
        return ' '.join(project_description.p.text.split(' ')[1:]).lstrip()

    '''
    if not project_description.span:
        print 'span'
        return project_description.p.text
    else:
        return project_description.p.text
        if project_description.span.text == ' ':
            return project_description.p.text
        else:
            print len(project_description.span.text)
            return project_description.span.text
    '''
def accession2substrate(accession, tc):
    import urllib2
    from bs4 import BeautifulSoup
    import re
    from urllib2 import URLError

    link = "http://www.tcdb.org/search/result.php?acc=%s&tc=%s" % (accession, tc)
    #print link
    req = urllib2.Request(link)
    try:
        page = urllib2.urlopen(req)
        #data = page.read().decode('utf-8').split('\n')
        soup = BeautifulSoup(page, 'html.parser')
    except:
        import time
        print 'connexion problem, trying again...'
        time.sleep(60)
    table = soup.find('table', attrs={'class':'proteins'})


    rows = table.find_all('tr')
    for row in rows:
        cols = row.find_all('td')
        cols = [ele.text.strip() for ele in cols]
        ths = row.find_all('th')
        th_text = [ele.text.strip() for ele in ths]
        if th_text[0] == 'Substrate':
            return cols[0]


def accession2species_and_product(accession, tc):
    import urllib2
    from bs4 import BeautifulSoup
    import re
    from urllib2 import URLError

    link = "http://www.tcdb.org/search/result.php?acc=%s&tc=%s" % (accession, tc)
    #print link
    req = urllib2.Request(link)
    try:
        page = urllib2.urlopen(req)
        #data = page.read().decode('utf-8').split('\n')
        soup = BeautifulSoup(page, 'html.parser')
    except:
        import time
        print 'connexion problem, trying again...'
        time.sleep(60)
    table = soup.find('table', attrs={'class':'proteins'})

    rows = table.find_all('tr')
    sp, nme = '',''
    for row in rows:
        cols = row.find_all('td')
        cols = [ele.text.strip() for ele in cols]
        ths = row.find_all('th')
        th_text = [ele.text.strip() for ele in ths]
        if th_text[0] == 'Species:':
            sp = cols[0]
        elif th_text[0] == 'Protein Name:':
            nme = cols[0]
        if sp and nme:
            return [sp, nme]
    if sp == '':
        rows = soup.findAll('tr')
        for row in rows:
            cols = row.find_all('td')
            cols = [ele.text.strip() for ele in cols]
            ths = row.find_all('th')
            th_text = [ele.text.strip() for ele in ths]
            if th_text[0] == 'Species:':
                sp = cols[0]
            elif th_text[0] == 'Protein Name:':
                nme = cols[0]
            if sp and nme:
                return [sp, nme]

if __name__ == '__main__':
    import argparse
    import sys
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",'--seq_proj_GOLD', default=False, type=str, help="GOLD sequencing project ID")

    args = parser.parse_args()


    #print accession2species_and_product("Q9Z7N9","2.A.1.4.6")
    #print accession2family("2.A.1.4.6")
    #print accession2family("2.A.2")
    #print accession2family("2.B.2")
    #print accession2family("2.A")
    #print accession2family("2")

    print accession2species_and_product("P71997","2.A.53.4.2")

    '''
    print "%s\t%s" % (accession2family("2.A.1"), accession2substrate("Q9Z7N9","2.A.1.4.6"))
    print "%s\t%s" % (accession2family("3.A.6"), accession2substrate("O84092","3.A.6.3.1"))
    print "%s\t%s" % (accession2family("2.A.12"), accession2substrate("Q6MEN4","2.A.12.1.7"))
    print "%s\t%s" % (accession2family("2.A.1"), accession2substrate("I7I1R6","2.A.1.53.14"))
    '''