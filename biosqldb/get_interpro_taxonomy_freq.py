#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-



def interpro_entry2taxonomy(accession):

    import urllib2
    from bs4 import BeautifulSoup
    import re


    c_bact = 0
    c_euk = 0
    c_arch = 0
    c_virus = 0

    '''
    link_euk = 'https://www.ebi.ac.uk/interpro/entry/%s/proteins-matched?taxonomy=2759&export=ids' % accession
    link_bact = 'https://www.ebi.ac.uk/interpro/entry/%s/proteins-matched?taxonomy=2&export=ids' % accession
    link_archae = 'https://www.ebi.ac.uk/interpro/entry/%s/proteins-matched?taxonomy=2157&export=ids' % accession

    req_euk = urllib2.Request(link_euk)
    req_bact = urllib2.Request(link_bact)
    req_archae = urllib2.Request(link_archae)

    response_euk = urllib2.urlopen(req_euk)
    response_bact = urllib2.urlopen(req_bact)
    response_archae = urllib2.urlopen(req_archae)

    c_euk = len([i for i in response_euk])
    c_bact = len([i for i in response_bact])
    c_archae = len([i for i in response_archae])
    '''

    link_entry = 'https://www.ebi.ac.uk/interpro/entry/%s/taxonomy' % accession
    req = urllib2.Request(link_entry)

    soup = BeautifulSoup(urllib2.urlopen(req), 'html.parser')
    data = soup.find_all('a')
    m = re.compile('open_([0-9]+)')

    m_cellular_org = re.compile('.*cellular organisms.*', re.MULTILINE)
    m_viruses = re.compile('.*Viruses.*', re.MULTILINE)

    counts_id_viruses = False
    counts_id_cellular_org = False

    for i in data:
        i = re.sub('\s+',' ', str(i))
        s = re.search(m, str(i))

        if s and re.match(m_viruses, str(i)) is not None:
            counts_id_viruses = s.group(1)
        if s and re.match(m_cellular_org, str(i)) is not None:
            counts_id_cellular_org = s.group(1)

    if counts_id_cellular_org:

        link_taxcounts = 'https://www.ebi.ac.uk/interpro/entry/childTax/%s/' % counts_id_cellular_org
        req_taxcounts = urllib2.Request(link_taxcounts)
        response_taxcounts = urllib2.urlopen(req_taxcounts)
        soup = BeautifulSoup(response_taxcounts, 'html.parser')
        data = soup.find_all('a')

        bact = re.compile('taxonomy=2">([0-9]+)')
        euk = re.compile('taxonomy=2759">([0-9]+)')
        arch = re.compile('taxonomy=2157">([0-9]+)')

        for i in data:
            s_bact = re.search(bact,str(i))
            s_euk = re.search(euk,str(i))
            s_arch = re.search(arch,str(i))

            if s_bact:
                c_bact = int(s_bact.group(1))
            if s_euk:
                c_euk = int(s_euk.group(1))
            if s_arch:
                c_arch = int(s_arch.group(1))
    if counts_id_viruses:
        link_taxcounts = 'https://www.ebi.ac.uk/interpro/entry/childTax/%s/' % counts_id_viruses
        req_taxcounts = urllib2.Request(link_taxcounts)
        response_taxcounts = urllib2.urlopen(req_taxcounts)
        soup = BeautifulSoup(response_taxcounts, 'html.parser')
        data = soup.find_all('a')

        # match any kind of viruses
        virus = re.compile('taxonomy=[0-9]+">([0-9]+)')

        for i in data:
            s_virus = re.search(virus, str(i))
            if s_virus:
                c_virus += int(s_virus.group(1))


    return (c_euk, c_bact, c_arch, c_virus)






def get_interpro_version():
    link = 'https://www.ebi.ac.uk/interpro/release_notes.html'
    import urllib2

    req = urllib2.Request(link)
    response = urllib2.urlopen(req)

    for line in response:
        if 'version_title_main' in line:
            return float(line.split('>')[1].split('<')[0].split(' ')[1])

def create_interpro_taxonomy_table(interpro_version):

    import MySQLdb

    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd="estrella3", # your password
                                db="interpro") # name of the data base
    cursor = conn.cursor()


    sql = 'CREATE table IF NOT EXISTS interpro_taxonomy_v_%s (interpro_id INT, eukaryote INT, bacteria INT, ' \
          ' archae INT, virus INT, total INT, p_eukaryote FLOAT, p_bacteria FLOAT, p_archae FLOAT, p_virus FLOAT, index interpro_id(interpro_id))' % int(interpro_version)
    print sql
    cursor.execute(sql,)
    conn.commit()

def get_whole_db_interpro_taxonomy():
    import MySQLdb
    import manipulate_biosqldb

    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd="estrella3", # your password
                                db="interpro") # name of the data base
    cursor = conn.cursor()


    interpro_vesrion = get_interpro_version()

    create_interpro_taxonomy_table(interpro_vesrion)

    sql = 'select name,interpro_id from interpro.entry'
    cursor.execute(sql,)
    interpro_accession2interpro_id = manipulate_biosqldb.to_dict(cursor.fetchall())
    for i, interpro_accession in enumerate(interpro_accession2interpro_id):
        print '%s / %s -- %s' % (i, len(interpro_accession2interpro_id), interpro_accession)
        euk, bact, arch, virus = interpro_entry2taxonomy(interpro_accession)

        total = float(euk + bact + arch + virus)
        sql = 'insert into interpro_taxonomy_v_%s values (%s,%s,%s,%s,%s,%s,%s,%s, %s, %s)' % (int(interpro_vesrion),
                                                                                    interpro_accession2interpro_id[interpro_accession],
                                                                                    euk,
                                                                                    bact,
                                                                                    arch,
                                                                                    virus,
                                                                                    total,
                                                                                    round((euk/total)*100,2),
                                                                                    round((bact/total)*100, 2),
                                                                                    round((arch/total)*100, 2),
                                                                                    round((virus/total)*100, 2))

        cursor.execute(sql,)
        conn.commit()

#print interpro_entry2taxonomy("IPR014000")
get_whole_db_interpro_taxonomy()
