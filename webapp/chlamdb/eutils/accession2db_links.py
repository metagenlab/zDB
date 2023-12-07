#!/usr/bin/env python

from Bio import Entrez


Entrez.email = "trestan.pillonel@unil.ch"


def uniprot2uniref100_with_refseq_or_embl_crossref(uniprot_id):
    import urllib
    import urllib2

    # to=EMBL_ID&query=A4D962&from=ACC&format=tab
    url = 'http://www.uniprot.org/uploadlists/'

    params = {
        'from': 'ACC',
        'to': 'NF100',
        'format': 'tab',
        'query': '%s' % uniprot_id
    }

    data = urllib.urlencode(params)
    request = urllib2.Request(url, data)
    # Please set your email address here to help us debug in case of problems.
    contact = "trestan.pillonel@gmail.com"
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib2.urlopen(request)
    page = response.read()
    result = page.split('\n')

    cluster_members = result[1].split('\t')[5].split(';')
    for member in cluster_members:
        if member != uniprot_id:
            print 'member', member
            acc = uniprot_accession2EMBL_GenBank_DDBJ_id(
                member, 'ACC', 'P_REFSEQ_AC')
            if acc is not None:
                return acc
            else:
                acc = uniprot_accession2EMBL_GenBank_DDBJ_id(
                    member, 'ACC', 'EMBL')
                if acc is not None:
                    return acc

    return None
    # acc = result[1].split('\t')[1]
    # return acc


def uniparc2active_refseq(uniparc_id):
    import urllib
    import urllib2
    url = 'http://www.uniprot.org/uniparc/%s.tab' % uniparc_id
    request = urllib2.Request(url)
    # Please set your email address here to help us debug in case of problems.
    contact = "trestan.pillonel@gmail.com"
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib2.urlopen(request)
    page = response.read()
    result = page.split('\n')
    for row in result:
        data = row.split('\t')
        if data[0] == 'RefSeq':
            print 'jackpot!'
            return data[1]
    return None


def uniparc2active_uniprotkb(uniparc_id):
    import urllib
    import urllib2
    url = 'http://www.uniprot.org/uniparc/%s.tab' % uniparc_id
    request = urllib2.Request(url)
    # Please set your email address here to help us debug in case of problems.
    contact = "trestan.pillonel@gmail.com"
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib2.urlopen(request)
    page = response.read()
    result = page.split('\n')
    for row in result:
        data = row.split('\t')
        if data[0] == 'UniProtKB/TrEMBL' and data[-1] == 'Yes':
            print 'jackpot UNIP!'
            return data[1]
    return None


def uniprot_accession2EMBL_GenBank_DDBJ_id(accession, dbfrom='ACC', dbto='EMBL'):
    import urllib
    import urllib2

    # to=EMBL_ID&query=A4D962&from=ACC&format=tab
    url = 'http://www.uniprot.org/uploadlists/'

    params = {
        'from': '%s' % dbfrom,
        'to': '%s' % dbto,
        'format': 'tab',
        'query': '%s' % accession
    }

    data = urllib.urlencode(params)
    request = urllib2.Request(url, data)
    # Please set your email address here to help us debug in case of problems.
    contact = "trestan.pillonel@gmail.com"
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib2.urlopen(request)
    page = response.read()
    result = page.split('\n')
    try:
        acc = result[1].split('\t')[1]
        return acc
    except:
        return None


def sequence2uniprot_id(sequence):

    link = "http://research.bioinformatics.udel.edu/peptidematch/webservices/peptidematch_rest?query=%s" % (
        sequence)

    from StringIO import StringIO
    import urllib2
    try:
        data = urllib2.urlopen(link).read().decode('utf-8')
    except urllib2.URLError:
        print 'echec', link
        return False
    for row in data.split('\n'):
        if len(row) == 0:
            continue
        elif row[0] == '#':
            continue
        elif len(row.split('\t')) == 3:
            return False
        else:
            return row.split('\t')[0]
    # rec = StringIO(data)
    # print rec


def gi2accession(ncbi_id):

    # https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=663070995,568815587&rettype=acc
    handle = Entrez.efetch(db='protein', id=ncbi_id,
                           rettype="acc", retmode='text')
    rows = [i.rstrip() for i in handle]
    if len(rows) == 1:
        return rows[0]
    else:
        return False


def uniprot_id2record(uniprot_accession, n_trial=0):

    link = "http://www.uniprot.org/uniprot/%s.xml" % (uniprot_accession)

    from StringIO import StringIO
    import urllib2
    try:
        data = urllib2.urlopen(link).read().decode('utf-8')
    except urllib2.URLError:
        print 'echec', link
        return False
    rec = StringIO(data)

    try:
        record = SeqIO.read(rec, 'uniprot-xml')
    except:
        import time
        print 'problem with %s, trying again, %s th trial' % (uniprot_accession, n_trial)
        time.sleep(50)
        n_trial += 1
        if n_trial < 5:
            record = uniprot_id2record(uniprot_accession, n_trial=n_trial)
        else:
            return False
    return record


def uniprot_record2annotations(record):

    data2info = {}
    if 'recommendedName_ecNumber' in record.annotations:
        data2info['ec_number'] = record.annotations['recommendedName_ecNumber'][0]
    else:
        data2info['ec_number'] = '-'
    if 'gene_name_primary' in record.annotations:
        data2info['gene'] = record.annotations['gene_name_primary']
    else:
        data2info['gene'] = '-'
    if 'comment_catalyticactivity' in record.annotations:
        data2info['comment_catalyticactivity'] = record.annotations['comment_catalyticactivity'][0]
    else:
        data2info['comment_catalyticactivity'] = '-'
    if 'comment_subunit' in record.annotations:
        data2info['comment_subunit'] = record.annotations['comment_subunit'][0]
    else:
        data2info['comment_subunit'] = '-'
    if 'comment_function' in record.annotations:
        data2info['comment_function'] = record.annotations['comment_function'][0]
    else:
        data2info['comment_function'] = '-'
    if 'recommendedName_fullName' in record.annotations:
        data2info['recommendedName_fullName'] = record.annotations['recommendedName_fullName'][0]
    else:
        data2info['recommendedName_fullName'] = '-'
    if 'comment_similarity' in record.annotations:
        data2info['comment_similarity'] = record.annotations['comment_similarity'][0]
    else:
        data2info['comment_similarity'] = '-'
    if 'proteinExistence' in record.annotations:
        data2info['proteinExistence'] = record.annotations['proteinExistence'][0]
    else:
        data2info['proteinExistence'] = '-'
    if 'keywords' in record.annotations:
        data2info['keywords'] = ';'.join(record.annotations['keywords'])
    else:
        data2info['keywords'] = '-'
    if 'comment_pathway' in record.annotations:
        data2info['comment_pathway'] = record.annotations['comment_pathway'][0]
    else:
        data2info['comment_pathway'] = '-'
    if 'comment_developmentalstage' in record.annotations:
        data2info['developmentalstage'] = record.annotations['comment_developmentalstage'][0]
    else:
        data2info['developmentalstage'] = '-'

    return data2info


def uniprot_record2db_refs(record):

    ref_name2ref_id = {}
    for ref in record.dbxrefs:
        ref_dat = ref.split(':')
        if len(ref_dat) > 2:
            if ref_dat[0] not in ref_name2ref_id:
                ref_name2ref_id[ref_dat[0]] = [(':').join(ref_dat[1:])]
            else:
                ref_name2ref_id[ref_dat[0]].append((':').join(ref_dat[1:]))
        else:
            if ref_dat[0] not in ref_name2ref_id:
                ref_name2ref_id[ref_dat[0]] = [ref_dat[1]]
            else:
                ref_name2ref_id[ref_dat[0]].append(ref_dat[1])
    return ref_name2ref_id


def ncbi_accession2uniprotid(ncbi_accession, gene=False, organism=False, genome_accession=False):
    '''
    :param genbank/refseq protein accession
    :return: unirpot accession
    '''
    print genome_accession

    import urllib2
    from urllib2 import URLError
    if not gene:
        if genome_accession is False:
            link = "http://www.uniprot.org/uniprot/?query=%s&format=tab" % (
                ncbi_accession)
        else:

            link = "http://www.uniprot.org/uniprot/?query=%s+%s&format=tab" % (
                genome_accession, ncbi_accession)
        print link
    else:
        if not organism:
            link = "http://www.uniprot.org/uniprot/?query=gene:%s&format=tab" % (
                ncbi_accession)
        else:
            link = "http://www.uniprot.org/uniprot/?query=organism:%s+gene:%s&format=tab" % (
                organism, ncbi_accession)
    req = urllib2.Request(link)
    try:
        page = urllib2.urlopen(req)
        data = page.read().decode('utf-8').split('\n')
        rows = [i.rstrip().split('\t') for i in data]
    except URLError, err:
        print 'echec'
        print link
        return False
    try:
        return rows[1][0]
    except IndexError:
        return False


def uniprot_accession2go_and_status(uniprot_accession):

    import urllib2
    from urllib2 import URLError

    link = "http://www.uniprot.org/uniprot/?query=%s&columns=go,annotation score,reviewed&format=tab" % (
        uniprot_accession)
    link = link.replace(' ', '%20')

    try:
        req = urllib2.Request(link)
        page = urllib2.urlopen(req)
        data = page.read().decode('utf-8').split('\n')
        rows = [i.rstrip().split('\t') for i in data]
    except URLError:
        success = False
        while not success:
            import time
            print 'connection problem, trying again...'
            time.sleep(10)
            try:
                req = urllib2.Request(link)
                page = urllib2.urlopen(req)
                data = page.read().decode('utf-8').split('\n')
                rows = [i.rstrip().split('\t') for i in data]
                success = True
            except:
                success = False

    go_id2description = {}
    if rows[1][0] != '':
        for go in rows[1][0].split(';'):
            go_id = go.split('[')[1][0:-1]
            go_description = go.split('[')[0][0:-1].lstrip()
            go_id2description[go_id] = go_description
    else:
        go_id2description = False
    return (rows[1][1][0], rows[1][2], go_id2description)


def get_whole_db_uniprot_crossref(biodb):

    # get gi from all database locus

    import MySQLdb
    from datetime import datetime
    import httplib
    import time
    import manipulate_biosqldb
    import re
    import urllib2
    import os
    # sqlpsw = os.environ['SQLPSW']
    from tempfile import NamedTemporaryFile
    conn = MySQLdb.connect(host="localhost",  # your host, usually localhost
                                user="root",  # your username
                                passwd="estrella3",  # your password
                                db="custom_tables")  # name of the data base

    cursor = conn.cursor()

    sql1 = 'CREATE TABLE IF NOT EXISTS uniprot_id2seqfeature_id_%s (seqfeature_id INT UNIQUE, uniprot_id INT AUTO_INCREMENT,' \
           ' uniprot_accession varchar(400), uniprot_status varchar(400), annotation_score INT, insert_date varchar(300), INDEX uniprot_id(uniprot_id))' % biodb

    sql2 = 'CREATE TABLE IF NOT EXISTS db_xref (db_xref_id INT AUTO_INCREMENT, db_xref_name varchar(200) UNIQUE, INDEX db_xref_id(db_xref_id))'

    sql3 = 'CREATE TABLE IF NOT EXISTS uniprot_db_xref_%s (uniprot_id INT, db_xref_id INT, db_accession varchar(200), ' \
           ' INDEX db_xref_id(db_xref_id), index uniprot_id(uniprot_id))' % biodb

    sql4 = 'CREATE TABLE IF NOT EXISTS uniprot_go_terms_%s (seqfeature_id INT, go_term_id varchar(400), go_description TEXT, ' \
           ' INDEX seqfeature_id(seqfeature_id))' % biodb

    sql5 = 'CREATE TABLE IF NOT EXISTS uniprot_annotation_%s (seqfeature_id INT, comment_function TEXT,' \
           ' ec_number TEXT,comment_similarity TEXT,comment_catalyticactivity TEXT,comment_pathway TEXT,keywords TEXT,' \
           ' comment_subunit TEXT, gene TEXT, recommendedName_fullName TEXT, proteinExistence TEXT, ' \
           ' developmentalstage TEXT, index seqfeature_id(seqfeature_id))' % biodb

    print sql1
    cursor.execute(sql1, )
    cursor.execute(sql2, )
    cursor.execute(sql3, )
    cursor.execute(sql4, )
    cursor.execute(sql5, )
    conn.commit()

    sql1 = 'select locus_tag, seqfeature_id from locus2seqfeature_id_%s' % biodb
    sql2 = 'select locus_tag, old_locus_tag from biosqldb.locus_tag2old_locus_tag'

    # attention EDIT!!!!!!!!!!!!!!
    sql3 = 'select locus_tag, protein_id from biosqldb.orthology_detail_%s where protein_id not like "%%%%CHUV%%%%"' % biodb

    sql4 = 'select locus_tag,t2.seqfeature_id from locus2seqfeature_id_%s t1 inner join uniprot_annotation_%s t2' \
           ' on t1.seqfeature_id=t2.seqfeature_id group by locus_tag;' % (
               biodb, biodb)
    sql5 = 'select locus_tag, organism from biosqldb.orthology_detail_%s' % biodb
    sql6 = 'select locus_tag, translation from biosqldb.orthology_detail_%s' % biodb

    sql7 = 'select locus_tag, accession from biosqldb.orthology_detail_%s' % biodb

    cursor.execute(sql1, )
    locus2seqfeature_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    cursor.execute(sql2, )
    locus2old_locus = manipulate_biosqldb.to_dict(cursor.fetchall())

    cursor.execute(sql3, )
    locus2protein_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    cursor.execute(sql4, )
    locus2uniprot_id = manipulate_biosqldb.to_dict(cursor.fetchall())

    cursor.execute(sql5, )
    locus2organism = manipulate_biosqldb.to_dict(cursor.fetchall())

    cursor.execute(sql6, )
    locus2sequence = manipulate_biosqldb.to_dict(cursor.fetchall())

    cursor.execute(sql7, )
    locus2genome_accession = manipulate_biosqldb.to_dict(cursor.fetchall())

    for i, locus in enumerate(locus2protein_id):
        print "%s -- %s : %s / %s" % (locus, locus2protein_id[locus], i, len(locus2protein_id))

        # already into database
        if locus in locus2uniprot_id:
            continue
        genome_accession = locus2genome_accession[locus]

        uniprot_id = ncbi_accession2uniprotid(
            locus2protein_id[locus], genome_accession=genome_accession)

        if not uniprot_id:
            uniprot_id = ncbi_accession2uniprotid(locus2protein_id[locus])

        if not uniprot_id:
            try:
                old_locus = locus2old_locus[locus]
            except KeyError:
                old_locus = False
            if old_locus:
                genus = locus2organism[locus].split(' ')[0]
                print 'trying with old_locus_tag'
                uniprot_id = ncbi_accession2uniprotid(
                    old_locus, gene=True, organism=genus)
            if not uniprot_id:
                print 'trying to match with sequence: %s' % locus
                try:
                    uniprot_id = sequence2uniprot_id(locus2sequence[locus])
                    print 'ok: %s' % uniprot_id
                except:
                    continue
        if uniprot_id:
            # insert uniprot_id into mysql table
            # 1. get seqfeatureid of the corresponding locus
            seqid = locus2seqfeature_id[locus]
            try:
                uniprot_score, uniprot_status, go_data = uniprot_accession2go_and_status(
                    uniprot_id)
            except:
                print 'echec, continue'
                continue
            # add go data
            if go_data:
                for one_go in go_data:
                    sql = 'insert into uniprot_go_terms_%s (seqfeature_id, go_term_id, go_description) ' \
                          'values(%s, "%s", "%s")' % (biodb,
                                                      seqid,
                                                      one_go,
                                                      go_data[one_go])
                    cursor.execute(sql, )
                    conn.commit()

            # insert uniprot_id
            now = datetime.now()
            str_date = "%s-%s-%s" % (now.year, now.month, now.day)

            sql = 'insert into uniprot_id2seqfeature_id_%s (seqfeature_id,uniprot_accession, uniprot_status, annotation_score,insert_date) ' \
                  ' values (%s, "%s", "%s", %s,"%s")' % (biodb,
                                                         seqid,
                                                         uniprot_id,
                                                         uniprot_status,
                                                         uniprot_score,
                                                         str_date)
            try:
                cursor.execute(sql, )
                conn.commit()
            # if seqfeature id already already inserted, no need to insert it again
            except conn.IntegrityError:
                print '%s already into uniprot_id2seqfeature_id_%s' % (seqid, biodb)
                pass
            sqlid = 'select t1.uniprot_id from uniprot_id2seqfeature_id_%s as t1 where t1.seqfeature_id=%s' % (biodb,
                                                                                                               locus2seqfeature_id[locus])
            # print sqlid
            cursor.execute(sqlid, )
            uniprot_db_id = cursor.fetchall()[0][0]
            # print 'uniprotdb id', uniprot_db_id

            uniprot_record = uniprot_id2record(uniprot_id)

            if not uniprot_record:
                import time
                time.sleep(5)
                uniprot_record = uniprot_id2record(uniprot_id)

            alldbref = uniprot_record2db_refs(uniprot_record)

            annotation = uniprot_record2annotations(uniprot_record)

            # add annotation
            sql = 'insert into uniprot_annotation_%s (seqfeature_id, comment_function,' \
                  ' ec_number,comment_similarity,comment_catalyticactivity,comment_pathway,keywords,' \
                  ' comment_subunit, gene, recommendedName_fullName, proteinExistence,developmentalstage) values' \
                  '  (%s, "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s", "%s")' % (biodb,
                                                                                                seqid,
                                                                                                re.sub(
                                                                                                    '"', '', annotation["comment_function"]),
                                                                                                re.sub(
                                                                                                    '"', '', annotation["ec_number"]),
                                                                                                re.sub(
                                                                                                    '"', '', annotation["comment_similarity"]),
                                                                                                re.sub(
                                                                                                    '"', '', annotation["comment_catalyticactivity"]),
                                                                                                re.sub(
                                                                                                    '"', '', annotation["comment_pathway"]),
                                                                                                re.sub(
                                                                                                    '"', '', annotation["keywords"]),
                                                                                                re.sub(
                                                                                                    '"', '', annotation["comment_subunit"]),
                                                                                                re.sub(
                                                                                                    '"', '', annotation["gene"]),
                                                                                                re.sub(
                                                                                                    '"', '', annotation["recommendedName_fullName"]),
                                                                                                re.sub(
                                                                                                    '"', '', annotation["proteinExistence"]),
                                                                                                re.sub('"', '', annotation["developmentalstage"]))
            cursor.execute(sql, )
            conn.commit()

            # add dbxrefs
            if alldbref:
                for database in alldbref:
                    # 1. check if cross ref database already in the database list
                    sql1 = 'select db_xref_id from db_xref where db_xref_name="%s"' % database
                    try:
                        cursor.execute(sql1, )
                        database_index = cursor.fetchall()[0][0]

                    except:
                        # insert new database name
                        sql2 = 'insert into db_xref (db_xref_name) values ("%s")' % database
                        cursor.execute(sql2, )
                        conn.commit()
                        cursor.execute(sql1, )
                        database_index = cursor.fetchall()[0][0]

                    for crossref in alldbref[database]:
                        # insert cross reference into database
                        sql3 = 'insert into uniprot_db_xref_%s (uniprot_id, db_xref_id, db_accession) values (%s, %s, "%s")' % (biodb,
                                                                                                                                uniprot_db_id,
                                                                                                                                database_index,
                                                                                                                                crossref)
                        # print sql3
                        cursor.execute(sql3, )
                        conn.commit()
            else:
                print 'echec ----------------'
                print 'echec ----------------'
        else:
            print 'UNIPRITID NOT FOUND'


if __name__ == '__main__':
    import argparse
    import sys
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--biodb', type=str, help="biodb_name")

    args = parser.parse_args()

    # seq = 'BMNKKTETGLRESMNKMRTSLGNFKNFRGVFSNDIGIDLGTANTLVYVRGKGIVLAEPSVVAVDSNTNEVLAVGHKAKAMLGRTPRKIHAVRPMKDGVIADFEVAEGMLKALIKRVTPSRSLFRPKILIAVPSGITGVEKRAVEDSALHAGAQEVILIEEPMAAAIGVDLPVHEPSANMIIDIGGGTTEIAIISLGGIVESRSLRIAGDEFDECIINYMRRTYNLMIGPRTAEEIKMTIGSAYPLGENELEMEVRGRDQVAGLPVTKRINSVEIRECLSEPIQQIIESVKLTLEKCPPELAADLVERGMVVAGGGALIKGLDKYLIKETGLPVILAQNPLLAVCLGTGKALEYLDKFKKRKATA'
    # print sequence2uniprot_id(seq)
    get_whole_db_uniprot_crossref(args.biodb)
    # uaccession = ncbi_accession2uniprotid('pc1161', True, organism='protochlamydia')
    # print uaccession
    # accession2db_xrefs('15604941')

    # urecord = uniprot_id2record('O84081')
    # uniprot_record2annotations(urecord)
    '''
    uaccession = ncbi_accession2uniprotid('BN1013_00847')
    urecord = uniprot_id2record('O84395')
    print uaccession
    annotation = uniprot_record2annotations(urecord)
    print annotation
    '''
    # print uniprot_accession2go_and_status(uaccession)
