#!/usr/bin/env python

from Bio import Entrez


Entrez.email = "trestan.pillonel@unil.ch"

def string_id2pubmed_id_list(accession):

    import urllib2
    link = 'http://string-db.org/api/tsv/abstractsList?identifiers=%s' % accession
    print link
    try:
        data = urllib2.urlopen(link).read().rstrip().decode('utf-8').split('\n')[1:]
    except urllib2.URLError:
        print 'echec', link
        return False
    pid_list = [row.split(':')[1] for row in data]
    print 'list', pid_list
    return pid_list


def string_id2connexions(accession):

    '''

    global_score    : score
    neighborhood    : n score
    gene_fusion     : f score
    cooccurence     : p score
    coexpression    : a score
    experiments     : e score
    databases       : d score
    textmining      : t score

    :return:
    '''
    import urllib2

    link = 'http://string-db.org/api/psi-mi-tab/interactionsList?identifiers=%s' % accession

    try:
        data = urllib2.urlopen(link).read().decode('utf-8').split('\n')
    except urllib2.URLError:
        print 'echec', link
        return False
    interaction_list = []
    for row in data:
        #print row
        row_data = row.rstrip().split('\t')
        if len(row_data)>1:
            interaction_list.append([row_data[0],
                                     row_data[1],
                                     row_data[2],
                                     row_data[3],
                                     row_data[14]])
    return interaction_list



def biodb2all_connections(biodb):

    import manipulate_biosqldb
    import time
    import re

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select db_accession from custom_tables.uniprot_id2seqfeature_id_%s t0 ' \
          ' inner join custom_tables.uniprot_db_xref_%s t1 on t0.uniprot_id=t1.uniprot_id ' \
          ' inner join custom_tables.db_xref t2 on t1.db_xref_id=t2.db_xref_id where db_xref_name="string" and db_accession like "%%%%CPn%%%%";' % (biodb, biodb)

    all_string_accessions = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    sql = 'select seqfeature_id, taxon_id from custom_tables.locus2seqfeature_id_%s' % biodb

    seqfeature_id2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag,seqfeature_id from custom_tables.locus2seqfeature_id_%s' % biodb

    new_locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select old_locus_tag,seqfeature_id from custom_tables.seqfeature_id2old_locus_tag_%s' % biodb

    old_locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'create table if not exists string.interactions_%s (taxon_id INT, ' \
          ' seqfeature_id_1 INT, ' \
          ' seqfeature_id_2 INT,' \
          ' old_locus_tag_1 varchar(400), ' \
          ' old_locus_tag_2 varchar (400), ' \
          ' label_1 varchar(400), ' \
          ' label_2 varchar (400), ' \
          ' global_score FLOAT,' \
          ' neighborhood FLOAT,' \
          ' gene_fusion FLOAT,' \
          ' cooccurence FLOAT,' \
          ' coexpression FLOAT,' \
          ' experiments FLOAT,' \
          ' biodatabases FLOAT,' \
          ' textmining FLOAT, ' \
          ' index seqfeature_id_1 (seqfeature_id_1),' \
          ' index seqfeature_id_2 (seqfeature_id_2),' \
          ' INDEX old_locus_tag_1 (old_locus_tag_1),' \
          ' index old_locus_tag_2 (old_locus_tag_2))' % biodb
    print sql
    #server.adaptor.execute(sql,)

    ref_locus_list = []

    for n, string_accession in enumerate(all_string_accessions):
        print "%s / %s" % (n, len(all_string_accessions))
        interactions = string_id2connexions(string_accession)

        if not interactions:
            while interactions is False:
                print 'trying again...'
                time.sleep(10)
                interactions = string_id2connexions(string_accession)

        for one_interaction in interactions:
            print string_accession, one_interaction
            gscore = 0
            fscore = 0
            pscore = 0
            nscore = 0
            ascore = 0
            escore = 0
            dscore = 0
            tscore = 0

            if string_accession in one_interaction[0]:
                ref_locus = one_interaction[0].split(':')[1].split('.')[1]
                link_locus = one_interaction[1].split(':')[1].split('.')[1]

            elif string_accession in one_interaction[1]:
                ref_locus = one_interaction[1].split(':')[1].split('.')[1]
                link_locus = one_interaction[0].split(':')[1].split('.')[1]
            else:
                # connection does not contain reference link, skiping
                continue
            ref_locus_list.append(ref_locus)
            if link_locus in ref_locus_list:
                # not a new connection
                continue

            label_1 = one_interaction[2]
            label_2 = one_interaction[3]

            # locus tag corresp between old and new RefSeq annotation
            try:
                ref_locus_seqfeature_id = old_locus_tag2seqfeature_id[ref_locus]
            except:
                # special case trachomatis
                try:
                    ref_locus = re.sub('CT','CT_',ref_locus)
                    ref_locus_seqfeature_id = new_locus_tag2seqfeature_id[ref_locus]
                except:
                    ref_locus_seqfeature_id = 'NULL'
            print 'ref_locus', ref_locus
            # locus tag corresp OK but pseudogene
            try:
                taxon_id = seqfeature_id2taxon_id[str(ref_locus_seqfeature_id)]
            except:
                taxon_id = 'NULL'
            if taxon_id is None:
                taxon_id = 'NULL'
            # locus tag corresp between old and new RefSeq annotation
            try:
                link_locus_seqfeature_id = old_locus_tag2seqfeature_id[link_locus]
            except:
                try:
                    link_locus = re.sub('CT','CT_',link_locus)
                    link_locus_seqfeature_id = new_locus_tag2seqfeature_id[link_locus]
                except:
                    link_locus_seqfeature_id = 'NULL'

            scores = one_interaction[4].split('|')

            for one_score in scores:
                score, value = one_score.split(':')
                #print ref_locus, link_locus, score, value
                if score == 'score':
                    gscore = value
                elif score == 'nscore':
                    nscore = value
                elif score == 'fscore':
                    fscore = value
                elif score == 'pscore':
                    pscore = value
                elif score == 'ascore':
                    ascore = value
                elif score == 'escore':
                    escore = value
                elif score == 'dscore':
                    dscore = value
                elif score == 'tscore':
                    tscore = value
                else:
                    print 'unkonwn score type', score, value
            # ref_locus, link_locus, ref_locus_seqfeature_id, link_locus_seqfeature_id, label_1, label_2, gscore, ncore, fscore, pscore, ascore, escore, dscore, tscore
            sql = 'insert into string.interactions_%s values ' \
                  ' (%s, %s, %s, "%s", "%s", "%s", "%s", %s, %s, %s, %s, %s, %s, %s, %s)' % (biodb,
                                                                                                 taxon_id,
                                                                                                 ref_locus_seqfeature_id,
                                                                                                 link_locus_seqfeature_id,
                                                                                                 ref_locus,
                                                                                                 link_locus,
                                                                                                 label_1,
                                                                                                 label_2,
                                                                                                 gscore,
                                                                                                 nscore,
                                                                                                 fscore,
                                                                                                 pscore,
                                                                                                 ascore,
                                                                                                 escore,
                                                                                                 dscore,
                                                                                                 tscore)
            print taxon_id, sql
            server.adaptor.execute(sql,)
        server.commit()


def biodb2string_pmid_data(biodb):

    import manipulate_biosqldb
    import pubmed_utils
    import time
    import re

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'select db_accession from custom_tables.uniprot_id2seqfeature_id_%s t0 ' \
          ' inner join custom_tables.uniprot_db_xref_%s t1 on t0.uniprot_id=t1.uniprot_id ' \
          ' inner join custom_tables.db_xref t2 on t1.db_xref_id=t2.db_xref_id where db_xref_name="string" and db_accession like "%%%%CPn%%%%";' % (biodb, biodb)

    all_string_accessions = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    sql = 'create table if not exists string.seqfeature_id2string_pmid_%s (taxon_id INT, ' \
          ' seqfeature_id INT, ' \
          ' pmid INT, ' \
          ' authors TEXT,' \
          ' title TEXT,' \
          ' abstract TEXT, ' \
          ' source TEXT,' \
          ' INDEX seqfeature_id(seqfeature_id))' % biodb

    server.adaptor.execute(sql,)
    server.commit()

    sql = 'select seqfeature_id, taxon_id from custom_tables.locus2seqfeature_id_%s' % biodb

    seqfeature_id2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select old_locus_tag, seqfeature_id from custom_tables.seqfeature_id2old_locus_tag_%s' % biodb

    old_locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    sql = 'select locus_tag,seqfeature_id from custom_tables.locus2seqfeature_id_%s' % biodb

    new_locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    for n, string_accession in enumerate(all_string_accessions):
        print "%s / %s" % (n, len(all_string_accessions))

        old_locus_tag = string_accession.split('.')[1]
        try:
            seqfeature_id = old_locus_tag2seqfeature_id[old_locus_tag]
        except:
            try:
                # special case trachomatis
                old_locus_tag = re.sub('CT', 'CT_', old_locus_tag)
                seqfeature_id = new_locus_tag2seqfeature_id[old_locus_tag]
            except:
                continue
        taxon_id = seqfeature_id2taxon_id[str(seqfeature_id)]
        if taxon_id is None:
            taxon_id = 'NULL'
        pmid_list = string_id2pubmed_id_list(string_accession)
        print 'miidjdjnjdhd', pmid_list
        if pmid_list is False:
            while pmid_list is False:
                print 'trying again'
                time.sleep(10)
                pmid_list = string_id2pubmed_id_list(string_accession)


        if len(pmid_list) == 0:
            print '0 pmid for', string_accession
            continue
        else:
            for one_pmid in pmid_list:
                abstract_data = pubmed_utils.pmid2abstract_info(one_pmid)
                print 'data', abstract_data
                abstract = re.sub("'", "", abstract_data['abstract'])
                abstract = re.sub("%", "%%%%", abstract)
                title = re.sub("'", "", abstract_data['title'])
                title = re.sub("%", "%%%%", title)
                source = re.sub("'", "", abstract_data['source'])
                source = re.sub("%", "%%%%", source)

                sql = '''insert into string.seqfeature_id2string_pmid_%s values (%s, %s, %s, '%s', '%s', '%s', '%s')''' % (biodb,
                                                                                                                 taxon_id,
                                                                                                                 seqfeature_id,
                                                                                                                 abstract_data['pmid'],
                                                                                                                 re.sub("'", "", str(abstract_data['authors'])),
                                                                                                                 title,
                                                                                                                 abstract,
                                                                                                                 source)
                print sql
                server.adaptor.execute(sql,)
            server.commit()

