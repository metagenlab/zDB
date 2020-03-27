#! /usr/bin/env python


def get_pfam_freq(biodb, 
                  db_version):

    import MySQLdb
    import os
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor
    
    sql_pfam_list = 'select hmm_id from pfam_summary_version_%s' % (db_version)
    cursor.execute(sql_pfam_list,)
    pfam_list = [i[0] for i in cursor.fetchall()]

    # get overall taxonomy statistics

    sql = 'select superkingdom, count(*) from refseq_ref_repres_genomes t1 ' \
          ' inner join blastnr_blastnr_taxonomy t2 on t1.taxid=t2.taxon_id group by superkingdom;'
    cursor.execute(sql, )
    superkingdom2count = manipulate_biosqldb.to_dict(cursor.fetchall())

    # create table
    # Archaea
    # Bacteria
    # Eukaryota
    # Viruses

    sql = 'create table pfam2superkingdom_frequency_%s (pfam_id INT primary key,' \
          ' archaea_freq FLOAT,' \
          ' bacteria_freq FLOAT,' \
          ' eukaryota_freq FLOAT,' \
          ' viruses_freq FLOAT,' \
          ' archaea_count INT,' \
          ' bacteria_count INT,' \
          ' eukaryota_count INT,' \
          ' viruses_count INT)' % db_version
    cursor.execute(sql, )
    conn.commit()

    # get domains taxonomy statistics
    for n, pfam in enumerate(pfam_list):
        print ('%s / %s' % (n, len(pfam_list)))
        sql = 'select B.superkingdom, count(*) from ' \
              ' (select t1.assembly_id, taxid from refseq_ref_repres_genomes_domains_pfam_%s t1 ' \
              ' inner join refseq_ref_repres_genomes t2 on t1.assembly_id=t2.assembly_id ' \
              ' where pfam_id=%s group by t1.assembly_id) A ' \
              ' inner join blastnr_blastnr_taxonomy B on A.taxid=B.taxon_id group by B.superkingdom;' % (db_version,
                                                                                                         pfam)
        cursor.execute(sql, )
        superkingdom2count_domain = manipulate_biosqldb.to_dict(cursor.fetchall())

        if 'Archaea' in superkingdom2count_domain:
            archaea_freq = int(superkingdom2count_domain['Archaea'])/int(superkingdom2count['Archaea'])
            archaea_count = int(superkingdom2count['Archaea'])
        else:
            archaea_freq = 0
            archaea_count = 0

        if 'Bacteria' in superkingdom2count_domain:
            bacteria_freq = int(superkingdom2count_domain['Bacteria'])/int(superkingdom2count['Bacteria'])
            bacteria_count = int(superkingdom2count['Bacteria'])

        else:
            bacteria_freq = 0
            bacteria_count = 0
        if 'Eukaryota' in superkingdom2count_domain:
            eukaryota_freq = int(superkingdom2count_domain['Eukaryota'])/int(superkingdom2count['Eukaryota'])
            eukaryota_count = int(superkingdom2count['Eukaryota'])

        else:
            eukaryota_freq = 0
            eukaryota_count = 0
        if 'Viruses' in superkingdom2count_domain:
            virus_freq = int(superkingdom2count_domain['Viruses'])/int(superkingdom2count['Viruses'])
            virus_count = int(superkingdom2count['Viruses'])
        else:
            virus_freq = 0
            virus_count = 0

        sql = 'insert into pfam2superkingdom_frequency_%s values (%s, %s, %s, %s, %s, %s, %s, %s, %s)' % (db_version,
                                                                                      pfam,
                                                                                      archaea_freq,
                                                                                      bacteria_freq,
                                                                                      eukaryota_freq,
                                                                                      virus_freq,
                                                                                      archaea_count,
                                                                                      bacteria_count,
                                                                                      eukaryota_count,
                                                                                      virus_count)

        cursor.execute(sql, )
    conn.commit()

def create_pfam_interpro_signature2pfam_id(biodb, pfam_db_id):

    '''

    TODO avoir une table complete des pfam signature presentes dans interpro

    :param biodb:
    :param pfam_db_id:
    :return:
    '''

    import MySQLdb
    import os
    from chlamdb.biosqldb import manipulate_biosqldb
    sqlpsw = os.environ['SQLPSW']

    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="interpro") # name of the data base
    cursor = conn.cursor()

    # create table
    sql = 'create table interpro_signature2pfam_id_%s (signature_id INT, pfam_id INT,' \
          ' index signature_id(signature_id),' \
          ' index pfam_hmm_id(pfam_id))' % biodb
    cursor.execute(sql,)
    conn.commit()

    # get cross referencing
    sql = 'select signature_id, signature_accession from analysis t1 inner join signature t2 on t1.analysis_id=t2.analysis_id where analysis_name="Pfam";'
    cursor.execute(sql,)
    interpro_signature2pfam_accession = manipulate_biosqldb.to_dict(cursor.fetchall())
    for interpro_signature in list(interpro_signature2pfam_accession.keys()):
        pfam_accession = interpro_signature2pfam_accession[interpro_signature]
        pf_id_sql = 'select hmm_id from pfam.pfam_summary_version_%s where hmm_accession like "%s%%%%"; ' % (pfam_db_id,
                                                                                                        pfam_accession)
        cursor.execute(pf_id_sql, )

        try:
            hmm_id = cursor.fetchall()[0][0]
            sql = 'insert into interpro_signature2pfam_id_%s values (%s, %s)' % (biodb,
                                                                                 interpro_signature,
                                                                                 hmm_id)
            cursor.execute(sql,)
        except IndexError:
            print ('echec!')
            print (pf_id_sql)
    conn.commit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", '--pfam_db_version', type=str, help="pfam db version number", required=True)
    parser.add_argument("-n", '--biodb', type=str, help="biodb name", required=True)

    args = parser.parse_args()

    create_pfam_interpro_signature2pfam_id('2017_06_29b_motile_chlamydiae', args.pfam_db_version)
    #get_pfam_freq(args.biodb,
    #               args.pfam_db_version)