#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-


# convert gbk file to faa
# manual writing with headers of the form
# gi|83716028|ref|YP_443839.1| matrix protein [Avian metapneumovirus]
# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2014
# ---------------------------------------------------------------------------

# '6 qacc sacc evalue nident pident positive gaps length qstart qend qcovs sstart send qseqid qgi qaccver '
# 25.08.16 changed identity cutoff from 30 to 25
# 01.09.16 changed identity cutoff from 25 to 30
# 21.12.16 changed identity cutoff from 25 to 20 and coverage from 60 to 50
def blast2COG(blast_file, coverage_cutoff=50, identity_cutoff=20):
    with open(blast_file, "r") as f:
        locus2hit_gi = {}
        locus2query_start = {}
        locus2query_end = {}
        locus2identity = {}
        locus2qcovs = {}
        locus2evalue = {}
        #TODO: add bitscore
        # 0 qgi
        # 1 qacc ok
        # 2 sgi
        # 3 sacc ok
        # 4 sscinames
        # 5 sskingdoms
        # 6 staxids
        # 7 evalue ok
        # 8 nident
        # 9 pident ok
        # 10 positive
        # 11 gaps
        # 12 length
        # 13 qstart ok
        # 14 qend ok
        # 15 qcovs ok
        # 16 sstart
        # 17 send
        # 18 sstrand
        # 19 stitle

        for line in f:
            data = line.rstrip().split('\t')
            try:
                locus_tag = data[1].split('|')[3]
            except IndexError:
                locus_tag = data[1]
                print (locus_tag)
            identity = float(data[9])
            query_coverage = float(data[15])

            hit_gi = data[3].split('|')[1]
            query_start = int(data[13])
            query_end = int(data[14])
            evalue = float(data[7])

            if identity >= identity_cutoff and query_coverage >= coverage_cutoff:
                locus2hit_gi[locus_tag] = hit_gi
                locus2query_start[locus_tag] = query_start
                locus2query_end[locus_tag] = query_end
                locus2identity[locus_tag] = identity
                locus2qcovs[locus_tag] = query_coverage
                locus2evalue[locus_tag] = evalue
        return locus2hit_gi, locus2query_start, locus2query_end, locus2identity, locus2qcovs, locus2evalue


def gi2COG(*protein_gi):
    from chlamdb.biosqldb import manipulate_biosqldb
    import os

    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn




    if len(protein_gi)>1:
        prot_id_filter='where protein_id in (%s' % protein_gi[0]
        for i in range(1,len(protein_gi)):
            prot_id_filter+=',%s' % protein_gi[i]
        prot_id_filter+=')'
    else:
        prot_id_filter = 'where protein_id="%s"' % protein_gi[0]

    sql ='select cog_2014.protein_id,cog_names_2014.COG_id,cog_names_2014.function, cog_names_2014.name ' \
         ' from cog_2014 inner join cog_names_2014 on cog_2014.COG_id=cog_names_2014.COG_id %s' % prot_id_filter

    cursor.execute(sql)

    return cursor.fetchall()

def load_locus2cog_into_sqldb(input_blast_files, biodb):
    import MySQLdb
    import os
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor
    #TODO: add bitscore
    # 0 qgi
    # 1 qacc ok
    # 2 sgi
    # 3 sacc ok
    # 4 sscinames
    # 5 sskingdoms
    # 6 staxids
    # 7 evalue ok
    # 8 nident
    # 9 pident ok
    # 10 positive
    # 11 gaps
    # 12 length
    # 13 qstart ok
    # 14 qend ok
    # 15 qcovs ok
    # 16 sstart
    # 17 send
    # 18 sstrand
    # 19 stitle

    # locus_tag2gi_hit_
    sql = 'create table COG_seqfeature_id2best_COG_hit (bioentry_id INT, ' \
          ' seqfeature_id INT, ' \
          ' hit_cog_id INT,' \
          ' hit_protein_id INT, ' \
          ' query_start int,' \
          ' query_end int,' \
          ' identity FLOAT,' \
          ' query_cov FLOAT,' \
          ' evalue FLOAT,' \
          ' index seqfeature_id (seqfeature_id), ' \
          ' index bioentry_id (bioentry_id),' \
          ' index hit_protein_id(hit_protein_id),' \
          ' index hit_cog_id(hit_cog_id))' % biodb

    cursor.execute(sql)
    conn.commit()
    sql = 'select locus_tag,bioentry_id from orthology_detail t1 ' \
          ' inner join biosqldb.bioentry as t2 on t1.accession=t2.accession ' \
          ' inner join biosqldb.biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
          ' where t3.name="%s"' % (biodb, biodb)
    sql2 = 'select protein_id, locus_tag from orthology_detail' % biodb
    sql3 = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % biodb
    sql4 = 'select protein_id,COG_id from COG.cog_2014'
    server, db = manipulate_biosqldb.load_db(biodb)
    locus2bioentry_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql))
    protein_id2locus_tag = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2))
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3))
    protein_id2COG_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql4))

    for input_blast in input_blast_files:
        print ('file', input_blast)
        locus2hit_gi, \
        locus2query_start, \
        locus2query_end, \
        locus2identity, \
        locus2qcovs, \
        locus2evalue = blast2COG(input_blast)

        for locus in locus2hit_gi:
            print ('locus', locus)
            sql = 'INSERT into seqfeature_id2best_COG_hit_%s (bioentry_id, ' \
                  ' seqfeature_id, ' \
                  ' hit_COG_id,' \
                  ' hit_protein_id,'\
                  ' query_start,' \
                  ' query_end,' \
                  ' identity,' \
                  ' query_cov,' \
                  ' evalue) VALUES (%s,%s, %s, %s, %s,%s, %s, %s, %s)' % (biodb,
                                                              locus2bioentry_id[locus],
                                                              locus_tag2seqfeature_id[locus],
                                                              protein_id2COG_id[str(locus2hit_gi[locus])],
                                                              locus2hit_gi[locus],
                                                              locus2query_start[locus],
                                                              locus2query_end[locus],
                                                              locus2identity[locus],
                                                              locus2qcovs[locus],
                                                              locus2evalue[locus])

            cursor.execute(sql)
    conn.commit()



def locus2function(input_blast_files, display_print=False,):

    import MySQLdb
    import os
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    locus2function_dico = {}
    for input_blast in input_blast_files:
        print ('file', input_blast)
        locus2hit_accession = blast2COG(input_blast)

        cogs = gi2COG(*locus2hit_accession.values())

        gi2cog_data = {}
        for cog in cogs:
            gi2cog_data[int(cog[0])] = cog[1:]

        for locus in locus2hit_accession:
            try:
                function = list(gi2cog_data[int(locus2hit_accession[locus])][1])
                locus2function_dico[locus] = function
            except:
                print ('problem with locus %s, skipping...' % locus)

    if display_print:
        for locus in locus2function_dico:
            for function in locus2function_dico[locus]:
                print ("%s\t%s" % (locus, function))
    else:
        return locus2function_dico




def investiguate_core_COGs(db_name, locus2function):
    from chlamdb.biosqldb import manipulate_biosqldb
    import biosql_own_sql_tables
    server, db = manipulate_biosqldb.load_db(db_name)
    sql = 'select taxon_id from bioentry' \
      ' inner join biodatabase on bioentry.biodatabase_id = biodatabase.biodatabase_id' \
      ' and biodatabase.name = "%s" group by taxon_id' % db_name

    all_accessions = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]

    all_accessions.pop(all_accessions.index(292))
    all_accessions.pop(all_accessions.index(125))
    all_accessions.pop(all_accessions.index(291))

    sql_include = ''

    for i in range(0, len(all_accessions)-1):
        sql_include += ' `%s` > 0 and ' % all_accessions[i]
    sql_include+='`%s` > 0' % all_accessions[-1]

    sql ='select orthogroup from orthology_%s where %s ' % (db_name, sql_include)

    core_groups = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
    #print len(core_groups)
    locus_tag2orthogroup = biosql_own_sql_tables.locus_tag2orthogroup(db_name)

    orthogroup2proteins = biosql_own_sql_tables.orthogroup2protein_id_list(db_name)

    COGs = ''
    core_cogs = {}
    for group in core_groups:
        proteins = orthogroup2proteins[group]
        core_cogs[group] = []
        for protein in proteins:
            try:
                COGs+= '%s (%s)\t' % (locus2function[protein], protein)
                core_cogs[group] += locus2function[protein]
            except:
                COGs+= '- (%s)\t' % protein
        COGs = COGs[0:-1] + '\n'
    #print COGs
    core_annot = 0
    core_non_annot = 0
    for i in core_cogs:
        if len(core_cogs[i])>0:
            core_annot+=1
        else:
            core_non_annot+=1
    print ('annot:', core_annot)
    print ('non annot', core_non_annot)

    with open('core_groups2cogs.tab', 'w') as f:

        for group in core_cogs:
            line = '%s\t' % group
            for one_cog in core_cogs[group]:
                line+='%s\t' % one_cog
            f.write(line[0:-1] + '\n')




    '''
    for locus in locus2function:
        group = locus_tag2orthogroup[locus]
        if group in core_groups:
            print group, locus, locus2function[locus]
    '''





if __name__ == '__main__':
    import argparse
    from Bio import SeqIO
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--input_blast', type=str, help="blast tab file", nargs='+')
    parser.add_argument("-o", '--outname', type=str, help="putput_name", default=False)
    parser.add_argument("-d", '--database_name', type=str, help="database name", default=False)


    args = parser.parse_args()

    load_locus2cog_into_sqldb(args.input_blast, args.database_name)
    '''
    locus2function_dico = locus2function(args.input_blast, display_print=False)
    print 'locus2f', len(locus2function_dico)
    investiguate_core_COGs('chlamydia_03_15', locus2function_dico)
    '''
