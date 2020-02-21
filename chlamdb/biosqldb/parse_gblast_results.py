#! /usr/bin/env python




def add_one_tc_id(biodb, tc_id):
    from chlamdb.biosqldb import manipulate_biosqldb
    import tcdb_utils
    import MySQLdb
    import re
    import os
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="transporters") # name of the data base
    cursor = conn.cursor()

    conn.set_character_set('utf8')
    cursor.execute('SET NAMES utf8;')
    cursor.execute('SET CHARACTER SET utf8;')
    cursor.execute('SET character_set_connection=utf8;')

    tc_db_id = check_if_tc_accession_already_in_db(biodb, tc_id)
    if not tc_db_id:
        print '%s not into database!' % tc_id
        description = re.sub('"', '', tcdb_utils.accession2family(tc_id))
        sql = 'insert into transporters.tc_table(tc_name, description) values ("%s", "%s")' % (tc_id, description)
        print sql
        cursor.execute(sql,)
        conn.commit()
        sql3 = 'select tc_id from transporters.tc_table where tc_name="%s"' % tc_id
        cursor.execute(sql3,)
        return cursor.fetchall()[0][0]
    else:
        return tc_db_id

def insert_new_tc_path(biodb, tc_name):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)



    separated_ids = tc_name.split('.')
    if len(separated_ids) != 5:
        print 'invalid id!'
        import sys
        sys.exit()
    else:
        id1 = separated_ids[0]
        id_superfamily = '.'.join(separated_ids[0:2])
        id_family = '.'.join(separated_ids[0:3])
        id_subfamily = '.'.join(separated_ids[0:4])
        id_complete = tc_name

        id_id1 = add_one_tc_id(biodb, id1)
        id_superfamily = add_one_tc_id(biodb, id_superfamily)
        id_family = add_one_tc_id(biodb, id_family)
        id_subfamily = add_one_tc_id(biodb, id_subfamily)
        id_complete = add_one_tc_id(biodb, id_complete)
    sql = 'select transporter_id from transporters.transporter_table where transporter_id=%s' % id_complete
    try:
        id = server.adaptor.execute_and_fetchall(sql,)[0][0]
        return id
    except:
        sql = 'insert into transporters.transporter_table (transporter_id, tc1,superfamily, family, subfamily) values ("%s",' \
              ' "%s","%s","%s","%s")' % (id_complete, id_id1, id_superfamily, id_family, id_subfamily)
        server.adaptor.execute(sql,)
        server.commit()
        return id_complete

def check_if_tc_accession_already_in_db(biodb, tc_name):
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'CREATE table if not exists transporters.transporter_table (transporter_id INT primary key, tc1 INT, superfamily INT,' \
          ' family INT, subfamily INT, index tc1 (tc1), index superfamily (superfamily),' \
          ' index family (family), index subfamily (subfamily))'
    print sql
    server.adaptor.execute(sql,)
    server.commit()

    sql2 = 'CREATE table if not EXISTS transporters.tc_table (tc_id INT primary key AUTO_INCREMENT, tc_name varchar(400), description TEXT)'
    server.adaptor.execute(sql2,)
    server.commit()

    sql3 = 'select tc_id from transporters.tc_table where tc_name="%s"' % tc_name
    try:
        tc_id = server.adaptor.execute_and_fetchall(sql3,)[0][0]
        return tc_id
    except:
        return False

def insert_uniprot_in_db(biodb, uniprot_accession, tc_id):
    from chlamdb.biosqldb import manipulate_biosqldb
    from chlamdb.biosqldb import blastswiss2sqltable
    import tcdb_utils
    import accession2taxon_id
    import MySQLdb
    import re
    import accession2taxon_id
    import sequence_id2scientific_classification
    import os
    sqlpsw = os.environ['SQLPSW']
    conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                user="root", # your username
                                passwd=sqlpsw, # your password
                                db="transporters") # name of the data base
    cursor = conn.cursor()

    conn.set_character_set('utf8')
    cursor.execute('SET NAMES utf8;')
    cursor.execute('SET CHARACTER SET utf8;')
    cursor.execute('SET character_set_connection=utf8;')

    sql = 'CREATE TABLE if not exists transporters.uniprot_table (uniprot_id INT primary key AUTO_INCREMENT, uniprot_accession varchar(400)' \
          ' ,substrate TEXT, taxon_id INT, ' \
          ' description TEXT,organism TEXT, uniprot_gene TEXT, uniprot_annotation_score INT)'
    cursor.execute(sql,)
    conn.commit()

    sql2 = 'select uniprot_id from transporters.uniprot_table where uniprot_accession="%s"' % (uniprot_accession)

    try:
        cursor.execute(sql2,)
        return cursor.fetchall()[0][0]
    except:
        print '%s not in database' % uniprot_accession
        annot_dico = blastswiss2sqltable.get_swissprot_annotation([uniprot_accession])
        try:
            # try with uniprot
            taxon_id = annot_dico[uniprot_accession][0]
            uniprot_description = annot_dico[uniprot_accession][2]
            annot_score = annot_dico[uniprot_accession][1].split(' ')[0]
            uniprot_gene = annot_dico[uniprot_accession][3]
            uniprot_organism = annot_dico[uniprot_accession][4]
        except:
            try:
                # try with NCBI (removed uniprot entries)
                taxon_id = accession2taxon_id.accession2taxon_id([uniprot_accession])[uniprot_accession]
                data = accession2taxon_id.accession2description([uniprot_accession])
                uniprot_description = data[uniprot_accession]['description']
                uniprot_organism = data[uniprot_accession]['source']
                annot_score = 0
                uniprot_gene = '-'
            except:
                try:
                    # try with tcdb itself
                    data = tcdb_utils.accession2species_and_product(uniprot_accession, tc_id)
                    taxon_id = data[0].split('[')[1].split(']')[0]
                    uniprot_description = data[1]
                    uniprot_organism = data[0].split('[')[0]
                    annot_score = 0
                    uniprot_gene = '-'
                except:
                    # try to search with genus and species id
                    data = tcdb_utils.accession2species_and_product(uniprot_accession, tc_id)
                    uniprot_description = data[1]
                    uniprot_organism = data[0].split('[')[0]
                    genus_and_species = ' '.join(uniprot_organism.split(' ')[0:2])
                    taxon_id = accession2taxon_id.species_name2taxon_id(genus_and_species)
                    if not taxon_id:
                        # set cellular organism
                        taxon_id = 131567
                    annot_score = 0
                    uniprot_gene = '-'

        uniprot_description = re.sub('"', '', uniprot_description)
        #product = tcdb_utils.accession2species_and_product(uniprot_accession, tc_id)
        substrate = tcdb_utils.accession2substrate(uniprot_accession, tc_id)
        sql = 'insert into transporters.uniprot_table (uniprot_accession,taxon_id, substrate,description, organism, ' \
              ' uniprot_gene, uniprot_annotation_score) values ("%s", %s, "%s", "%s","%s","%s",%s)' % (uniprot_accession,
                                                                                                       taxon_id,
                                                                                                       substrate,
                                                                                                       uniprot_description,
                                                                                                       uniprot_organism,
                                                                                                       uniprot_gene,
                                                                                                       annot_score)
        cursor.execute(sql,)
        conn.commit()
        sql2 = 'select uniprot_id from transporters.uniprot_table where uniprot_accession="%s"' % (uniprot_accession)
        cursor.execute(sql2,)
        return cursor.fetchall()[0][0]


def import_annot(gblast_file, biodb, fasta_file_query, fasta_file_db, xml_dir):
    from chlamdb.biosqldb import manipulate_biosqldb
    from Bio import SeqIO
    from bs4 import BeautifulSoup
    import tcdb_utils
    import os
    from Bio.Blast import NCBIXML
    import re




    '''
    ajouter:
    - query length
    - hit length
    - query cov
    - hit cov

    Tables
    TC avec description, family id, family et substrats
    TCDB hit (accession)
    hits_chlamydia_04_16 avec details du hit (score,...)

    '''

    with open(fasta_file_db, 'r') as f:
        db_accession2record = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
    with open(fasta_file_query, 'r') as f:
        query_accession2record = SeqIO.to_dict(SeqIO.parse(f, 'fasta'))
    blast_data = {}
    xml = os.listdir(xml_dir)
    print 'n xml files', len(xml)
    for n,one_xml in enumerate(xml):
        if n % 100 == 0:
            print n, len(xml)
        #print os.path.join(xml_dir, one_xml)
        result_handle = open(os.path.join(xml_dir, one_xml), 'r')
        blast_records = NCBIXML.parse(result_handle)
        for record in blast_records:
            locus = record.query.split(' ')[0]
            blast_data[locus] = record

    server, db = manipulate_biosqldb.load_db(biodb)

    sql = 'create table IF NOT EXISTS transporters_transporters (seqfeature_id INT primary key,' \
          ' taxon_id INT,' \
          ' hit_uniprot_id INT,' \
          ' transporter_id INT,' \
          ' align_length INT,' \
          ' n_hsps INT,' \
          ' evalue FLOAT,' \
          ' bitscore_first_hsp FLOAT,' \
          ' identity FLOAT,' \
          ' query_TMS INT,' \
          ' hit_TMS INT,' \
          ' TM_overlap_score FLOAT,' \
          ' query_len INT,' \
          ' hit_len INT,' \
          ' query_cov FLOAT,' \
          ' hit_cov FLOAT,' \
          ' index transporter_id(transporter_id),' \
          ' index hit_uniprot_id(hit_uniprot_id));' % biodb

    server.adaptor.execute(sql,)
    server.commit()

    sql1 = 'select locus_tag, taxon_id from orthology_detail' % biodb
    sql2 = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id' % biodb

    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql1,))
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

    with open(gblast_file, 'r') as f:

        soup = BeautifulSoup(f, 'html.parser')
        table = soup.find('table')

        rows = table.find_all('tr')
        for n, one_row in enumerate(rows):
            cols = one_row.find_all('td')
            cols = [ele.text.strip() for ele in cols]



            if n == 0:
                pass
            else:

                locus_tag = cols[0]
                hit_uniprot_accession = cols[1]
                hit_tcid = cols[2]
                hit_description = re.sub('"', '', cols[3])
                hit_accession = hit_description.split(' ')[1]

                align_length = int(cols[4])
                evalue = cols[5]
                identity = cols[6]
                query_TMS = cols[7]
                hit_TMS = cols[8]
                TM_overlap_score = cols[9]
                if TM_overlap_score == "None":
                    print 'none!!!!!!!!!!'
                    TM_overlap_score = 0
                family_abrv_= cols[10]

                first_hsp = blast_data[locus_tag].alignments[0].hsps[0]

                #print blast_data[locus_tag].effective_query_length
                #print blast_data[locus_tag].effective_database_length
                query_align_length = first_hsp.query_end-first_hsp.query_start
                query_length = len(query_accession2record[locus_tag])
                hit_length = len(db_accession2record[hit_accession])

                hit_align_length = first_hsp.sbjct_end-first_hsp.sbjct_start
                query_cov = round(query_align_length/float(query_length), 2)
                hit_cov = round(hit_align_length/float(hit_length), 2)
                n_hsps = len(blast_data[locus_tag].alignments[0].hsps)
                bitscore_first_hsp = first_hsp.bits

                transporter_id = insert_new_tc_path(biodb, hit_tcid)
                hit_uniprot_id = insert_uniprot_in_db(biodb, hit_uniprot_accession, hit_tcid)


                sql = 'insert into transporters_transporters (seqfeature_id,' \
                      ' taxon_id,' \
                      ' hit_uniprot_id,' \
                      ' transporter_id,' \
                      ' align_length,' \
                      ' n_hsps,' \
                      ' evalue,' \
                      ' bitscore_first_hsp,' \
                      ' identity,' \
                      ' query_TMS,' \
                      ' hit_TMS,' \
                      ' TM_overlap_score,' \
                      ' query_len,' \
                      ' hit_len,' \
                      ' query_cov,' \
                      ' hit_cov) values (%s, %s, %s, %s, %s, %s, %s,%s,%s,%s,%s,%s,%s,%s,%s,%s)' % (biodb,
                                                                                                     locus_tag2seqfeature_id[locus_tag],
                                                                                                     locus_tag2taxon_id[locus_tag],
                                                                                                     hit_uniprot_id,
                                                                                                     transporter_id,
                                                                                                     align_length,
                                                                                                     n_hsps,
                                                                                                     evalue,
                                                                                                     bitscore_first_hsp,
                                                                                                     identity,
                                                                                                     query_TMS,
                                                                                                     hit_TMS,
                                                                                                     TM_overlap_score,
                                                                                                     query_length,
                                                                                                     hit_length,
                                                                                                     query_cov,
                                                                                                     hit_cov)

                print sql
                server.adaptor.execute(sql,)
        server.commit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--gblast_html', type=str, help="gblast html file")
    parser.add_argument("-b", '--fasta_database', type=str, help="fasta db (TCDB)")
    parser.add_argument("-f", '--fasta_query', type=str, help="fasta query")
    parser.add_argument("-x", '--xml_dir', type=str, help="xml directory (relative path)")
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)

    args = parser.parse_args()


    #insert_new_tc_path("chlamydia_04_16", "2.A.1.4.6")
    #print insert_uniprot_in_db("chlamydia_04_16", "F8L7K6", "9.B.163.1.5")
    import_annot(args.gblast_html, args.db_name,args.fasta_query, args.fasta_database, args.xml_dir)
