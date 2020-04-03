#!/usr/bin/env python

def import_annot(gblast_file,
                 biodb,
                 fasta_file_query,
                 fasta_file_db,
                 xml_dir,
                 hash2locus_list):

    from chlamdb.biosqldb import manipulate_biosqldb
    from Bio import SeqIO
    from bs4 import BeautifulSoup
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
    print ('n xml files', len(xml))
    for n, one_xml in enumerate(xml):
        if n % 100 == 0:
            print (n, len(xml))
        result_handle = open(os.path.join(xml_dir, one_xml), 'r')
        blast_records = NCBIXML.parse(result_handle)
        for record in blast_records:
            locus = record.query.split(' ')[0]
            blast_data[locus] = record

    server, db = manipulate_biosqldb.load_db(biodb)
    
    sql = 'select protein_accession, protein_id from transporters_protein_entry'
    
    protein_accession2protein_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    
    # add accession without version to the dictionnary
    for protein in list(protein_accession2protein_id.keys()):
        if '.' in protein:
            accession_no_version = protein.split(".")[0]
            protein_accession2protein_id[accession_no_version] = protein_accession2protein_id[protein]

    sql = 'create table IF NOT EXISTS transporters_transporters_BBH (seqfeature_id INTEGER primary key,' \
          ' taxon_id INT,' \
          ' hit_protein_id INT,' \
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
          ' hit_cov FLOAT);'


    server.adaptor.execute(sql,)
    server.commit()

    sql1 = 'select locus_tag, taxon_id from orthology_detail'
    sql2 = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id'

    locus_tag2taxon_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql1,))
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

    with open(gblast_file, 'r') as f:

        soup = BeautifulSoup(f, 'html.parser')
        table = soup.find('table')

        rows = table.find_all('tr')
        for n, one_row in enumerate(rows):
            
            '''
            0 Query ID	
            1 Hit ID	
            2 Hit TCID	
            3 Hit Description	
            4 Match Len	
            5 e-Val	
            6 % Identity	
            7 Query Length	
            8 Hit Length	
            9 Query Coverage	
            10 Hit Coverage	
            11 Query TMS	
            12 Hit TMS	
            13 TM-Overlap Score	
            14 Family Abrv.	
            15 Predicted Substrate
            '''
            
            cols = one_row.find_all('td')
            cols = [ele.text.strip() for ele in cols]

            if n == 0:
                pass
            else:

                hash = cols[0]
                hit_uniprot_accession = cols[1]
                #hit_tcid = cols[2]
                hit_description = re.sub('"', '', cols[3])
                hit_accession = hit_description.split(' ')[1]

                align_length = int(cols[4])
                evalue = cols[5]
                identity = cols[6]
                query_TMS = cols[11]
                hit_TMS = cols[12]
                TM_overlap_score = cols[13]
                
                if TM_overlap_score == "None":
                    TM_overlap_score = 0
                family_abrv = cols[14]

                for locus_tag in hash2locus_list[hash]:

                    first_hsp = blast_data[hash].alignments[0].hsps[0]
                    query_align_length = first_hsp.query_end-first_hsp.query_start
                    query_length = len(query_accession2record[locus_tag])
                    hit_length = len(db_accession2record[hit_accession])

                    hit_align_length = first_hsp.sbjct_end-first_hsp.sbjct_start
                    query_cov = round(query_align_length/float(query_length), 2)
                    hit_cov = round(hit_align_length/float(hit_length), 2)
                    n_hsps = len(blast_data[hash].alignments[0].hsps)
                    bitscore_first_hsp = first_hsp.bits
                    
                    print(baba)
                    
                    
                    try:
                        hit_protein_id = protein_accession2protein_id[hit_uniprot_accession]
                    except KeyError:
                        print("WARNING accession %s was not found in the database, might have been removed from TCDB, skipping" % hit_uniprot_accession)
                        continue

                    sql = 'insert into transporters_transporters_BBH (seqfeature_id,' \
                            ' taxon_id,' \
                            ' hit_protein_id,' \
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
                            ' hit_cov) values (%s, %s, %s, %s, %s, %s, %s,%s,%s,%s,%s,%s,%s,%s,%s)' % (locus_tag2seqfeature_id[locus_tag],
                                                                                                        locus_tag2taxon_id[locus_tag],
                                                                                                        hit_protein_id,
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

                    server.adaptor.execute(sql,)
    server.commit()
    sql_index1 = 'create index tttid on transporters_transporters_BBH(hit_protein_id)'
    sql_index2 = 'create index ttthuid on transporters_transporters_BBH(hit_protein_id)'
    try:
        server.adaptor.execute(sql_index1,)
        server.adaptor.execute(sql_index2,)
    except:
        pass
    server.commit()



if __name__ == '__main__':
    import argparse
    import chlamdb_setup_utils

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--gblast_html', type=str, help="gblast html file")
    parser.add_argument("-b", '--fasta_database', type=str, help="fasta db (TCDB)")
    parser.add_argument("-f", '--fasta_query', type=str, help="fasta query")
    parser.add_argument("-x", '--xml_dir', type=str, help="xml directory (relative path)")
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)
    parser.add_argument("-c", '--corresp_table', type=str, help="hash to locus correspondance table")

    args = parser.parse_args()

    hash2locus_list = chlamdb_setup_utils.get_hash2locus_list(args.corresp_table)

    import_annot(args.gblast_html,
                 args.db_name,
                 args.fasta_query,
                 args.fasta_database,
                 args.xml_dir,
                 hash2locus_list)
