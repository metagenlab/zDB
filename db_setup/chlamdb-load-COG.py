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
def blast2COG(blast_file,
              hash2locus_tag_list,
              cdd_id2cog_id,
              cog_id2length,
              locus_tag2protein_length,
              coverage_cutoff=50,
              identity_cutoff=20):

    with open(blast_file, "r") as f:
        locus2data = {}
        #TODO: add bitscore
        '''
         0. 	 qseqid 	 query (e.g., gene) sequence id
         1. 	 sseqid 	 subject (e.g., reference genome) sequence id
         2. 	 pident 	 percentage of identical matches
         3. 	 length 	 alignment length
         4. 	 mismatch 	 number of mismatches
         5. 	 gapopen 	 number of gap openings
         6. 	 qstart 	 start of alignment in query
         7. 	 qend 	     end of alignment in query
         8. 	 sstart 	 start of alignment in subject
         9. 	 send 	     end of alignment in subject
         10. 	 evalue 	 expect value
         11. 	 bitscore 	 bit score
        '''
        for line in f:
            data = line.rstrip().split('\t')
            locus_tag_list = hash2locus_tag_list[data[0]]
            hit_cdd_id = data[1].split(":")[1]
            cog_id = cdd_id2cog_id[hit_cdd_id]
            identity = float(data[2])

            evalue = float(data[10])
            bitscore = float(data[11])

            query_start = int(data[6])
            query_end = int(data[7])
            hit_start = int(data[8])
            hit_end = int(data[9])

            query_align_length = (query_end-query_start) + 1
            hit_align_length = (hit_end-hit_start) + 1

            # locus aa sequence
            query_coverage = round((int(query_align_length)/int(locus_tag2protein_length[locus_tag_list[0]])) * 100, 2)

            # COG profile
            hit_coverage = round((int(hit_align_length)/int(cog_id2length[cog_id])) * 100, 2)

            if identity >= identity_cutoff and query_coverage >= coverage_cutoff and hit_coverage >= coverage_cutoff:
                # iter the list of locus corresponding to the same hash
                for locus_tag in locus_tag_list:
                    if locus_tag not in locus2data:
                        locus2data[locus_tag] = {}
                        locus2data[locus_tag]["hit_cdd_id"] = hit_cdd_id
                        locus2data[locus_tag]["cog_id"] = cog_id
                        locus2data[locus_tag]["cdd_id"] = hit_cdd_id
                        locus2data[locus_tag]["identity"] = identity
                        locus2data[locus_tag]["evalue"] = evalue
                        locus2data[locus_tag]["bitscore"] = bitscore
                        locus2data[locus_tag]["query_start"] = query_start
                        locus2data[locus_tag]["query_end"] = query_end
                        locus2data[locus_tag]["hit_start"] = hit_start
                        locus2data[locus_tag]["hit_start"] = hit_start
                        locus2data[locus_tag]["hit_end"] = hit_end
                        locus2data[locus_tag]["query_coverage"] = query_coverage
                        locus2data[locus_tag]["hit_coverage"] = hit_coverage

        return locus2data


def load_locus2cog_into_sqldb(input_blast_files,
                              biodb,
                              hash2locus_tag_list,
                              cdd_id2cog_id,
                              cog_id2length):
    import MySQLdb
    import os
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    '''
    locus2data[locus]["hit_cdd_id"] = hit_cdd_id
    locus2data[locus]["cog_id"] = cog_id
    locus2data[locus]["identity"] = identity
    locus2data[locus]["evalue"] = evalue
    locus2data[locus]["bitscore"] = bitscore
    locus2data[locus]["query_start"] = query_start
    locus2data[locus]["query_end"] = query_end
    locus2data[locus]["hit_start"] = hit_start
    locus2data[locus]["hit_end"] = hit_end
    locus2data[locus]["query_coverage"] = query_coverage
    locus2data[locus]["hit_coverage"] = hit_coverage
    '''
    # locus_tag2gi_hit_
    sql = 'create table COG_seqfeature_id2best_COG_hit (bioentry_id INT, ' \
          ' seqfeature_id INT, ' \
          ' hit_cog_id INT,' \
          ' cdd_id varchar(200), ' \
          ' query_start int,' \
          ' query_end int,' \
          ' hit_start int,' \
          ' hit_end int,' \
          ' query_coverage FLOAT,' \
          ' hit_coverage FLOAT,' \
          ' identity FLOAT,' \
          ' evalue FLOAT,' \
          ' bitscore FLOAT)'
    
    cursor.execute(sql)
    conn.commit()

    sql = 'select locus_tag,bioentry_id from orthology_detail t1 ' \
          ' inner join bioentry as t2 on t1.accession=t2.accession ' \
          ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
          ' where t3.name="%s"' % (biodb)
    sql2 = 'select protein_id, locus_tag from orthology_detail'
    sql3 = 'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id'
    sql4 = 'select COG_name,COG_id from COG_cog_names_2014'
    sql5 = 'select locus_tag,length(translation) from orthology_detail;'

    server, db = manipulate_biosqldb.load_db(biodb)
    locus2bioentry_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql))
    protein_id2locus_tag = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2))
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3))
    COG_name2COG_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql4))
    locus_tag2protein_length = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql5))

    for input_blast in input_blast_files:
        print ('file', input_blast)
        locus2data = blast2COG(input_blast,
                               hash2locus_tag_list,
                               cdd_id2cog_id,
                               cog_id2length,
                               locus_tag2protein_length)

        for locus in locus2data:

            if locus2data[locus]["cog_id"] not in COG_name2COG_id:
                print("COG %s not in COG_name2COG_id, probably removed..." % locus2data[locus]["cog_id"])
            else:
                sql = 'INSERT into COG_seqfeature_id2best_COG_hit (bioentry_id, ' \
                      ' seqfeature_id, ' \
                      ' hit_cog_id,' \
                      ' cdd_id, ' \
                      ' query_start,' \
                      ' query_end,' \
                      ' hit_start,' \
                      ' hit_end,' \
                      ' query_coverage,' \
                      ' hit_coverage,' \
                      ' identity,' \
                      ' evalue,' \
                      ' bitscore) VALUES (%s, %s, %s, "%s", %s, %s, %s, %s, %s, %s, %s, %s, %s)' % (locus2bioentry_id[locus],
                                                                                                    locus_tag2seqfeature_id[locus],
                                                                                                    COG_name2COG_id[locus2data[locus]["cog_id"]],
                                                                                                    locus2data[locus]["cdd_id"],
                                                                                                    locus2data[locus]["query_start"],
                                                                                                    locus2data[locus]["query_end"],
                                                                                                    locus2data[locus]["hit_start"],
                                                                                                    locus2data[locus]["hit_end"],
                                                                                                    locus2data[locus]["query_coverage"],
                                                                                                    locus2data[locus]["hit_coverage"],
                                                                                                    locus2data[locus]["identity"],
                                                                                                    locus2data[locus]["evalue"],
                                                                                                    locus2data[locus]["bitscore"])

            cursor.execute(sql)
    conn.commit()
    sql_index1 = 'create index csidbchs on COG_seqfeature_id2best_COG_hit(seqfeature_id)'
    sql_index2 = 'create index csidbchb on COG_seqfeature_id2best_COG_hit(bioentry_id)'
    sql_index3 = 'create index csidbchc on COG_seqfeature_id2best_COG_hit(cdd_id)'
    sql_index4 = 'create index csidbchh on COG_seqfeature_id2best_COG_hit(hit_cog_id)'

    cursor.execute(sql_index1)
    cursor.execute(sql_index2)
    cursor.execute(sql_index3)
    cursor.execute(sql_index4)
    conn.commit()


def load_locus2cog_into_sqldb_legacy(input_blast_files,
                                     biodb,
                                     cdd_id2cog_id,
                                     hash2locus_tag_list):
    import MySQLdb
    import os
    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(biodb)
    conn = server.adaptor.conn
    cursor = server.adaptor.cursor

    sql = 'create table COG_locus_tag2gi_hit (accession varchar(100), locus_tag varchar(100), gi INT, COG_id varchar(100))'
    
    cursor.execute(sql)
    conn.commit()
    sql = 'select locus_tag, accession from orthology_detail'
    sql2 = 'select protein_id, locus_tag from orthology_detail'
    sql5 = 'select locus_tag,length(translation) from orthology_detail;'

    server, db = manipulate_biosqldb.load_db(biodb)
    locus2genome_accession = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql))
    protein_id2locus_tag = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2))
    locus_tag2protein_length = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql5))

    for input_blast in input_blast_files:

        locus2data = blast2COG(input_blast,
                               hash2locus_tag_list,
                               cdd_id2cog_id,
                               cog_id2length,
                               locus_tag2protein_length)

        for locus in locus2data:
            sql = 'INSERT into COG_locus_tag2gi_hit (accession, locus_tag, gi, COG_id) VALUES ("%s", "%s", %s, "%s")' % (locus2genome_accession[locus],
                                                                                                                         locus,
                                                                                                                         0, # no gi anymore
                                                                                                                         locus2data[locus]["cog_id"])
            cursor.execute(sql)
    conn.commit()
    sql_index1 = 'create index cltghlt on COG_locus_tag2gi_hit(locus_tag)'
    sql_index2 = 'create index cltgha on COG_locus_tag2gi_hit(accession)'
    cursor.execute(sql_index1)
    cursor.execute(sql_index2)
    conn.commit()
    

if __name__ == '__main__':
    import argparse
    import chlamdb_setup_utils
    from Bio import SeqIO
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", '--input_blast', type=str, help="blast tab file", nargs='+')
    parser.add_argument("-d", '--database_name', type=str, help="database name", default=False)
    parser.add_argument("-u", '--hash2locus_tag', type=str, help="Tab separated file with correspondance between sequence hashes and locus tags")
    parser.add_argument("-cc", '--cog2cdd', type=str, help="cog2cdd")
    parser.add_argument("-cl", '--cdd2length', type=str, help="cdd2length")

    args = parser.parse_args()

    cdd_id2cog_id = {}
    with open(args.cog2cdd, 'r') as f:
        for row in f:
            data = row.rstrip().split("\t")
            cdd_id2cog_id[data[1]] = data[0]

    cog_id2length = {}
    with open(args.cdd2length, 'r') as f:
        for row in f:
            data = row.rstrip().split("\t")
            cog_id2length[data[0]] = data[1]

    hash2locus_list = chlamdb_setup_utils.get_hash2locus_list(args.hash2locus_tag)

    load_locus2cog_into_sqldb(args.input_blast,
                              args.database_name,
                              hash2locus_list,
                              cdd_id2cog_id,
                              cog_id2length)

    load_locus2cog_into_sqldb_legacy(args.input_blast,
                                        args.database_name,
                                        cdd_id2cog_id,
                                        hash2locus_list)

    manipulate_biosqldb.update_config_table(args.database_name, "COG_data")