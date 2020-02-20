#! /usr/bin/env python


def gbk2taxid(gbk_files, db_name):
    from Bio import SeqIO
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(db_name)
    file_name2taxid = {}

    for gbk in gbk_files:
        records = [i for i in SeqIO.parse(gbk, 'genbank')]
        sql = 'select taxon_id from bioentry t1 ' \
              ' inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
              ' where t1.accession="%s" and t2.name="%s";' % (records[0].name, db_name)
        print(sql)
        file_name2taxid[gbk.split("/")[-1].split(".")[0]] = int(server.adaptor.execute_and_fetchall(sql,)[0][0])
    return file_name2taxid


def import_checkm(checkm_file,
                  biodb, 
                  gbk_files):
    
    from chlamdb.biosqldb import manipulate_biosqldb
    import re

    server, db = manipulate_biosqldb.load_db(biodb)
    
    file_name2taxid = gbk2taxid(gbk_files, biodb)

    #   Bin Id             Marker lineage   # genomes   # markers   # marker sets   0     1    2    3   4   5+   Completeness   Contamination   Strain heterogeneity
    sql = 'create table IF NOT EXISTS custom_tables_checkm (taxon_id INT, ' \
          ' n_missing INT,' \
          ' n_1 INT,' \
          ' n_2 INT,' \
          ' n_3 INT,' \
          ' n_4 INT,' \
          ' n_5_plus INT,' \
          ' n_total INT,' \
          ' completeness FLOAT,' \
          ' contamination FLOAT,' \
          ' heterogeneity FLOAT);'
          
    server.adaptor.execute(sql,)
    with open(checkm_file, 'r') as f:
        for n, one_row in enumerate(f):
            data = re.split('\s+', one_row.rstrip())
            '''
            0    Bin Id	
            1    Marker lineage	
            2    # genomes	
            3    # markers	
            4    # marker sets	
            5    Completeness	
            6    Contamination	
            7    Strain heterogeneity	
            8    Genome size (bp)	
            9    # ambiguous bases	
            10    # scaffolds	
            11    # contigs	
            12    N50 (scaffolds)	
            13    N50 (contigs)	
            14    Mean scaffold length (bp)
            15    Mean contig length (bp)	
            16    Longest scaffold (bp)	
            17    Longest contig (bp)
            18    GC	GC std (scaffolds > 1kbp)	
            19    Coding density	
            20    Translation table	
            21    # predicted genes	
            22    0	
            23    1	
            24    2	
            25    3	
            26    4	
            27    5+
            
            '''
            if len(data) != 1 and data[0] != 'Bin':
                accession = data[0]
                n_missing = data[23]
                n_1 = data[24]
                n_2 = data[25]
                n_3 = data[26]
                n_4 = data[27]
                n_5_plus = data[28]
                completeness = data[5]
                contamination = data[6]
                heterogeneity = data[7]
                n_sum = int(n_2) + int(n_3)+ int(n_4) + int(n_5_plus)

                sql = 'insert into  custom_tables_checkm values (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s);' % (file_name2taxid[accession],
                                                                                                                  n_missing,
                                                                                                                  n_1,
                                                                                                                  n_2,
                                                                                                                  n_3,
                                                                                                                  n_4,
                                                                                                                  n_5_plus,
                                                                                                                  n_sum,
                                                                                                                  completeness,
                                                                                                                  contamination,
                                                                                                                  heterogeneity)

                server.adaptor.execute(sql,)
        server.commit()

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--checkm_table', type=str, help="check table")
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)
    parser.add_argument("-g", '--gbk_files', type=str, help="GBK files", nargs='+')

    args = parser.parse_args()

    import_checkm(args.checkm_table, 
                  args.db_name, 
                  args.gbk_files)
