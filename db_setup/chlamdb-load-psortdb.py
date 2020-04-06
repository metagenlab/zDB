#!/usr/bin/env python


def parse_psort_results(out_file):
    import re 
    
    '''
    SeqID: CRC-5EEC22EB4852E47C ELAC_p0001 Estrella lausannensis CRIB-30 plasmid 1
  Analysis Report:
    CMSVM-            Unknown                       [No details]
    CytoSVM-          Unknown                       [No details]
    ECSVM-            Unknown                       [No details]
    ModHMM-           Unknown                       [No internal helices found]
    Motif-            Unknown                       [No motifs found]
    OMPMotif-         Unknown                       [No motifs found]
    OMSVM-            Unknown                       [No details]
    PPSVM-            Unknown                       [No details]
    Profile-          Unknown                       [No matches to profiles found]
    SCL-BLAST-        Unknown                       [No matches against database]
    SCL-BLASTe-       Unknown                       [No matches against database]
    Signal-           Unknown                       [No signal peptide detected]
  Localization Scores:
    Cytoplasmic            2.00
    OuterMembrane          2.00
    Extracellular          2.00
    CytoplasmicMembrane    2.00
    Periplasmic            2.00
  Final Prediction:
    Unknown
    '''
   
    pattern_seqid = re.compile("SeqID: ([A-Za-z0-9\-]+) .*") 
    pattern_score = re.compile("Localization Scores:") 
    pattern_final_prediction = re.compile("Final Prediction:")
    seq_hash2data = {}
    with open(out_file, 'r') as f:
        rows = [i.rstrip() for i in f]
        for row in rows:
            data = row.split()
            if len(data) == 0:
                continue
            if pattern_seqid.search(row):
                seq_hash = pattern_seqid.search(row).group(1)
                final_prediction = False
                localization = False
            if pattern_score.search(row):
                localization = True 
                seq_hash2data[seq_hash] = {}
                continue
            if pattern_final_prediction.search(row):
                localization = False 
                final_prediction = True
                continue
            if localization:
                seq_hash2data[seq_hash][data[0]] = data[1]
            if final_prediction:
                seq_hash2data[seq_hash]["final_prediction"] = data[0]
                final_prediction = False
    return seq_hash2data


def load_sport_data(psortdb_output_files,
                    db_name,
                    hash2locus_list):

    from chlamdb.biosqldb import manipulate_biosqldb
    server, db = manipulate_biosqldb.load_db(db_name)

    hash2psort_results = {}
    for psortdb_output in psortdb_output_files:
        hash2psort_results.update(parse_psort_results(psortdb_output))

    sql = f'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id'
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    
    sql = f'create table if not exists custom_tables_seqfeature_id2psortb (seqfeature_id INT, final_prediction varchar(200), score FLOAT)'
    server.adaptor.execute(sql,)
    
    sql_template = f'insert into custom_tables_seqfeature_id2psortb values (%s, %s, %s)'
    
    for hash in hash2psort_results:
        for locus in hash2locus_list[hash]:
            final_pred = hash2psort_results[hash]['final_prediction']
            if final_pred == 'Unknown':
                score = None
            else:
                score = hash2psort_results[hash][final_pred]
            server.adaptor.execute(sql_template, (locus_tag2seqfeature_id[locus],
                                                  final_pred,
                                                  score))
    
    sql = f'create index ctsipsid on custom_tables_seqfeature_id2psortb (seqfeature_id)'                 
    server.adaptor.execute(sql,)
    server.commit()


if __name__ == '__main__':
    import argparse
    import chlamdb_setup_utils
    from chlamdb.biosqldb import manipulate_biosqldb
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--psortdb_output', type=str, help="BLAST output file(s)", required=True, nargs="+")
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)
    parser.add_argument("-c", '--corresp_table',  type=str, help="hash to locus correspondance table")



    args = parser.parse_args()
    
    hash2locus_list = chlamdb_setup_utils.get_hash2locus_list(args.corresp_table)

    load_sport_data(args.psortdb_output, 
                    args.db_name,
                    hash2locus_list)
    
    manipulate_biosqldb.update_config_table(args.db_name, "psortb_data")