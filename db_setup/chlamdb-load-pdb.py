#!/usr/bin/env python


def retrieve_pdb_data():

    import urllib.request

    '''
    0 IDCODE, 
    1 HEADER, 
    2 ACCESSION DATE, 
    3 COMPOUND, 
    4 SOURCE, 
    5 AUTHOR LIST, 
    6 RESOLUTION, 
    7 EXPERIMENT TYPE (IF NOT X-RAY)
    '''
    
    url = 'ftp://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx'
    
    handle = urllib.request.urlopen(url)
    pdb_table = handle.read().decode("UTF-8")
    pdb2data = {}
    for row in pdb_table.split('\n'):
        if len(row) == 0:
            continue
        if row[0] == '-':
            continue
        data = row.split("\t")
        if len(data) != 8:
            continue
        pdb2data[data[0]] = {}
        pdb2data[data[0]]["HEADER"] = data[1]
        pdb2data[data[0]]["accession_date"] = data[2]
        pdb2data[data[0]]["COMPOUND"] = data[3]
        pdb2data[data[0]]["SOURCE"] = data[4]
        pdb2data[data[0]]["RESOLUTION"] = data[6]
        pdb2data[data[0]]["METHOD"] = data[7]
    return pdb2data


def load_pdb_data(blast_output_files,
                  db_name,
                  hash2locus_list):
    
    from chlamdb.biosqldb import manipulate_biosqldb

    print("db conn:", db_name)
    server, db = manipulate_biosqldb.load_db(db_name)
    
    sql = f'select locus_tag, seqfeature_id from custom_tables_locus2seqfeature_id'
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    
    sql = f'create table if not exists custom_tables_seqfeature_id2pdb_BBH (seqfeature_id INT, pdb_accession varchar(200), header TEXT, compound TEXT, source TEXT, resolution FLOAT, method TEXT, identity FLOAT, score FLOAT, evalue FLOAT)'
    server.adaptor.execute(sql,)
    
    sql_template = f'insert into custom_tables_seqfeature_id2pdb_BBH values (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)'

    pdb2data = retrieve_pdb_data()

    problem_list = []
    
    for blast_file in blast_output_files:
        with open(blast_file, "r") as f:
            rows = [i.rstrip().split("\t") for i in f]
            count = 1
            for n, row in enumerate(rows):
                #print(row)
                if n % 1000000 == 0:
                    print(n)
                query_hash = row[0]
                previous_query_hash = rows[n-1][0]
                if query_hash == previous_query_hash:
                    if not match_problem:
                        continue
                    else:
                        match_problem = False
            
                hit_accession = row[1].split("_")[0].upper()
                e_value = row[10]
                score = row[11]
                identity = row[2]
                try:
                    header = pdb2data[hit_accession]["HEADER"]
                    compound = pdb2data[hit_accession]["COMPOUND"]
                    source = pdb2data[hit_accession]["SOURCE"]
                    resolution = pdb2data[hit_accession]["RESOLUTION"]
                    method = pdb2data[hit_accession]["METHOD"]
                    match_problem = False
                # in case of mapping problem, take the next BBH
                except KeyError:
                    print("problem with", hit_accession)
                    if hit_accession not in problem_list:
                        problem_list.append(hit_accession)
                    match_problem = True
                    continue
                try:
                    float(resolution)
                except:
                    resolution = None
                
                for locus in hash2locus_list[query_hash]:
                    server.adaptor.execute(sql_template, (locus_tag2seqfeature_id[locus],
                                                          hit_accession,
                                                          header,
                                                          compound,
                                                          source,
                                                          resolution,
                                                          method,
                                                          identity,
                                                          score,
                                                          e_value))
                    
    sql = f'create index ctsipsid on custom_tables_seqfeature_id2pdb_BBH (seqfeature_id)'                 
    server.adaptor.execute(sql,)
    server.commit()
    print(len(problem_list), problem_list)


if __name__ == '__main__':
    import argparse
    import chlamdb_setup_utils
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", '--blast_output', type=str, help="BLAST output file(s)", required=True, nargs="+")
    parser.add_argument("-d", '--db_name', type=str, help="db name", required=True)
    parser.add_argument("-c", '--corresp_table',  type=str, help="hash to locus correspondance table")



    args = parser.parse_args()
    
    hash2locus_list = chlamdb_setup_utils.get_hash2locus_list(args.corresp_table)


    load_pdb_data(args.blast_output, 
                    args.db_name,
                    hash2locus_list)
    
    manipulate_biosqldb.update_config_table(args.db_name, "PDB_data")
