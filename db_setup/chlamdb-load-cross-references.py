#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

# Author: Trestan Pillonel (trestan.pillonel[]gmail.com)
# Date: 2019
# ---------------------------------------------------------------------------


def parse_idmapping_crossrefs(table):
    # CRC-84A6E563630AFDE9    A0A317JCC6      UniProtKB-ID    A0A317JCC6_9BACT
    uniprot_accession2crosserfs = {}
    with open(table, 'r') as f:
        for row in f:
            data = row.rstrip().split("\t")
            uniprot_accession = data[1]
            database_name = data[2]
            accession = data[3]
            if uniprot_accession not in uniprot_accession2crosserfs:
                uniprot_accession2crosserfs[uniprot_accession] = [[database_name, accession]]
            else:
                uniprot_accession2crosserfs[uniprot_accession].append([database_name, accession])
    return uniprot_accession2crosserfs


def create_db(biodatabase):

    server, db = manipulate_biosqldb.load_db(biodatabase)

    # create crossrefs table
    sql = f'CREATE TABLE IF NOT EXISTS biosqldb.cross_references_{biodatabase} (seqfeature_id INT,'\
          f' db_name VARCHAR(200),' \
          f' accession varchar(200),' \
          f' index seqfeature_id (seqfeature_id), index accession (accession));'

    server.adaptor.execute_and_fetchall(sql,)

    # disallow ducpliucate seqfeature id - accession pairs
    sql = f'ALTER TABLE biosqldb.cross_references_{biodatabase} ADD UNIQUE (seqfeature_id, accession);'
    server.adaptor.execute_and_fetchall(sql,) 
     

def unirpot_crossrefs(biodatabase, 
                      hash2locus_list, 
                      uniprot_crossrefs_table):

    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodatabase)
          
    create_db(biodatabase)
    
    # get uniprot accession list and locus2seqfeature_id dictionnary
    sql2 = f'select seqfeature_id,uniprot_accession from custom_tables.uniprot_id2seqfeature_id_{biodatabase};'
    sql3 = f'select locus_tag,seqfeature_id from annotation.seqfeature_id2locus_{biodatabase};'
    sql4 = f'select seqfeature_id,protein_id from annotation.seqfeature_id2CDS_annotation_{biodatabase}'

    seqfeature_id2uniprot_accession = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2))
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3))
    seqfeature_id2protein_accession = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql4))

    # parse crossrefs table
    uniprot_accession2crosserfs = parse_idmapping_crossrefs(uniprot_crossrefs_table)
    
    for hash in hash2locus_list:
        for locus in list(set(hash2locus_list[hash])):
            seqfeature_id = locus_tag2seqfeature_id[locus]
            # add seqfeature_id itself
            sql = f'insert into biosqldb.cross_references_{biodatabase} (seqfeature_id, db_name, accession) values ({seqfeature_id}, "seqfeature_id", "{seqfeature_id}")'
            server.adaptor.execute(sql,)          
            # add locus
            sql = f'insert into biosqldb.cross_references_{biodatabase} (seqfeature_id, db_name, accession) values ({seqfeature_id}, "Locus tag", "{locus}")'
            server.adaptor.execute(sql,)
            # add protein accession
            prot = seqfeature_id2protein_accession[str(seqfeature_id)]
            
            # can happen with prokka annotation
            if prot != locus:
                sql = f'insert into biosqldb.cross_references_{biodatabase} (seqfeature_id, db_name, accession) values ({seqfeature_id}, "Protein Accession", "{prot}")'
                server.adaptor.execute(sql,)
    
            # get uniprot crossrefs and insert
            try:
                uniprot_accession = seqfeature_id2uniprot_accession[str(seqfeature_id)]
            except KeyError:
                continue
            uniprot_crossrefs = uniprot_accession2crosserfs[uniprot_accession]
            
            # add alll crossrefs
            for uniprot_crossref in uniprot_crossrefs:
                db_name = uniprot_crossref[0]
                db_acc = uniprot_crossref[1]

                # if match to any of the locus tag of the db (primary key), skip
                if db_acc in locus_tag2seqfeature_id:
                    continue
                
                # skip if the accession is the same as Genbank
                if db_acc == prot:
                    continue

                sql = f'insert into biosqldb.cross_references_{biodatabase} (seqfeature_id, db_name, accession) values ({seqfeature_id}, "{db_name}", "{db_acc}")'
                try:
                    server.adaptor.execute(sql,)
                except:
                    print("Duplicate key:", seqfeature_id, db_acc)
                # if format accession.version add accession as well
                accession_and_version = db_acc.split(".")
                if len(accession_and_version) == 2 and db_name != "STRING":
                    accession_no_version = accession_and_version[0]
                    sql = f'insert into biosqldb.cross_references_{biodatabase} (seqfeature_id, db_name, accession) values ({seqfeature_id}, "{db_name}", "{accession_no_version}")'
                    try:
                        server.adaptor.execute(sql,)
                    except:
                        print("problem with", sql)
                        continue
            # add uniprot accession itself
            try:
                sql = f'insert into biosqldb.cross_references_{biodatabase} (seqfeature_id, db_name, accession) values ({seqfeature_id}, "UniProtKB-Accession", "{uniprot_accession}")'
                server.adaptor.execute(sql,)
            except:
                pass
            # add locus_tag             

        server.commit()


def refseq_crossrefs(biodatabase,
                     refseq_table):
    # "old genbank locus tag"      "refseq locus_tag"    "refseq protein id"
    #  cpL17_0003	                cpL17_RS00015	      WP_081194193.1
    from chlamdb.biosqldb import manipulate_biosqldb

    server, db = manipulate_biosqldb.load_db(biodatabase)
    
    create_db(biodatabase)

    # get locus2seqfeature_id dictionnary
    sql = f'select locus_tag,seqfeature_id from annotation.seqfeature_id2locus_{biodatabase};'
    locus_tag2seqfeature_id = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql))

    with open(refseq_table, "r") as f:
        for row in f:
            data = row.rstrip().split("\t")
            biodb_locus = data[0]
            try:
                seqfeature_id = locus_tag2seqfeature_id[data[0]]
            except:
                print("Pseudogene?", seqfeature_id, biodb_locus)
            refseq_locus = data[1]
            refseq_protein_accession = data[2]
            refseq_protein_accession_no_version = data[2].split(".")[0]
            
            sql1 =  f'insert into biosqldb.cross_references_{biodatabase} (seqfeature_id, db_name, accession) values ({seqfeature_id}, "RefSeq locus_tag", "{refseq_locus}")'
            sql2 =  f'insert into biosqldb.cross_references_{biodatabase} (seqfeature_id, db_name, accession) values ({seqfeature_id}, "RefSeq", "{refseq_protein_accession}")'
            sql3 =  f'insert into biosqldb.cross_references_{biodatabase} (seqfeature_id, db_name, accession) values ({seqfeature_id}, "RefSeq", "{refseq_protein_accession_no_version}")'
            try:
                server.adaptor.execute(sql1,)
            except:
                print("Duplicate key:", seqfeature_id, refseq_locus)
            try:
                server.adaptor.execute(sql2,) 
            except:
                print("Duplicate key:", seqfeature_id, refseq_protein_accession)
            try:
                server.adaptor.execute(sql3,) 
            except:
                print("Duplicate key:", seqfeature_id, refseq_protein_accession_no_version)
    server.commit()


if __name__ == '__main__':
    import argparse
    from chlamdb.biosqldb import manipulate_biosqldb
    import chlamdb_setup_utils

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", '--idmapping', type=str, help="idmapping crossrefs")
    parser.add_argument("-r", '--refseq', type=str, help="RefSeq old_locus mapping")
    parser.add_argument("-d", '--database_name', type=str, help="database name")
    parser.add_argument("-c", '--corresp_table', type=str, help="hash to locus correspondance table")

    args = parser.parse_args()

    hash2locus_list = chlamdb_setup_utils.get_hash2locus_list(args.corresp_table)

    if args.idmapping:
        unirpot_crossrefs(args.database_name, 
                        hash2locus_list, 
                        args.idmapping)

    if args.refseq:
        refseq_crossrefs(args.database_name, 
                         args.refseq)       