import os
import sys
import re

import chlamdb
import pandas as pd

from Bio import SeqIO
from Bio import AlignIO

# assumes orthofinder named: OG000N
# returns the N as int
def get_og_id(string):
    return int(string[2:])

def parse_orthofinder_output_file(output_file):
    protein_id2orthogroup_id = {}
    parsing = open(output_file, 'r')

    for line in parsing:
        tokens = line.strip().split(' ')

        # Skips the ":" at the end of the orthgroup id
        group = get_og_id(tokens[0][:-1])
        for locus in tokens[1:]:
            protein_id2orthogroup_id[locus] = group
    return protein_id2orthogroup_id

# TODO: import the two following functions into the chlamdb file to remove
# all database code from this file
def setup_biodb(kwargs):
    sqlpsw = os.environ['SQLPSW']
    db_type = kwargs["chlamdb.db_type"]
    db_name = kwargs["chlamdb.db_name"]
    schema_dir = kwargs["chlamdb.biosql_schema_dir"]
    err_code = 0

    if db_type=="sqlite":
        import sqlite3
        conn = sqlite3.connect(db_name)
        cursor = conn.cursor()
        url_biosql_scheme = 'biosqldb-sqlite.sql'
        err_code = os.system(f"sqlite3 {db_name} < {schema_dir}/{url_biosql_scheme}")
    else:
        import MySQLdb
        conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                    user="root", # your username
                                    passwd=sqlpsw) # name of the data base
        cursor = conn.cursor()
        sql_db = f'CREATE DATABASE IF NOT EXISTS {db_name};'
        cursor.execute(sql_db,)
        conn.commit()
        cursor.execute(f"use {db_name};",)
        url_biosql_scheme = 'biosqldb-mysql.sql'
        err_code = os.system(f"mysql -uroot -p{sqlpsw} {db_name} < {schema_dir}/{url_biosql_scheme}")

    # not really logical to me, but creating a database
    # from the biosql is necessary
    db = chlamdb.DB.load_db(kwargs)
    db.create_biosql_database(kwargs)
    db.commit()
    if err_code != 0:
        raise IOError("Problem loading sql schema:", err_code)

def create_data_table(kwargs):
    db_type = kwargs["chlamdb.db_type"]
    db_name = kwargs["chlamdb.db_name"]
    if db_type=="sqlite":
        import sqlite3
        conn = sqlite3.connect(db_name)
        cursor = conn.cursor()
    else:
        import os
        import MySQLdb
        sqlpsw = os.environ['SQLPSW']

        conn = MySQLdb.connect(host="localhost", # your host, usually localhost
                                    user="root", # your username
                                    passwd=sqlpsw, # your password
                                    db=db_name) # name of the data base
        cursor = conn.cursor()
    entry_list = [
        ("minimal", "mandatory", False),
        ("gbk_files", "mandatory", False),
        ("orthology_data", "mandatory", False),
        ("orthology_comparative", "mandatory", False),
        ("orthology_comparative_accession", "mandatory", False),
        ("orthology_consensus_annotation", "mandatory", False),
        ("orthogroup_alignments", "mandatory", False),
        ("old_locus_table", "mandatory", False),
        ("reference_phylogeny", "mandatory", False),
        ("taxonomy_table", "mandatory", False),
        ("genome_statistics", "mandatory", False),
        ("BLAST_database", "optional", False),
        ("gene_phylogenies", "optional", False),
        ("interpro_data", "optional", False),
        ("interpro_comparative", "optional", False),
        ("interpro_comparative_accession", "optional", False),
        ("priam_data", "optional", False),
        ("priam_comparative", "optional", False),
        ("priam_comparative_accession", "optional", False),
        ("COG_data", "optional", False),
        ("COG_comparative", "optional", False),
        ("COG_comparative_accession", "optional", False),
        ("KEGG_data", "optional", False),
        ("KEGG_comparative", "optional", False),
        ("KEGG_comparative_accession", "optional", False),
        ("pfam_comparative", "optional", False),
        ("pfam_comparative_accession", "optional", False),       
        ("TCDB_data", "optional", False),
        ("psortb_data", "optional", False),
        ("T3SS_data", "optional", False),
        ("PDB_data", "optional", False),
        ("BLAST_refseq", "optional", False),
        ("BLAST_swissprot", "optional", False),
        ("BBH_phylogenies", "optional", False),
        ("GC_statistics", "optional", False),
        ("gene_clusters", "optional", False),
        ("phylogenetic_profile", "optional", False),
        ("synonymous_table", "optional", False),
        ("interpro_taxonomy", "optional", False), # interpro taxnonomy statistics
        ("pfam_taxonomy", "optional", False), #  taxnonomy statistics
        ("COG_taxonomy", "optional", False) # COG taxnonomy statistics
    ]
    
    sql = 'create table biodb_config (name varchar(200), type varchar(200), status BOOLEAN)'
    
    cursor.execute(sql)
    conn.commit()
    
    sql = 'insert into biodb_config values ("%s", "%s", %s)'
    for row in entry_list:
        cursor.execute(sql % (row[0], row[1], row[2]),)
    conn.commit()
    
def insert_gbk(db, one_gbk):
    input_file = open(one_gbk, "r")
    records = SeqIO.parse(input_file, 'genbank')
    for record in records:
        params = {}
        params["accession"] = record.id.split('.')[0]
        params["taxon_id"] = db.get_taxid_from_accession(params["accession"])
        for feature in record.features:
            if feature.type == 'CDS' and not 'pseudo' in feature.qualifiers:
                params["start"] = re.sub('>|<','', str(feature.location.start))
                params["end"] = re.sub('>|<','', str(feature.location.end))
                params["strand"] = feature.strand
                params["gene"] = feature.qualifiers.get("gene", ['-'])[0]
                params["product"] = feature.qualifiers.get("product", ['-'])[0]
                params["locus_tag"] = feature.qualifiers["locus_tag"][0]
                params["old_locus_tag"] = feature.qualifiers.get("old_locus_tag", [False])[0]
                params["protein_id"] = feature.qualifiers.get("protein_id", ['.'])[0]
                params["translation"] = feature.qualifiers.get("translation", ['-'])[0]
                db.insert_cds(params)
            elif feature.type == 'rRNA':
                params["start"] = re.sub('>|<','', str(feature.location.start))
                params["end"] = re.sub('>|<','', str(feature.location.end))
                params["strand"] = feature.strand
                params["product"] = feature.qualifiers.get("product", ["-"])[0]
                db.insert_rRNA_in_feature_table(params)

def setup_chlamdb(**kwargs):
    setup_biodb(kwargs)
    create_data_table(kwargs)

def load_gbk(gbks, args):
    db = chlamdb.DB.load_db(args)

    db.create_cds_tables()
    for gbk in gbks:
        records = [i for i in SeqIO.parse(gbk, 'genbank')]
        db.load_gbk_wrapper(records)
        assert(len(records) == 1)
        insert_gbk(db, gbk)
    db.set_status_in_config_table("gbk_files", 1)
    db.create_indices_on_cds()
    db.commit()


def load_orthofinder_results(orthofinder_output, args):
    db = chlamdb.DB.load_db(args)
    orthogroup_feature_id = db.setup_orthology_table()
    hsh_prot_to_group = parse_orthofinder_output_file(orthofinder_output)
    hsh_locus_to_feature_id = db.get_hsh_locus_to_seqfeature_id()
    db.add_orthogroups_to_seq(hsh_prot_to_group, hsh_locus_to_feature_id, orthogroup_feature_id)

    arr_cnt_tables = db.get_orthogroup_count_table()

    db.create_orthology_table(arr_cnt_tables)
    db.load_orthology_table(arr_cnt_tables)
    # db.create_locus_to_feature_table()
    db.set_status_in_config_table("orthology_data", 1)
    db.commit()

# Note: as this is an alignment, the lengths are the same
def get_identity(seq1, seq2):
    identity = 0
    aligned = 0
    identical = 0
    gaps_1 = 0
    gaps_2 = 0

    assert(len(seq1) == len(seq2))
    for i in range(len(seq1)):
        if seq1[i]=="-":
            gaps_1 += 1
        if seq2[i]=="-":
            gaps_2 += 1
        if seq1[i]=="-" or seq2[i]=="-":
            continue
        if seq1[i]==seq2[i]:
            identical += 1
        aligned += 1

    if aligned/(len(seq1)-gaps_1) < 0.3 or aligned/(len(seq2)-gaps_2) < 0.3:
        return 0
    return 100*(identical/float(aligned))

def load_refseq_matches(args, diamond_tsvs):
    db = chlamdb.DB.load_db(args)
    db.create_refseq_hits_table()

    # map accessions to id
    sseqid_hsh = {}
    sseqid_id = 0

    # Note : add multithreading here?
    query_hash_to_seqfeature_id = db.hash_to_seqfeature()
    for tsv in diamond_tsvs:
        hit_table = pd.read_csv(tsv, sep="\t")
        hit_count = 0
        seqfeature_id = None
        query_hash = None
        data = []
        for index, row in hit_table.iterrows():
            # remove version number
            match_accession = row[1].split(".")[0]
            if match_accession not in sseqid_hsh:
                sseqid_hsh[match_accession] = sseqid_id
                sseqid_id += 1

            if query_hash==row[0]:
                hit_count += 1
            else:
                hit_count = 0
                query_hash = row[0]
                seqfeature_id = query_hash_to_seqfeature_id[query_hash]
            lst_args = row.tolist()
            data.append([hit_count, seqfeature_id, sseqid_hsh[match_accession]] + lst_args[2:])
        db.load_refseq_hits(data)
    # see which indices are better fitted
    # db.create_refseq_hits_indices()
    db.commit()
    return sseqid_hsh

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]

def load_refseq_matches_infos(args, hsh_sseqids):
    db = chlamdb.DB.load_db(args)
    db.create_diamond_refseq_match_id()

    # ugly
    for chunk in chunks(list(hsh_sseqids.keys()), 2000):
        hsh_accession_to_taxid = db.get_accession_to_taxid(chunk, args)
        hsh_accession_to_prot = db.get_accession_to_prot(chunk, args)

        data = []
        # Note: need to add some error handling in case the accession is not found
        # in the databases
        for key in chunk:
            taxid = hsh_accession_to_taxid[key]
            description, length = hsh_accession_to_prot[key]
            data.append([hsh_sseqids[key], key, taxid, description, length])
        db.load_diamond_refseq_match_id(data)
    db.commit()

def load_seq_hashes(args, nr_mapping, nr_fasta):
    db = chlamdb.DB.load_db(args)
    fasta_dict = SeqIO.to_dict(SeqIO.parse(nr_fasta, "fasta"))
    hsh_locus_to_id = db.get_hsh_locus_to_seqfeature_id()

    to_load = []
    for line in open(nr_mapping, "r"):
        record_id, hsh, genome = line.split("\t")
        sequence = str(fasta_dict[hsh].seq)
        seqfeature_id = hsh_locus_to_id[record_id]
        to_load.append( (seqfeature_id, hsh, sequence) )
    db.create_hash_table(to_load)
    db.commit()

def load_alignments_results(args, alignment_files):
    db = chlamdb.DB.load_db(args)
    db.create_new_og_matrix()
    locus_to_feature_id = db.get_hsh_locus_to_seqfeature_id()

    # assumes filename of the format OG00N_mafft.faa, with the orthogroup
    # being the integer following the OG string
    matrix = []
    averages = []
    for alignment in alignment_files:
        align = AlignIO.read(alignment, "fasta")
        orthogroup = get_og_id(alignment.split("_")[0])
        ttl = 0.0
        n = 0
        for i in range(len(align)):
            for j in range(i+1, len(align)):
                alignment_1 = align[i]
                alignment_2 = align[j]
                id_1 = locus_to_feature_id[alignment_1.name]
                id_2 = locus_to_feature_id[alignment_2.name]
                identity = get_identity(alignment_1, alignment_2)
                ttl += identity
                n += 1
                matrix.append( (orthogroup, id_1, id_2, identity) )
        averages.append((orthogroup, float(ttl)/ n))
    db.load_og_matrix(matrix)
    db.load_og_averages(averages)
    db.create_og_matrix_indices()
    db.set_status_in_config_table("orthogroup_alignments", 1)
    db.commit()
