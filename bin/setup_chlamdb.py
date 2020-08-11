import os
import sys
import re

from chlamdb_utils import db_utils

import pandas as pd

from Bio import SeqIO
from Bio import AlignIO
from Bio import SeqUtils


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
        conn.execute("pragma journal_mode=wal")
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
    db = db_utils.DB.load_db(kwargs)
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
        # Done
        ("gbk_files", "mandatory", False),

        # Done
        ("orthology_data", "mandatory", False),

        # Done: may be removed as the data is already in the database and the result
        # can be generated on the fly.
        # Note: not all data was generated as it can be generated on the fly.
        ("orthology_comparative", "mandatory", False),

        # this is one load tables that can generated on the fly
        # will need to write the queries
        ("orthology_comparative_accession", "mandatory", False),

        # this is one load tables that can generated on the fly
        # will need to write the queries
        ("orthology_consensus_annotation", "mandatory", False),

        # Done
        ("orthogroup_alignments", "mandatory", False),

        # wut?
        ("old_locus_table", "mandatory", False),

        # Done
        ("reference_phylogeny", "mandatory", False),

        # Done
        ("taxonomy_table", "mandatory", False),

        # this is one load tables that can generated on the fly
        # will need to write the queries/functions
        ("genome_statistics", "mandatory", False),

        ############# Optional ###################
        ("BLAST_database", "optional", False),

        # Done
        ("gene_phylogenies", "optional", False),
        ("interpro_data", "optional", False),
        ("interpro_comparative", "optional", False),
        ("interpro_comparative_accession", "optional", False),
        ("priam_data", "optional", False),
        ("priam_comparative", "optional", False),
        ("priam_comparative_accession", "optional", False),

        # WIP
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

        # Done
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

def setup_chlamdb(**kwargs):
    setup_biodb(kwargs)
    create_data_table(kwargs)

def load_gbk(gbks, args):
    db = db_utils.DB.load_db(args)
    data = []
    for gbk in gbks:
        records = [i for i in SeqIO.parse(gbk, 'genbank')]

        # hack to link the bioentry to the filename, useful later for parsing and
        # storing checkM results in the dtb.
        for record in records:
            db.load_gbk_wrapper([record])
            bioentry_id = db.server.adaptor.last_id("bioentry")
            data.append( (bioentry_id, gbk.replace(".gbk", "")) )

    db.load_filenames(data)
    db.set_status_in_config_table("gbk_files", 1)
    db.commit()

def load_orthofinder_results(orthofinder_output, args):
    db = db_utils.DB.load_db(args)
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
    db = db_utils.DB.load_db(args)
    db.create_refseq_hits_table()

    # map accessions to id
    sseqid_hsh = {}
    sseqid_id = 0
    for tsv in diamond_tsvs:
        hit_table = pd.read_csv(tsv, sep="\t")
        hit_count = 0
        seqfeature_id = None
        query_hash = None
        query_hash_64b = None
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
                query_hash_64b = hsh_from_s(row[0][len("CRC-"):])
                query_hash = row[0]
            lst_args = row.tolist()
            data.append([hit_count, query_hash_64b, sseqid_hsh[match_accession]] + lst_args[2:])
        db.load_refseq_hits(data)
    # see which indices are better fitted
    # db.create_refseq_hits_indices()
    db.commit()
    return sseqid_hsh

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]

def load_refseq_matches_linear_taxonomy(args):
    db = db_utils.DB.load_db(args)
    db.create_refseq_hits_taxonomy()

    taxids = db.get_all_taxids()
    hsh_taxid_to_name = {}
    hsh_taxo_key = db.hsh_taxo_key
    for chunk in chunks(taxids, 2000):
        query_results = db.get_linear_taxonomy(args, chunk)
        results_only_taxids = []
        for result in query_results:
            only_taxids = []
            only_taxids.append(result[0])
            for rank, (idx_name, idx_taxid) in hsh_taxo_key.items():
                if result[idx_taxid] not in hsh_taxid_to_name:
                    hsh_taxid_to_name[result[idx_taxid]] = (rank, result[idx_name]) 
                only_taxids.append(result[idx_taxid])
            results_only_taxids.append(only_taxids)
        db.load_refseq_hits_taxonomy(results_only_taxids)

    db.create_taxonomy_mapping(hsh_taxid_to_name)
    db.create_refseq_hits_taxonomy_indices()
    db.commit()

# hsh_sseqids: maps accession to 
def load_refseq_matches_infos(args, hsh_sseqids):
    db = db_utils.DB.load_db(args)
    db.create_diamond_refseq_match_id()

    # ugly and very slow, to be improved
    for iteration, chunk in enumerate(chunks(list(hsh_sseqids.keys()), 5000)):
        print(iteration, flush=True)
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
    db.create_diamond_refseq_match_id_indices()
    db.commit()

# This is a hack to be able to store 64bit unsigned values into 
# sqlite3's 64 signed value. Values higher than 0x7FFFFFFFFFFFFFFF could not
# be inserted into sqlite3 if unsigned as sqlite3 integers are 64bits signed integer.
def hsh_from_s(s):
    v = int(s, 16)
    if v > 0x7FFFFFFFFFFFFFFF:
        return v-0x10000000000000000
    return v

def load_seq_hashes(args, nr_mapping):
    db = db_utils.DB.load_db(args)
    hsh_locus_to_id = db.get_hsh_locus_to_seqfeature_id()

    to_load_hsh_to_seqid = {}
    for line in open(nr_mapping, "r"):
        record_id, hsh, genome = line.split("\t")
        seqfeature_id = hsh_locus_to_id[record_id]

        short_hsh = hsh[len("CRC-"):]
        int_from_64b_hash = hsh_from_s(short_hsh)
        if int_from_64b_hash not in to_load_hsh_to_seqid:
            to_load_hsh_to_seqid[int_from_64b_hash] = [seqfeature_id]
        else:
            print("Multiple entries for", int_from_64b_hash)
            to_load_hsh_to_seqid[int_from_64b_hash].append(seqfeature_id)

    to_load = []
    for hsh_64b, seqids in to_load_hsh_to_seqid.items():
        for seqid in seqids:
            to_load.append( (hsh_64b, seqid) )

    db.create_seq_hash_to_seqid(to_load)
    db.commit()

def load_alignments_results(args, alignment_files):
    db = db_utils.DB.load_db(args)
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

def load_cog(params, cog_filename):
    db = db_utils.DB.load_db(params)

    # TODO: avoid hardcoding file names. May be worth it to explicitely
    # add them to the nextflow process for more clarity.
    cog2cdd_file = open(params["databases_dir"]+"/COG/cog_corresp.tab", "r")
    cog2length_file = open(params["databases_dir"]+"/COG/cog_length.tab", "r")
    fun_names_file = open(params["databases_dir"]+"/COG/fun2003-2014.tab")
    cog_names_file = open(params["databases_dir"]+"/COG/cognames2003-2014.tab") 

    hsh_cdd_to_cog = {}
    for line in cog2cdd_file:
        tokens = line.split("\t")
        hsh_cdd_to_cog[tokens[1].strip()] = int(tokens[0][3:])

    # maybe an overkill: it seems like COG are ordered
    # it should be possible to only use an array instead of 
    # an hash table
    hsh_cog_to_length = {}
    for line in cog2length_file:
        tokens = line.split("\t")
        hsh_cog_to_length[int(tokens[0][3:])] = int(tokens[1])

    cog_ref_data = []
    # necessary to track, as some cogs listed in the CDD to COG mapping table
    # are not present in those descriptors
    hsh_cog_ids = {} 
    for line_no, line in enumerate(cog_names_file):
        # pass header
        if line_no == 0:
            continue
        tokens = line.split("\t")
        cog_id = int(tokens[0][3:])
        hsh_cog_ids[cog_id] = True
        fun = tokens[1]
        description = tokens[2].strip()
        cog_ref_data.append( (cog_id, fun, description) )
    db.load_cog_ref_data(cog_ref_data)

    cog_fun_data = []
    for line in fun_names_file:
        function, description = line.split("\t") 
        cog_fun_data.append( (function, description) )
    db.load_cog_fun_data(cog_fun_data)
    db.commit()
    
    hsh_seq_to_length = db.seqfeature_to_prot_length()
    cogs_hits = pd.read_csv(cog_filename, sep="\t")
    data = []

    # TODO: check with Trestan whether it is necessary to keep only the 
    # best cog hit for every feature or if all of them should be included in the db.
    # This might mess the cog count in the main page of chlamdb.
    for index, row in cogs_hits.iterrows():
        hsh = hsh_from_s(row[0][len("CRC-"):])
        #  cdd in the form cdd:N
        cog = int(hsh_cdd_to_cog[row[1].split(":")[1]])
        if cog not in hsh_cog_ids:
            print("Unknown cog id: ", cog, " skipping this entry")
            continue
        identity = float(row[2])
        evalue = float(row[10])
        bitscore = float(row[11])
        query_start = int(row[6])
        query_end = int(row[7])
        hit_start = int(row[8])
        hit_end = int(row[9])

        query_coverage = round((query_end-query_start)/hsh_seq_to_length[seqfeature], 2)
        hit_coverage = round((hit_end-hit_start)/hsh_cog_to_length[cog], 2)
        entry = [hsh, cog, query_start, query_end, hit_start]
        entry += [hit_end, query_coverage, hit_coverage, identity]
        entry += [evalue, bitscore]
        data.append(entry)
    db.load_cog_hits(data)
    db.set_status_in_config_table("COG_data", 1)
    db.commit()

# Note: the trees are stored in files with name formatted as:
# OGN_nr_hits_mafft.nwk. To retrieve the orthogroup, parse the filename
# and convert it to int.
#
# Note2: from a database design perspective, may be worth to put all the 
# phylogenies in the same table (BBH/gene and reference) and reference them
# on the orthogroup id and/or a term_id
def load_BBH_phylogenies(kwargs, lst_orthogroups):
    import ete3

    db = db_utils.DB.load_db(kwargs)
    data = []

    for tree in lst_orthogroups:
        t = ete3.Tree(tree)
        og_id = int(tree.split("_")[0][2:])
        data.append( (og_id, t.write()) )
    db.create_BBH_phylogeny_table(data)
    db.set_status_in_config_table("BBH_phylogenies", 1)
    db.commit()

def load_gene_phylogenies(kwargs, lst_orthogroups):
    import ete3

    db = db_utils.DB.load_db(kwargs)
    data = []
    for tree in lst_orthogroups:
        t = ete3.Tree(tree)
        og_id = int(tree.split("_")[0][2:])
        data.append( (og_id, t.write()) )
    db.create_gene_phylogeny_table(data)
    db.set_status_in_config_table("gene_phylogenies", 1)
    db.commit()

def load_reference_phylogeny(kwargs, tree):
    import ete3
    db = db_utils.DB.load_db(kwargs)

    newick_file = open(tree, "r")
    newick_string = newick_file.readline()
    hsh_filename_to_bioentry = db.get_filenames_to_bioentry()

    # convert leaf names to bioentry_id instead of filename
    tree = ete3.Tree(newick_string)
    for leaf in tree.iter_leaves():
        leaf.name = hsh_filename_to_bioentry[leaf.name]

    db.load_reference_phylogeny(tree.write())
    db.set_status_in_config_table("reference_phylogeny", 1)
    db.commit()

# Several values will be inserted in the seqfeature_qualifier_value table, under different
# term_id
# - the GC of the genomes, under the gc term_id
# - the length, under length term_id
# - the number of contigs, under the n_contigs term_id
# - the number of contigs without BBH hits to Chlamydiae
# - the coding density
def load_genomes_summary(kwargs):
    db = db_utils.DB.load_db(kwargs)
    hsh_seq = db.get_genomes_sequences()
    hsh_seq_to_length_coding = db.get_coding_region_total_length()

    gcs = []
    lengths = []
    no_contigs = []
    coding_densities = []
    for entry_id, sequence in hsh_seq.items():
        sequence_without_N = sequence.replace("N", "")

        # During the annotation procedure, contigs are merged, with the insertion
        # of 200 "N"s between them.
        no_contig = sequence.count("N" * 200) + 1
        length = len(sequence_without_N)
        gc = SeqUtils.GC(sequence_without_N)
        coding_density = round(100*hsh_seq_to_length_coding[entry_id]/length, 2)

        gcs.append( (entry_id, gc) )
        lengths.append( (entry_id, length) )
        no_contigs.append( (entry_id, no_contig))
        coding_densities.append( (entry_id, coding_density) )
    db.load_genomes_gc(gcs)
    db.load_genomes_lengths(lengths)
    db.load_genomes_n_contigs(no_contigs)
    db.load_genomes_coding_density(coding_densities)
    db.commit()

def load_checkm_results(params, checkm_results):
    db = db_utils.DB.load_db(params)
    tab = pd.read_table(checkm_results)

    hsh_filename_to_bioentry = db.get_filenames_to_bioentry()
    data = []
    for index, row in tab.iterrows():
        values = (hsh_filename_to_bioentry[row["Bin Id"]], row["Completeness"], row["Contamination"])
        data.append(values)

    db.load_checkm_results(data)
    db.commit()
