import os
import sys
import re

from metagenlab_libs import db_utils
from metagenlab_libs.chlamdb import search_bar

import pandas as pd

from Bio import SeqIO
from Bio import AlignIO
from Bio import SeqUtils

from collections import namedtuple
from collections import defaultdict


# assumes orthofinder named: OG000N
# returns the N as int
def get_og_id(string):
    return int(string[2:])

def get_ko_id(string):
    return int(string[1:])

def parse_orthofinder_output_file(output_file):
    protein_id2orthogroup_id = {}
    parsing = open(output_file, 'r')

    for line in parsing:
        tokens = line.strip().split(' ')

        # Skips the ":" at the end of the orthgroup id
        group = get_og_id(tokens[0][:-1])
        for locus in tokens[1:]:
            assert locus not in protein_id2orthogroup_id
            protein_id2orthogroup_id[locus] = group
    return protein_id2orthogroup_id


def is_plasmid(record):
    # XXX: after discussion with Trestan, chose to keep it
    # as it allows to catch malformated genbanks that do not
    # specify plasmids correctly. An alternative way would be 
    # to just ignore it.
    if "plasmid" in record.description:
        return True

    for fet in record.features:
        if fet.type == "source" and "plasmid" in fet.qualifiers:
            return True
    return False


def load_gbk(gbks, args, db_file):
    """
    NOTE:

    This function assumes that only one organism is present in each 
    genbank file. An arbitraty taxon_id is assigned to each genbank file
    to differentiate the organisms in chlamdb. It has nothing to do
    with the taxon id from NCBI. It just is an identifier to group
    all the bioentries from a bacteria together (its chromosomes and plasmids).
    """
    db = db_utils.DB.load_db(db_file, args)
    data = []
     
    bioentry_plasmids = []
    bioentry_to_taxid = []
    for taxon_id, gbk in enumerate(gbks):
        records = [i for i in SeqIO.parse(gbk, 'genbank')]

        for record in records:
            db.load_gbk_wrapper([record])
            bioentry_id = db.server.adaptor.last_id("bioentry")
            bioentry_plasmids.append( (bioentry_id, is_plasmid(record)))
            bioentry_to_taxid.append( (bioentry_id, taxon_id+1) )

        # hack to link the bioentry to the filename, useful later for parsing and
        # storing checkM results in the dtb.
        data.append( (taxon_id+1, gbk.replace(".gbk", "")) )

    db.load_filenames(data)
    db.commit()
    db.update_taxon_ids(bioentry_to_taxid)
    db.update_plasmid_status(bioentry_plasmids)
    db.set_status_in_config_table("gbk_files", 1)
    db.commit()


def load_orthofinder_results(orthofinder_output, args, db_file):
    db = db_utils.DB.load_db(db_file, args)
    hsh_prot_to_group = parse_orthofinder_output_file(orthofinder_output)
    hsh_locus_to_feature_id = db.get_hsh_locus_to_seqfeature_id(only_CDS=True)
    hits_to_load = [(hsh_locus_to_feature_id[locus], group) for locus, group in hsh_prot_to_group.items()]
    db.load_og_hits(hits_to_load)
    db.set_status_in_config_table("orthology", 1)
    db.commit()


# Note: as this is an alignment, the lengths are the same
# to be replaced by zip and to be tested
def get_identity(seq1, seq2):
    if len(seq1)!=len(seq2):
        raise RuntimeException("The lengths of two aligned sequences should be identical")

    identity = 0
    aligned = 0
    identical = 0
    gaps_1 = 0
    gaps_2 = 0

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


# will need to rewrite it using iterators instead of 
# copying the whole list
def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]


def get_prot(refseq_file, hsh_accession):
    # Brute-force scanning of the whole refseq file.
    # This is however faster than serial SQL queries.

    hsh_results = {}
    for record in SeqIO.parse(refseq_file, "fasta"):
        accession = remove_accession_version(record.name)
        if not accession in hsh_accession:
            continue
        hsh_results[accession] = record
    return hsh_results


def get_taxids(ncbi_tax_file, hsh_accession_to_taxid):
    line_iter = open(ncbi_tax_file, "r")
    
    # pass header
    next(line_iter)
    for line in open(ncbi_tax_file, "r"):
        end_of_accession = line.index("\t")
        accession = line[:end_of_accession]
        if not accession in hsh_accession_to_taxid:
            continue
        tokens = line.split("\t")
        hsh_accession_to_taxid[accession] = int(tokens[2])


def remove_accession_version(accession):
    return accession.split(".")[0]


def simplify_hash(hsh):
    return hsh_from_s(hsh[len("CRC-"):])


def parse_description(description):
    end_accession_index = description.index(' ')
    # NOTE: this assumes that no brackets are found in the
    # name of the protein
    beg_tax_index = description.index('[')
    # also skip the white spaces
    return description[end_accession_index+1:beg_tax_index-1]


# Heavy-duty function, with several things performed to avoid database
# queries and running too many times acrosse refseq_merged.faa :
#  * get the mapping between accession to taxid
#  * get the linear taxonomy for all taxids that were extracted at the previous step
#  * extract the protein informations for 
#  * for all protein of all orthogroup, get the non-PVC hits and output them in OG.fasta
#
# This function will not scale for big sets of genomes: it will be necessary
# to run it several time for a given chunk of results.
#
# Alternatively, already pre-process the results directly after diamond finishes. This should
# be possible as the files are accessed as read-only.
def load_refseq_matches_infos(args, lst_diamond_files, db_file):
    db = db_utils.DB.load_db(db_file, args)
    columns=["str_hsh", "accession", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

    print("Reading tsvs", flush=True)
    all_data = pd.DataFrame(columns=columns)
    for tsv in lst_diamond_files:
        hit_table = pd.read_csv(tsv, sep="\t", names=columns, header=None)
        all_data = all_data.append(hit_table)

    print("Initial size: ", len(all_data.index), flush=True)
    all_data.accession = all_data.accession.map(remove_accession_version)
    all_data.str_hsh = all_data.str_hsh.map(simplify_hash)
    hsh_accession_to_taxid = {accession:None for index, accession in all_data.accession.iteritems()}

    print("Obtaining accession to taxid", flush=True)
    get_taxids(args["ncbi_taxonomy"], hsh_accession_to_taxid)
    all_data["taxid"] = all_data.accession.map(hsh_accession_to_taxid)

    taxid_set = set(hsh_accession_to_taxid.values())
    non_pvc_taxids = set()
    pvc = args.get("refseq_diamond_BBH_phylogeny_phylum_filter", [])
    to_avoid = set()
    hsh_taxid_to_name = {}
    results_only_taxids = []

    print("Getting linear taxonomy", flush=True)
    for chunk in chunks(list(taxid_set), 5000):
        query_results = db.get_linear_taxonomy(args, chunk)
        for result in query_results:
            # Not necessary to keep those in the database
            if result.phylum() in pvc:
                continue
            non_pvc_taxids.add(result.taxid())
            results_only_taxids.append(result.get_all_taxids())
            result.update_hash(hsh_taxid_to_name)

    db.create_refseq_hits_taxonomy()
    db.load_refseq_hits_taxonomy(results_only_taxids)
    db.create_taxonomy_mapping(hsh_taxid_to_name)
    db.create_refseq_hits_taxonomy_indices()

    refseq = args["databases_dir"] + "/refseq/merged.faa"
    print("Extracting records for refseq", flush=True)
    hsh_accession_to_record = get_prot(refseq, hsh_accession_to_taxid)
    hsh_accession_to_match_id = {}

    db.create_diamond_refseq_match_id()
    refseq_match_id = []
    print("Loading refseq matches id", flush=True)

    for sseqid, (accession, taxid) in enumerate(hsh_accession_to_taxid.items()):
        hsh_accession_to_match_id[accession] = sseqid
        record = hsh_accession_to_record[accession]
        data = [sseqid, accession, taxid, parse_description(record.description), len(record)]
        refseq_match_id.append(data)
    db.load_diamond_refseq_match_id(refseq_match_id)
    db.create_diamond_refseq_match_id_indices()

    print("Loading refseq matches", flush=True)
    db.create_refseq_hits_table()
    all_data = all_data[all_data.taxid.isin(non_pvc_taxids)]
    all_data["count"] = all_data.groupby("str_hsh").cumcount()
    all_data["match_id"] = all_data["accession"].map(hsh_accession_to_match_id)

    print("after filtering: ", len(all_data), flush=True)
    to_load = all_data[["count", "str_hsh", "match_id", "pident",
        "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
        "send", "evalue", "bitscore"]]
    db.load_refseq_hits(to_load.values.tolist())
    db.create_refseq_hits_indices()

    if args.get("refseq_diamond_BBH_phylogeny", True):
        max_hits = args.get("refseq_diamond_BBH_phylogeny_top_n_hits", 4)
        all_og = db.get_all_orthogroups(min_size=3)
        for og in all_og:
            refseq_matches = db.get_diamond_match_for_og(og)
            to_keep = []
            cur_accesion, cur_count = None, 0
            for accession, taxid in refseq_matches:
                if taxid not in non_pvc_taxids:
                    continue
                if accession != cur_accesion:
                    cur_accesion = accession
                    cur_count = 0
                if cur_count == max_hits:
                    continue
                to_keep.append(hsh_accession_to_record[accession])
                cur_count += 1
            sequences = db.get_all_sequences_for_orthogroup(og)
            SeqIO.write(to_keep + sequences, f"{og}_nr_hits.faa", "fasta")
    else:
        # just create a empty file to avoid a tantrum of nextflow for a missing file
        f = open("null_nr_hits.faa", "w")
        f.write("My hovercraft is full of eels")
        f.close()
    db.set_status_in_config_table("BLAST_database", 1)
    db.commit()

# This is a hack to be able to store 64bit unsigned values into 
# sqlite3's 64 signed value. Values higher than 0x7FFFFFFFFFFFFFFF could not
# be inserted into sqlite3 if unsigned as sqlite3 integers are 64bits signed integer.
def hsh_from_s(s):
    v = int(s, 16)
    if v > 0x7FFFFFFFFFFFFFFF:
        return v-0x10000000000000000
    return v

def load_seq_hashes(args, nr_mapping, db_file):
    db = db_utils.DB.load_db(db_file, args)
    hsh_locus_to_id = db.get_hsh_locus_to_seqfeature_id(only_CDS=True)

    to_load_hsh_to_seqid = {}
    for line in open(nr_mapping, "r"):
        record_id, hsh, genome = line.split("\t")
        seqfeature_id = hsh_locus_to_id[record_id]

        short_hsh = hsh[len("CRC-"):]
        int_from_64b_hash = hsh_from_s(short_hsh)
        if int_from_64b_hash not in to_load_hsh_to_seqid:
            to_load_hsh_to_seqid[int_from_64b_hash] = [seqfeature_id]
        else:
            to_load_hsh_to_seqid[int_from_64b_hash].append(seqfeature_id)

    to_load = []
    for hsh_64b, seqids in to_load_hsh_to_seqid.items():
        for seqid in seqids:
            to_load.append( (hsh_64b, seqid) )

    db.create_seq_hash_to_seqid(to_load)
    db.commit()


def load_alignments_results(args, alignment_files, db_file):
    db = db_utils.DB.load_db(db_file, args)
    db.create_new_og_matrix()
    locus_to_feature_id = db.get_hsh_locus_to_seqfeature_id(only_CDS=True)

    # assumes filename of the format OG00N_mafft.faa, with the orthogroup
    # being the integer following the OG string
    matrix = []
    for alignment in alignment_files:
        align = AlignIO.read(alignment, "fasta")
        orthogroup = get_og_id(alignment.split("_")[0])
        for i in range(len(align)):
            for j in range(i+1, len(align)):
                alignment_1 = align[i]
                alignment_2 = align[j]
                id_1 = locus_to_feature_id[alignment_1.name]
                id_2 = locus_to_feature_id[alignment_2.name]
                identity = get_identity(alignment_1, alignment_2)
                matrix.append( (orthogroup, id_1, id_2, identity) )
    db.load_og_matrix(matrix)
    db.create_og_matrix_indices()
    db.set_status_in_config_table("orthogroup_alignments", 1)
    db.commit()


def load_cog(params, filelist, db_file):
    db = db_utils.DB.load_db(db_file, params)
    hsh_cdd_to_cog = db.get_cdd_to_cog()

    data = []
    for chunk in filelist:
        cogs_hits = pd.read_csv(chunk, sep="\t", header=None,
            names=["seq_hsh", "cdd", "pident", "length", "mismatch", "gapopen", "qstart",
                "qend", "sstart", "send", "evalue", "bitscore"])

        # Select only the best hits: using pandas clearly is an overkill here
        min_hits = cogs_hits.groupby("seq_hsh")[["cdd", "evalue"]].min()
        for index, row in min_hits.iterrows():
            hsh = hsh_from_s(index[len("CRC-"):])
            #  cdd in the form cdd:N
            cog = hsh_cdd_to_cog[int(row["cdd"].split(":")[1])]
            evalue = float(row["evalue"])
            entry = [hsh, cog, evalue]
            data.append(entry)
    db.load_cog_hits(data)
    db.set_status_in_config_table("COG", 1)
    db.commit()


# Note: the trees are stored in files with name formatted as:
# OGN_nr_hits_mafft.nwk. To retrieve the orthogroup, parse the filename
# and convert it to int.
#
# Note2: from a database design perspective, may be worth to put all the 
# phylogenies in the same table (BBH/gene and reference) and reference them
# on the orthogroup id and/or a term_id
def load_BBH_phylogenies(kwargs, lst_orthogroups, db_file):
    import ete3

    db = db_utils.DB.load_db(db_file, kwargs)
    data = []

    for tree in lst_orthogroups:
        t = ete3.Tree(tree)
        og_id = int(tree.split("_")[0])
        data.append( (og_id, t.write()) )
    db.create_BBH_phylogeny_table(data)
    db.set_status_in_config_table("BBH_phylogenies", 1)
    db.commit()


def load_gene_phylogenies(kwargs, lst_orthogroups, db_file):
    """
    NOTE: the leafs of those tree are locus tags
    """
    import ete3

    db = db_utils.DB.load_db(db_file, kwargs)
    data = []
    for tree in lst_orthogroups:
        t = ete3.Tree(tree)
        og_id = int(tree.split("_")[0][2:])
        data.append( (og_id, t.write()) )
    db.create_gene_phylogeny_table(data)
    db.set_status_in_config_table("gene_phylogenies", 1)
    db.commit()


def load_reference_phylogeny(kwargs, tree, db_file):
    import ete3
    db = db_utils.DB.load_db(db_file, kwargs)

    newick_file = open(tree, "r")
    newick_string = newick_file.readline()
    hsh_filename_to_taxon = db.get_filenames_to_taxon_id()

    # convert leaf names to taxon_id instead of filename
    tree = ete3.Tree(newick_string)
    for leaf in tree.iter_leaves():
        leaf.name = hsh_filename_to_taxon[leaf.name]

    db.load_reference_phylogeny(tree.write())
    db.set_status_in_config_table("reference_phylogeny", 1)
    db.commit()


def get_gen_stats(gbk_list):
    hsh_gen_stats = {}

    for gbk_file in gbk_list:
        ttl_length = 0
        gc_cum = 0
        cds_length = 0
        for record in SeqIO.parse(gbk_file, "genbank"):
            ttl_length += len(record)
            gc_cum += SeqUtils.GC(record.seq)*len(record)
            for fet in record.features:
                if fet.type=="CDS":
                    location = fet.location
                    cds_length += location.end-location.start
        gbk_shortened = gbk_file.replace(".gbk", "")
        hsh_gen_stats[gbk_shortened] = (float(gc_cum)/ttl_length, float(cds_length)/ttl_length, ttl_length)
    return hsh_gen_stats


def load_genomes_info(kwargs, gbk_list, checkm_results, db_file):
    db = db_utils.DB.load_db(db_file, kwargs)
    tab = pd.read_table(checkm_results)

    hsh_taxid_to_gen_stats = get_gen_stats(gbk_list)

    hsh_filename_to_taxid = db.get_filenames_to_taxon_id()
    data = []
    for index, row in tab.iterrows():
        taxon_id = hsh_filename_to_taxid[row["Bin Id"]]
        completeness = row.Completeness
        contamination = row.Contamination
        gc, coding_density, length = hsh_taxid_to_gen_stats[row["Bin Id"]]
        values = [taxon_id, completeness, contamination, gc, length, coding_density]
        data.append(values)
        
    db.load_genomes_info(data)
    db.set_status_in_config_table("genome_statistics", 1)
    db.commit()

PfamEntry = namedtuple("PfamEntry", ["accession", "description"])


def parse_pfam_entry(file_iter):
    accession, description = None, None 
    for line in file_iter:
        if line.startswith("//"):
            yield PfamEntry(accession=accession, description=description)
            accession, description = None, None
        if line.startswith("#=GF AC"):
            accession_offset = len("#=GF AC   PF") 
            accession_length = 5
            accession = int(line[accession_offset:accession_offset+accession_length])
        elif line.startswith("#=GF DE"):
            description = line[len("#=GF DE   "):-1]


def load_pfam(params, pfam_files, db, pfam_def_file):
    db = db_utils.DB.load_db(db, params)

    db.create_pfam_hits_table()
    pfam_ids = set()
    for pfam in pfam_files:
        entries = []
        for line in open(pfam, "r"):
            if len(line)<len("CRC-") or line[0:len("CRC-")]!="CRC-":
                continue

            tokens = line.split()
            hsh_i = simplify_hash(tokens[0])
            start = int(tokens[1])
            end = int(tokens[2])
            pfam_raw_str = tokens[5].split(".")[0]
            pfam_i = int(pfam_raw_str[len("PF"):])
            entries.append((hsh_i, pfam_i, start, end))
            pfam_ids.add(pfam_i)

        db.load_data_into_table("pfam_hits", entries)
    db.commit()

    pfam_entries = []
    pfam_def_file_iter = open(pfam_def_file, "r")
    for entry in parse_pfam_entry(pfam_def_file_iter):
        if not entry.accession in pfam_ids:
            continue
        pfam_entries.append([entry.accession, entry.description])
    db.create_pfam_def_table(pfam_entries)
    db.set_status_in_config_table("pfam", 1)
    db.commit()


class SwissProtIdCount(defaultdict):
    def __init__(self):
        super().__init__()
        self.id_count = 0

    def __missing__(self, key):
        tmp = self.id_count
        self[key] = self.id_count
        self.id_count += 1
        return tmp

def parse_swissprot_entry(description):
    tokens = iter(description.split())

    # skip entry id
    next(tokens)
    acc = []
    gene = None
    for curr_token in tokens:
        if curr_token.startswith("OS="):
            organism = curr_token[len("OS="):]
            prot_name = " ".join(acc)
            acc = [organism]
        elif curr_token.startswith("OX="):
            taxid = int(curr_token[len("OX="):])
        elif curr_token.startswith("GN="):
            gene = curr_token[len("GN="):]
        elif curr_token.startswith("PE="):
            pe = int(curr_token[len("PE="):])
        elif curr_token.startswith("SV="):
            version = curr_token[len("SV="):]
        else:
            acc.append(curr_token)
    return prot_name, taxid, " ".join(acc), gene, pe, version


# Swissprot id are in the format
# db|UniqueIdentifier|EntryName
def parse_swissprot_id(to_parse):
    db, ident, name = to_parse.split("|")
    return db, ident, name


def load_swissprot(params, blast_results, db_name, swissprot_fasta):
    db = db_utils.DB.load_db(db_name, params)
    hsh_swissprot_id = SwissProtIdCount()
    db.create_swissprot_tables()

    # Note: this is not really scalable to x genomes, as 
    # it necessitates to keep the prot id in memory
    # may need to process the blast results in batch instead (slower but
    # spares memory).
    for blast_file in blast_results:
        data = []
        for line in open(blast_file, "r"):
            crc, prot_id, perid, leng, n_mis, n_gap, qs, qe, ss, se, e, score = line.split()
            hsh = simplify_hash(crc)
            # swissprot accession in the format x|prot_id|org
            _, prot_id, _ = parse_swissprot_id(prot_id)
            db_prot_id = hsh_swissprot_id[prot_id]
            data.append((hsh, db_prot_id, float(e), int(float(score)),
                int(float(perid)), int(n_gap), int(leng)))
        db.load_swissprot_hits(data)

    swiss_prot_defs = []
    for record in SeqIO.parse(swissprot_fasta, "fasta"):
        _, prot_id, _ = parse_swissprot_id(record.name)
        if prot_id not in hsh_swissprot_id:
            continue
        db_prot_id = hsh_swissprot_id[prot_id]
        descr, taxid, org, gene, pe, version = parse_swissprot_entry(record.description)
        swiss_prot_defs.append((db_prot_id, prot_id, descr, taxid, org, gene, pe))
    db.load_swissprot_defs(swiss_prot_defs)
    db.set_status_in_config_table("BLAST_swissprot", 1)
    db.commit()


def simplify_ko(raw_ko):
    return int(raw_ko[len("ko:K"):])


# NOTE:
# Several KO marked as significant can be assigned to the same locus
# only take the hit with the lowest evalue (the first in the list)
def load_KO(params, ko_files, db_name):
    db = db_utils.DB.load_db(db_name, params)
    data = []
    for ko_file in ko_files:
        curr_hsh = None
        for ko_line in open(ko_file, "r"):
            tokens = ko_line.split()
            # ignore all but the best hits
            if tokens[0] != "*":
                continue
            crc_raw, ko_str, thrs_str, score_str, evalue_str, *descr = tokens[1:]
            hsh = simplify_hash(crc_raw)
            if hsh == curr_hsh:
                # skip the entries that were classified as significant, but
                # with a higher e-value
                continue
            else:
                curr_hsh = hsh
            ko = get_ko_id(ko_str)
            thrs = float(thrs_str)
            score = float(score_str)
            evalue = float(evalue_str)
            entry = [hsh, ko, thrs, score, evalue]
            data.append(entry)
    db.load_ko_hits(data)
    db.set_status_in_config_table("KEGG", 1)
    db.commit()


def load_blastDBs(params, db_name):
    db = db_utils.DB.load_db(db_name, params)
    db.set_status_in_config_table("BLAST_database", 1)
    db.commit()
    

def format_cog(cog_n):
    if pd.isna(cog_n):
        return None
    return f"COG{int(cog_n):04}"

def format_ko(ko_n):
    if pd.isna(ko_n):
        return None
    return f"K{int(ko_n):05}"

def format_og(og_n):
    return f"group_{int(og_n)}"

def format_pfam(pfam_n):
    if pd.isna(pfam_n):
        return None
    return f"PF{int(pfam_n):05}"


def setup_chlamdb_search_index(params, db_name, index_name):
    db = db_utils.DB.load_db(db_name, params)
    os.mkdir(index_name)

    has_cog = params.get("cog", False)
    has_ko = params.get("ko", False)
    has_pfam = params.get("pfam_scan", False)

    genomes = db.get_genomes_description()
    index = search_bar.ChlamdbIndex.new_index(index_name)

    for taxid, data in genomes.iterrows():
        all_infos = db.get_proteins_info([taxid], search_on="taxid",
                as_df=True, inc_non_CDS=True, inc_pseudo=True)
        for seqid, data in all_infos.iterrows():
            product, locus_tag, gene = data[["product", "locus_tag", "gene"]]
            if pd.isna(gene):
                gene = None
            if product=="hypothetical protein":
                product = None

            organism = genomes.loc[taxid].description
            index.add(locus_tag=locus_tag,
                    name=gene, description=product, organism=organism, 
                    entry_type = search_bar.EntryTypes.Gene)

    if has_cog:
        cog_data = db.get_cog_summaries(cog_ids=None, only_cog_desc=True)
        for cog, (func, descr) in cog_data.items():
            index.add(name=format_cog(cog),
                    description=descr, entry_type=search_bar.EntryTypes.COG)

    if has_ko:
        ko_data = db.get_ko_desc(ko_ids=None)
        for ko, descr in ko_data.items():
            index.add(name=format_ko(ko), description=descr,
                    entry_type=search_bar.EntryTypes.KO)

    if has_pfam:
        pfam_data = db.get_pfam_def(pfam_ids=None)
        for pfam, data in pfam_data.iterrows():
            index.add(name=format_pfam(pfam),
                    description=data["def"],
                    entry_type=search_bar.EntryTypes.PFAM)
    index.done_adding()

