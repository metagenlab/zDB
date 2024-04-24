import os
import re
from collections import defaultdict, namedtuple
from datetime import datetime

import pandas as pd
import xlrd
from Bio import SeqIO, SeqUtils
from lib import KO_module, search_bar
from lib.db_utils import DB

from annotations import InputHandler


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
    parsing.close()
    return protein_id2orthogroup_id


def is_plasmid(record):
    # it allows to catch malformated genbanks that do not
    # specify plasmids correctly. An alternative way would be
    # to just ignore it.
    if "plasmid" in record.description:
        return True

    for fet in record.features:
        if fet.type == "source" and "plasmid" in fet.qualifiers:
            return True
    return False


def load_gbk(gbks, args, db_file):
    db = DB.load_db(db_file, args)
    data = []

    bioentry_plasmids = []
    for gbk in gbks:
        records = [i for i in SeqIO.parse(gbk, 'genbank')]
        taxon_id = None

        for record in records:
            db.load_gbk_wrapper([record])
            bioentry_id = db.server.adaptor.last_id("bioentry")

            # this assumes that the taxon was not already in the database
            # as long as we enforce the "unique organism per file" in check_gbk,
            # this should work
            taxon_id = db.server.adaptor.last_id("taxon")
            bioentry_plasmids.append((bioentry_id, is_plasmid(record)))

        # hack to link the bioentry to the filename, useful later for parsing and
        # storing checkM results in the dtb.
        data.append((taxon_id, os.path.splitext(gbk)[0]))

    db.load_filenames(data)
    db.commit()
    db.update_plasmid_status(bioentry_plasmids)
    db.set_status_in_config_table("gbk_files", 1)
    db.commit()


def load_groups(input_file, kwargs, db_file):
    db = DB.load_db(db_file, kwargs)
    csv_entries, group_names = InputHandler.parse_csv(input_file)
    filenames_to_taxid = db.get_filenames_to_taxon_id()
    group_taxon = []
    for entry in csv_entries:
        taxon_id = filenames_to_taxid[os.path.splitext(entry.file)[0]]
        for group in entry.groups:
            group_taxon.append((group, taxon_id))

    db.load_groups([[name] for name in group_names.values()], group_taxon)
    db.commit()


def load_orthofinder_results(orthofinder_output, args, db_file):
    db = DB.load_db(db_file, args)
    hsh_prot_to_group = parse_orthofinder_output_file(orthofinder_output)
    hsh_locus_to_feature_id = db.get_hsh_locus_to_seqfeature_id(only_CDS=True)
    hits_to_load = [(hsh_locus_to_feature_id[locus], group)
                    for locus, group in hsh_prot_to_group.items()]
    db.load_og_hits(hits_to_load)
    db.set_status_in_config_table("orthology", 1)
    db.commit()


# will need to rewrite it using iterators instead of
# copying the whole list
def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def get_prot(refseq_file, hsh_accession):
    # Brute-force scanning of the whole refseq file.
    # This is however faster than serial SQL queries.

    for record in SeqIO.parse(refseq_file, "fasta"):
        accession = remove_accession_version(record.name)
        if accession not in hsh_accession:
            continue
        hsh_accession[accession] = record


def remove_accession_version(accession):
    return accession.split(".")[0]


def simplify_hash(hsh):
    return hsh_from_s(hsh[len("CRC-"):])


def parse_record(record):
    description = record.description
    end_accession_index = description.index(' ')
    # NOTE: this assumes that no brackets are found in the
    # name of the protein
    beg_tax_index = description.index('[')
    end_tax_index = description.index(']')
    # also skip the white spaces

    prot_descr = description[end_accession_index + 1:beg_tax_index - 1]
    organism = description[beg_tax_index + 1:end_tax_index]
    return prot_descr, organism


def load_refseq_matches_infos(args, lst_diamond_files, db_file):
    db = DB.load_db(db_file, args)
    columns = ["str_hsh", "accession", "pident", "length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

    print("Reading tsvs", flush=True)
    all_data = pd.DataFrame(columns=columns)
    for tsv in lst_diamond_files:
        hit_table = pd.read_csv(tsv, sep="\t", names=columns, header=None)
        all_data = all_data.append(hit_table)
    all_data.accession = all_data.accession.map(remove_accession_version)
    all_data.str_hsh = all_data.str_hsh.map(simplify_hash)

    hsh_accession_to_record = {
        accesion: None for accesion in all_data.accession.tolist()}

    print("Extracting records for refseq", flush=True)
    refseq = args["refseq_db"] + "/merged.faa"
    get_prot(refseq, hsh_accession_to_record)
    hsh_accession_to_match_id = {}

    db.create_diamond_refseq_match_id()
    refseq_match_id = []

    print("Loading refseq matches id", flush=True)
    for sseqid, (accession, record) in enumerate(hsh_accession_to_record.items()):
        hsh_accession_to_match_id[accession] = sseqid
        record = hsh_accession_to_record[accession]
        descr, organism = parse_record(record)
        data = [sseqid, accession, organism, descr, len(record)]
        refseq_match_id.append(data)
    db.load_diamond_refseq_match_id(refseq_match_id)
    db.create_diamond_refseq_match_id_indices()

    print("Loading refseq matches", flush=True)
    db.create_refseq_hits_table()
    all_data["count"] = all_data.groupby("str_hsh").cumcount()
    all_data["match_id"] = all_data["accession"].map(hsh_accession_to_match_id)
    to_load = all_data[["count", "str_hsh", "match_id", "pident",
                        "length", "mismatch", "gapopen", "qstart", "qend", "sstart",
                        "send", "evalue", "bitscore"]]
    db.load_refseq_hits(to_load.values.tolist())
    db.create_refseq_hits_indices()
    db.commit()

    max_hits = args.get("refseq_diamond_BBH_phylogeny_top_n_hits", 100)
    all_og = db.get_all_orthogroups()
    for og, og_size in all_og:
        if og_size < 3:
            continue

        refseq_matches = db.get_diamond_match_for_og(og, sort_by_evalue=True)
        sequences = db.get_all_sequences_for_orthogroup(og)
        locus_set = {seq.id for seq in sequences}
        n_hits = 0
        to_keep = []
        for idx, data in refseq_matches.iterrows():
            if data.accession in locus_set:
                # if the genome was downloaded, it is likely that the same
                # protein will be present in both og and refseq hits
                continue
            to_keep.append(hsh_accession_to_record[data.accession])
            n_hits += 1
            if n_hits == max_hits:
                break

        SeqIO.write(to_keep + sequences, f"{og}_nr_hits.faa", "fasta")
    db.set_status_in_config_table("BLAST_database", 1)
    db.commit()


# This is a hack to be able to store 64bit unsigned values into
# sqlite3's 64 signed value. Values higher than 0x7FFFFFFFFFFFFFFF could not
# be inserted into sqlite3 if unsigned as sqlite3 integers are 64bits signed integer.
def hsh_from_s(s):
    v = int(s, 16)
    if v > 0x7FFFFFFFFFFFFFFF:
        return v - 0x10000000000000000
    return v


def load_seq_hashes(args, nr_mapping, db_file):
    db = DB.load_db(db_file, args)
    hsh_locus_to_id = db.get_hsh_locus_to_seqfeature_id(only_CDS=True)

    to_load_hsh_to_seqid = {}
    with open(nr_mapping, "r") as nr_file:
        for line in nr_file:
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
            to_load.append((hsh_64b, seqid))

    db.create_seq_hash_to_seqid(to_load)
    db.commit()


def load_alignments_results(args, identity_csvs, db_file):
    db = DB.load_db(db_file, args)
    db.create_new_og_matrix()
    locus_to_feature_id = db.get_hsh_locus_to_seqfeature_id(only_CDS=True)

    # assumes filename of the format OG00N_mafft.faa, with the orthogroup
    # being the integer following the OG string
    for identity_csv in identity_csvs:
        orthogroup = get_og_id(identity_csv.split("_")[0])
        matrix = []
        with open(identity_csv, "r") as csv_file:
            for line in csv_file:
                lt1, lt2, ident, le = line.split(",")
                id1 = locus_to_feature_id[lt1]
                id2 = locus_to_feature_id[lt2]
                matrix.append((orthogroup, id1, id2, float(ident)))
        db.load_og_matrix(matrix)
    db.create_og_matrix_indices()
    db.set_status_in_config_table("orthogroup_alignments", 1)
    db.commit()


def load_cog(params, filelist, db_file, cdd_to_cog, cog_db_dir):
    db = DB.load_db(db_file, params)
    hsh_cdd_to_cog = {}
    with open(cdd_to_cog, "r") as cdd_to_cog_file:
        for line in cdd_to_cog_file:
            cog, cdd = line.split()
            hsh_cdd_to_cog[int(cdd)] = int(cog)

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

    # Determine reference DB version
    with open(os.path.join(cog_db_dir, "cdd.info")) as fh:
        dbname, _, version = fh.readline().split()
        if dbname != "cdd":
            raise ValueError("Could not determine Cog DB version")
    db.load_data_into_table("versions", [("CDD", version.strip())])
    db.commit()


def amr_hit_to_db_entry(hit):
    columns = ["Gene symbol",
               "Sequence name",
               "Scope",
               "Element type",
               "Element subtype",
               "Class",
               "Subclass",
               "% Coverage of reference sequence",
               "% Identity to reference sequence",
               "Accession of closest sequence",
               "Name of closest sequence",
               "HMM id"]
    entry = [hsh_from_s(hit["Protein identifier"][len("CRC-"):])]
    entry.extend([hit[column] for column in columns])
    return entry


def load_amr(params, filelist, db_file, version_file):
    db = DB.load_db(db_file, params)

    data = []
    for chunk in filelist:
        amr_hits = pd.read_csv(chunk, sep="\t", header=0)
        data.extend(amr_hit_to_db_entry(hit) for i, hit in amr_hits.iterrows())
    db.load_amr_hits(data)
    db.set_status_in_config_table("AMR", 1)
    db.commit()

    # Determine software and reference DB version
    with open(version_file) as fh:
        for line in fh:
            if line.startswith("Software version"):
                soft_version = line.rsplit(":", 1)[1].strip()
            elif line.startswith("Database version"):
                db_version = line.rsplit(":", 1)[1].strip()
    db.load_data_into_table("versions", [("AMRFinderSoftware", soft_version),
                                         ("AMRFinderDB", db_version)])
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

    db = DB.load_db(db_file, kwargs)
    data = []

    for tree in lst_orthogroups:
        t = ete3.Tree(tree)
        og_id = int(tree.split("_")[0])
        data.append((og_id, t.write()))
    db.create_BBH_phylogeny_table(data)
    db.set_status_in_config_table("BBH_phylogenies", 1)
    db.commit()


def load_gene_phylogenies(kwargs, og_summary, lst_orthogroups, db_file):
    """
    NOTE: the leafs of those tree are locus tags
    """
    import ete3

    db = DB.load_db(db_file, kwargs)
    hsh_ogs = {}
    for tree in lst_orthogroups:
        t = ete3.Tree(tree)
        og_id = int(tree.split("_")[0][2:])
        hsh_ogs[og_id] = t.write()

    data = []
    with open(og_summary, "r") as og_file:
        for line in og_file:
            og, is_core, og_size, num_genomes = line.split("\t")
            og_id = int(og.split("_")[0][2:])
            newick = hsh_ogs.get(og_id, "")
            data.append((og_id, newick, int(is_core),
                        int(og_size), int(num_genomes)))

    db.create_gene_phylogeny_table(data)
    db.set_status_in_config_table("gene_phylogenies", 1)
    db.commit()


def load_reference_phylogeny(kwargs, tree, db_file):
    import ete3
    db = DB.load_db(db_file, kwargs)

    with open(tree, "r") as newick_file:
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
    # NOTE: for now, the coding density do not take overlapping genes
    # into account. Depending on how many of them are present in a genome,
    # this may cause an overestimation of the coding density, as each
    # CDS will be accounted for separately (and a same region will be counted
    # several times).

    hsh_gen_stats = {}

    for gbk_file in gbk_list:
        ttl_length = 0
        gc_cum = 0
        cds_length = 0
        for record in SeqIO.parse(gbk_file, "genbank"):
            ttl_length += len(record)
            gc_cum += SeqUtils.gc_fraction(record.seq) * len(record)
            for fet in record.features:
                if fet.type in ["CDS", "tmRNA", "rRNA", "ncRNA", "tRNA"]:
                    if "pseudo" in fet.qualifiers:
                        continue
                    location = fet.location

                    # allow to take compoundlocation into account
                    for part in location.parts:
                        cds_length += part.end - part.start
        gbk_shortened = os.path.splitext(gbk_file)[0]
        hsh_gen_stats[gbk_shortened] = (float(gc_cum) / ttl_length,
                                        float(cds_length) / ttl_length, ttl_length)
    return hsh_gen_stats


def load_genomes_info(kwargs, gbk_list, checkm_results, db_file):
    db = DB.load_db(db_file, kwargs)
    tab = pd.read_table(checkm_results)

    hsh_taxid_to_gen_stats = get_gen_stats(gbk_list)

    hsh_filename_to_taxid = db.get_filenames_to_taxon_id()
    data = []
    for index, row in tab.iterrows():
        taxon_id = hsh_filename_to_taxid[row["Bin Id"]]
        completeness = row.Completeness
        contamination = row.Contamination
        gc, coding_density, length = hsh_taxid_to_gen_stats[row["Bin Id"]]
        values = [taxon_id, completeness,
                  contamination, gc, length, coding_density]
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
            accession = int(
                line[accession_offset:accession_offset + accession_length])
        elif line.startswith("#=GF DE"):
            description = line[len("#=GF DE   "):-1]


def load_pfam(params, pfam_files, db, pfam_def_file, ref_db_dir):
    db = DB.load_db(db, params)

    db.create_pfam_hits_table()
    pfam_ids = set()
    for pfam in pfam_files:
        entries = []
        with open(pfam, "r") as pfam_file:
            for line in pfam_file:
                if len(line) < len("CRC-") or line[0:len("CRC-")] != "CRC-":
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
        if entry.accession not in pfam_ids:
            continue
        pfam_entries.append([entry.accession, entry.description])
    pfam_def_file_iter.close()
    db.create_pfam_def_table(pfam_entries)
    db.set_status_in_config_table("pfam", 1)
    db.commit()

    # Determine reference DB version
    with open(os.path.join(ref_db_dir, "Pfam.version")) as fh:
        title, version = fh.readline().split(":")
        if "Pfam release" not in title:
            raise ValueError("Could not determine Pfam DB version")
    db.load_data_into_table("versions", [("Pfam", version.strip())])
    db.commit()


class ProtIdCounter(defaultdict):
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


vf_gene_id_expr = re.compile(r"(.*)\(gb\|(.*)\)")


def parse_vf_gene_id(to_parse):
    """IDs in vfdb.fasta are either of the form VFID(gb_accession) or VFID
    """
    if "(" in to_parse:
        return vf_gene_id_expr.match(to_parse).groups()
    return to_parse, None


vfdb_descr_expr = re.compile(r"\(.*?\) (.*?) \[.*?\((VF.*?)\) - (.*?) \((VFC.*?)\)\] \[(.*?)\]")


def parse_vfdb_entry(description):
    description = description.split(" ", 1)[1]
    prot_name, vfid, category, cat_id, organism = vfdb_descr_expr.match(
        description).groups()
    return prot_name, vfid, category, cat_id, organism


def load_swissprot(params, blast_results, db_name, swissprot_fasta, swissprot_db_dir):
    db = DB.load_db(db_name, params)
    hsh_swissprot_id = ProtIdCounter()
    db.create_swissprot_tables()

    # Note: this is not really scalable to x genomes, as
    # it necessitates to keep the prot id in memory
    # may need to process the blast results in batch instead (slower but
    # spares memory).
    for blast_file in blast_results:
        data = []
        with open(blast_file, "r") as blast_fh:
            for line in blast_fh:
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
        descr, taxid, org, gene, pe, version = parse_swissprot_entry(
            record.description)
        swiss_prot_defs.append(
            (db_prot_id, prot_id, descr, taxid, org, gene, pe))
    db.load_swissprot_defs(swiss_prot_defs)
    db.set_status_in_config_table("BLAST_swissprot", 1)
    db.commit()

    # Determine reference DB version
    with open(os.path.join(swissprot_db_dir, "relnotes.txt")) as fh:
        dbname, _,  version = fh.readline().split()
        if dbname != "UniProt":
            raise ValueError("Could not determine SwissProt DB version")
    db.load_data_into_table("versions", [("SwissProt", version.strip())])
    db.commit()


def load_vfdb_hits(params, blast_results, db_name, vfdb_fasta, vfdb_defs, min_seqid):
    db = DB.load_db(db_name, params)
    included_vf_genes = set()
    db.create_vf_tables()

    for blast_file in blast_results:
        data = []
        hits = pd.read_csv(
            blast_file, sep="\t", header=None,
            names=["crc", "gene_id", "seqid", "leng", "evalue", "score", "qcov"])

        # filter out hits with seqid too low
        hits = hits[hits["seqid"] >= min_seqid]
        # Select only the best hits: using pandas clearly is an overkill here
        hits = hits.sort_values(
            ["evalue", "seqid", "qcov"],
            ascending=[True, False, False]).drop_duplicates("crc")

        for index, row in hits.iterrows():
            hsh = simplify_hash(row.crc)
            vf_gene_id, _ = parse_vf_gene_id(row.gene_id)
            included_vf_genes.add(vf_gene_id)
            data.append((hsh, vf_gene_id, row.evalue, int(row.score),
                         row.seqid, row.leng, row.qcov))

        db.load_vf_hits(data)

    # Definitions are constructed from two data sources.
    # vfdb_fasta contains matching entries for every hit, with a
    # VFID relating to more information in VFs.xls
    vf_defs = pd.read_excel(vfdb_defs, header=1)
    vf_defs = vf_defs.set_index("VFID")

    vfdb_prot_defs = []
    for record in SeqIO.parse(vfdb_fasta, "fasta"):
        vf_gene_id, gb_accession = parse_vf_gene_id(record.name)
        if vf_gene_id not in included_vf_genes:
            continue
        prot_name, vfid, category, cat_id, organism = parse_vfdb_entry(
            record.description)
        # Get info from definitions.
        # Note that not all vfids have an entry in the definitions table
        if vfid in vf_defs.index:
            vf_data = vf_defs.loc[vfid]
        else:
            vf_data = {}
        vfdb_prot_defs.append(
            (vf_gene_id, gb_accession, prot_name, vfid, category, cat_id,
             vf_data.get("Characteristics"), vf_data.get("Structure"),
             vf_data.get("Function"), vf_data.get("Mechanism")))
    db.load_vf_defs(vfdb_prot_defs)
    db.set_status_in_config_table("BLAST_vfdb", 1)
    db.commit()

    # There is no version of the VFdb, only a download date
    # As there is a new version of the DB every week, the download date
    # will have to do.
    book = xlrd.open_workbook(vfdb_defs)
    sheet = book.sheet_by_index(0)
    download_date = re.match(r".*\[(.*)\]", sheet.row(0)[3].value).groups()[0]
    download_date = datetime.strptime(download_date, "%a %b %d %H:%M:%S %Y")
    db.load_data_into_table("versions", [("VFDB", download_date.date().isoformat())])
    db.commit()


# NOTE:
# Several KO marked as significant can be assigned to the same locus
# only take the hit with the lowest evalue (the first in the list)
def load_KO(params, ko_files, db_name, ko_db_dir):
    db = DB.load_db(db_name, params)
    data = []
    for ko_file in ko_files:
        curr_hsh = None
        with open(ko_file, "r") as ko_fh:
            for ko_line in ko_fh:
                tokens = ko_line.split()
                # ignore all but the best hits
                if tokens[0] != "*":
                    continue
                crc_raw, ko_str, thrs_str, score_str, evalue_str, * \
                    descr = tokens[1:]
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

    # Determine reference DB version
    with open(os.path.join(ko_db_dir, "version.txt")) as fh:
        version = fh.readline()
    db.load_data_into_table("versions", [("Ko", version.strip())])
    db.commit()


def load_module_completeness(params, db_name):
    db = DB.load_db(db_name, params)
    all_genomes = db.get_genomes_description().description.to_dict()

    db.create_module_completeness_table()

    for taxid, _ in all_genomes.items():
        complete_modules = []
        ko_count = db.get_ko_hits([taxid])

        ko_list = ko_count.index.unique().tolist()
        if len(ko_list) == 0:
            continue

        modules = db.get_ko_modules(ko_list, as_pandas=True)
        module_def = db.get_modules_info(modules.module_id.unique().tolist())

        for module_id, descr, definit, _, _ in module_def:
            parser = KO_module.ModuleParser(definit)
            ko_set = {ko: 1 for ko in ko_list}
            expr_tree = parser.parse()
            if expr_tree.get_n_missing(ko_set) == 0:
                complete_modules.append((module_id, taxid))
        db.load_module_completeness(complete_modules)
    db.commit()


def format_og(og_n):
    return f"group_{int(og_n)}"


def setup_chlamdb_search_index(params, db_name, index_name):
    db = DB.load_db(db_name, params)
    os.mkdir(index_name)

    has_cog = params.get("cog", False)
    has_ko = params.get("ko", False)
    has_pfam = params.get("pfam_scan", False)
    has_amr = params.get("amr", False)
    has_vf = params.get("vfdb", False)

    genomes = db.get_genomes_description()
    index = search_bar.ChlamdbIndex.new_index(index_name)

    for taxid, data in genomes.iterrows():
        all_infos = db.get_proteins_info([taxid], search_on="taxid",
                                         as_df=True, inc_non_CDS=True, inc_pseudo=True)
        ogs = db.get_og_count(all_infos.index.tolist(), search_on="seqid")
        all_infos = all_infos.join(ogs)
        for seqid, data in all_infos.iterrows():
            locus_tag = data.locus_tag
            gene, product = None, None
            if "product" in data:
                product = data["product"]
                if pd.isna(product):
                    product = None
                elif product == "hypothetical protein":
                    product = None

            gene = None
            if "gene" in data:
                gene = data["gene"]
                if pd.isna(gene):
                    gene = None

            og = None
            if not pd.isna(data.orthogroup):
                og = format_og(data.orthogroup)

            organism = genomes.loc[taxid].description
            search_bar.GeneEntry().add_to_index(
                index, locus_tag, gene, product, organism, og)

    if has_cog:
        cog_data = db.get_cog_summaries(cog_ids=None, only_cog_desc=True)
        for cog, (func, descr) in cog_data.items():
            search_bar.CogEntry().add_to_index(index, cog, descr)

    if has_ko:
        ko_data = db.get_ko_desc(ko_ids=None)
        for ko, descr in ko_data.items():
            search_bar.KoEntry().add_to_index(index, ko, descr)

        mod_data = db.get_modules_info(ids=None, search_on=None)
        for mod_id, mod_desc, _, _, _ in mod_data:
            search_bar.ModuleEntry().add_to_index(index, mod_id, mod_desc)

        pat_data = db.get_pathways()
        for pat_id, path_desc in pat_data:
            search_bar.PathwayEntry().add_to_index(index, pat_id, path_desc)

    if has_pfam:
        pfam_data = db.get_pfam_def(pfam_ids=None)
        for pfam, data in pfam_data.iterrows():
            search_bar.PfamEntry().add_to_index(index, pfam, data["def"])

    if has_amr:
        amr_data = db.get_amr_descriptions()
        for amr, data in amr_data.iterrows():
            search_bar.AmrEntry().add_to_index(index, data["gene"], data["seq_name"])

    if has_vf:
        vf_data = db.vf.get_hit_descriptions(hit_ids=None)
        for vf, data in vf_data.iterrows():
            search_bar.VfEntry().add_to_index(index, data["vf_gene_id"], data["prot_name"])

    index.done_adding()
