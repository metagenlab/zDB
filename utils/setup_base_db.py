#!/usr/bin/env python3
"""
This script is used to generate the skeleton of the database
used by zdb. The entries will be filled by the annotation
pipeline.

For now, this script creates a new mysql database and creates
the tables of the biosql schema.
It then pulls all kegg orthologs, their related modules and pathways
and enters them in the database.
Same with COGs.

##############################################################################
NOTE: for this script to work, you either need to disable multithreading by
setting N_THREAD to 1 or you need to patch REST.py to use urllib3 (thread safe).


TODO: as the REST API does not seem to be under active development, it would
be a solution to just import the source files into this directory to solve this
problem.
##############################################################################

Bastian Marquis (bastian.marquis@protonmail.com)
Date: 29.10.2020
"""

import argparse
import os
import sys
from collections import defaultdict
from urllib.request import urljoin
from urllib.request import urlretrieve

# to be removed in favor of a local version
from Bio.KEGG import REST

zdbdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(zdbdir, "webapp"))
from lib import db_utils  # noqa

# from REST documentation, can get a max of 10 queries
# in kegg_get
MAX_N_QUERIES = 10

# number of threads that will be used simultaneously
# to download the ko genes. 5 seems a good compromise
# between being blacklisted and speed
DEFAULT_N_THREADS = 5
DEFAULT_KO_DIR = "tmp"

N_RETRY = 5


def create_data_table(db):
    entry_list = [
        ("orthology", "mandatory", False),
        ("orthogroup_alignments", "mandatory", False),
        ("reference_phylogeny", "mandatory", False),
        ("taxonomy_table", "mandatory", False),
        ("genome_statistics", "mandatory", False),
        ("gene_phylogenies", "optional", False),
        ("COG", "optional", False),
        ("KEGG", "optional", False),
        ("pfam", "optional", False),
        ("BLAST_refseq", "optional", False),
        ("BLAST_swissprot", "optional", False),
        ("BBH_phylogenies", "optional", False),
        ("GC_statistics", "optional", False),
        ("phylogenetic_profile", "optional", False),
        ("gbk_files", "optional", False),
        ("AMR", "optional", False),
        ("BLAST_vfdb", "optional", False),
        ("GIS", "optional", False),
    ]
    db.load_chlamdb_config_tables(entry_list)
    db.commit()


def create_versions_table(db):
    db.create_versions_table()
    db.commit()


def setup_biodb(kwargs):
    sqlpsw = kwargs["db_psswd"]
    db_type = kwargs["db_type"]
    db_name = kwargs["db_name"]
    sql_def_file = kwargs.get("biosql_schema", "")
    if not sql_def_file:
        schema_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "biosql_schema"
        )
        sql_def_file = [os.path.join(schema_dir, f"biosqldb-{db_type}.sql")]

    if db_type == "sqlite":
        import sqlite3

        conn = sqlite3.connect(db_name)
        cursor = conn.cursor()
        err_code = os.system(f"sqlite3 {db_name} < {sql_def_file[0]}")
        conn.execute("pragma journal_mode=wal")
    else:
        # NOTE: this part is untested!
        import MySQLdb

        conn = MySQLdb.connect(
            host="localhost",  # your host, usually localhost
            user="root",  # your username
            passwd=sqlpsw,  # name of the data base
        )
        cursor = conn.cursor()
        sql_db = f"CREATE DATABASE IF NOT EXISTS {db_name};"
        cursor.execute(
            sql_db,
        )
        conn.commit()
        cursor.execute(
            f"use {db_name};",
        )
        url_biosql_scheme = "biosqldb-mysql.sql"
        schema_dir = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "biosql_schema"
        )
        err_code = os.system(
            f"mysql -uroot -p{sqlpsw} {db_name} < {schema_dir}/{url_biosql_scheme}"
        )

    if err_code != 0:
        raise IOError("Problem loading sql schema:", err_code)

    # not really logical to me, but creating a database
    # from the biosql is necessary
    chlamdb_args = {"zdb.db_type": db_type, "zdb.db_name": db_name, "zdb.psswd": sqlpsw}
    db = db_utils.DB.load_db(db_name, chlamdb_args)
    db.create_biosql_database(chlamdb_args)
    db.commit()
    return db


class Record(object):
    def __init__(self, entry, name, modules, pathways):
        self.entry = entry
        self.name = name
        self.modules = modules
        self.pathways = pathways

    def simplified_entry(self):
        return int(self.entry[len("K") :])

    def __str__(self):
        acc = self.entry + " " + self.definition
        return acc


class Module(object):
    def __init__(self):
        self.entry = None
        self.descr = None
        self.classes = []
        self.definition = []

    def simplified_entry(self):
        return int(self.entry[len("M") :])

    def get_definition(self):
        if len(self.definition) == 1:
            return self.definition[0]
        # Module definitions may come on multiple lines, join them
        # with an AND
        return " ".join(f"{part}" for part in self.definition)

    def sub_category(self):
        return self.classes[-1]

    def category(self):
        return self.classes[-2]

    def is_signature(self):
        return self.classes[0] == "Signature modules"

    def __str__(self):
        return self.entry + ": " + self.descr

    def __hash__(self):
        return hash(self.entry)

    def __eq__(self, other):
        return self.entry == other.entry


class Pathway(object):
    def __init__(self, entry, descr):
        self.entry = entry
        self.descr = descr

    def simplified_entry(self):
        return int(self.entry[len("map") :])

    def __hash__(self):
        return hash(self.entry)

    def __eq__(self, other):
        return self.entry == other.entry

    def __str__(self):
        return self.entry + ": " + self.descr


def parse_module(handle):
    module = Module()
    for line in handle:
        if line[:3] == "///":
            yield module
            module = Module()
            continue
        if line[:12] != "            ":
            keyword = line[:12].strip()
        data = line[12:].strip()
        if keyword == "ENTRY":
            module.entry = data.split()[0]
        if keyword == "NAME":
            module.descr = data
        if keyword == "CLASS":
            tokens = data.split(";")
            module.classes = [token.strip() for token in tokens]
        if keyword == "DEFINITION":
            module.definition.append(data)


def chunk_modules(modules, size=MAX_N_QUERIES):
    curr = []
    for i, module in enumerate(modules):
        curr.append(module)
        if i % MAX_N_QUERIES == 0:
            yield curr
            curr = []
    if len(curr) > 0:
        yield curr


def get_ko_modules_mapping():
    data = REST.kegg_link("module", "ko")
    mapping = defaultdict(list)
    for line in data:
        entry, module = line.split("\t")
        entry = entry.split("ko:")[1].strip()
        module = module.split("md:")[1].strip()
        # it seems there are a few dupplicate
        if module not in mapping[entry]:
            mapping[entry].append(module)
    return mapping


def get_ko_pathways_mapping():
    pathway_mapping = {}
    for line in REST.kegg_list("pathway"):
        entry, descr = line.split("\t")
        entry = entry.strip()
        descr = descr.strip()
        pathway_mapping[entry] = Pathway(entry, descr)

    data = REST.kegg_link("pathway", "ko")
    mapping = defaultdict(list)
    for line in data:
        entry, pathway_id = line.split("\t")
        entry = entry.split("ko:")[1].strip()
        pathway_id = pathway_id.split("path:")[1].strip()
        if pathway_id.startswith("map"):
            pathway = pathway_mapping[pathway_id]
            # Some entries are dupplicate
            if pathway not in mapping[entry]:
                mapping[entry].append(pathway)
    return mapping


def load_KO_references(db, params, ko_dir=DEFAULT_KO_DIR):
    ko_modules = get_ko_modules_mapping()
    ko_pathways = get_ko_pathways_mapping()

    genes = []
    for line in REST.kegg_list("KO"):
        entry, rest = line.split("\t", 1)
        if ";" in rest:
            name = rest.split(";", 1)[1]  # Record for K17548 has two semicola.
        else:
            # record K23479 does not have a semicolon
            name = rest.split("/")[1]
        entry = entry.strip()
        name = name.strip()
        genes.append(Record(entry, name, ko_modules[entry], ko_pathways[entry]))

    # get only the module present in the genes
    modules_set = {module for gene in genes for module in gene.modules}
    hsh_modules = {}
    print("Downloading module data")
    for i, module_list in enumerate(chunk_modules(modules_set)):
        if i % 10 == 0:
            print(i * MAX_N_QUERIES)
        mod_data = REST.kegg_get(module_list)
        for module in parse_module(mod_data):
            hsh_modules[module.entry] = module
    print("Done")

    pathway_set = set()
    category_set = set()
    for entry, module in hsh_modules.items():
        if len(module.classes) > 0:
            category_set.add(module.sub_category())
            category_set.add(module.category())
    for gene in genes:
        for pathway in gene.pathways:
            pathway_set.add(pathway)

    # really ugly, but spares some web requests
    module_classes = []
    hsh_class_to_id = {}
    for cat_id, cat in enumerate(category_set):
        module_classes.append((cat_id, cat))
        hsh_class_to_id[cat] = cat_id
    for entry, module in hsh_modules.items():
        cat_id = hsh_class_to_id[module.category()]
        subcat_id = hsh_class_to_id[module.sub_category()]
        module.cat_id = cat_id
        module.subcat_id = subcat_id

    db.load_ko_module_classes(module_classes)
    db.load_ko_module(
        [
            (
                m.simplified_entry(),
                m.descr,
                m.get_definition(),
                m.is_signature(),
                m.cat_id,
                m.subcat_id,
            )
            for m in hsh_modules.values()
        ]
    )
    db.load_ko_pathway([(p.simplified_entry(), p.descr) for p in pathway_set])
    db.load_ko_def([(gene.simplified_entry(), gene.name) for gene in genes])

    ko_to_path = []
    ko_to_module = []
    for gene in genes:
        simp = gene.simplified_entry()
        ko_to_path.extend((simp, path.simplified_entry()) for path in gene.pathways)
        ko_to_module.extend(
            (simp, hsh_modules[mod].simplified_entry()) for mod in gene.modules
        )
    db.load_ko_to_pathway(ko_to_path)
    db.load_ko_to_module(ko_to_module)
    db.commit()
    return genes


def download_cog_files(cog_dir):
    base_url = "ftp://ftp.ncbi.nih.gov/pub/COG/COG2020/data/"
    files = ["cog-20.def.tab", "fun-20.tab"]
    for file in files:
        urlretrieve(urljoin(base_url, file), os.path.join(cog_dir, file))


def setup_cog(db, cog_dir):
    fun_names_file = open(cog_dir + "/fun-20.tab", "r")
    fun_names_file.reconfigure(encoding="Latin-1")
    cog_names_file = open(cog_dir + "/cog-20.def.tab", "r")
    cog_names_file.reconfigure(encoding="Latin-1")

    cog_ref_data = []
    for line in cog_names_file:
        tokens = line.strip().split("\t")
        cog_id = int(tokens[0][3:])
        fun = tokens[1].strip()
        description = tokens[2].strip()
        cog_ref_data.append((cog_id, fun, description))
    db.load_cog_ref_data(cog_ref_data)

    cog_fun_data = []
    for line in fun_names_file:
        function, random_str, description = line.split("\t")
        cog_fun_data.append((function.strip(), description.strip()))
    db.load_cog_fun_data(cog_fun_data)
    db.commit()
    fun_names_file.close()
    cog_names_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Creates a chlamdb database skeleton")

    parser.add_argument(
        "--db_name",
        nargs="?",
        default="zdb_base",
        help="name of the database (default name zdb_base)",
    )

    parser.add_argument(
        "--load_cog", action="store_true", help="load cog definitions (default no)"
    )

    parser.add_argument(
        "--cog_dir",
        nargs="?",
        default="./",
        help="directory where the cog definitions files are "
        "(default current directory)",
    )

    parser.add_argument(
        "--load_kegg",
        action="store_true",
        help="load kegg definitions (default no, must specify ko genes dir)",
    )

    parser.add_argument(
        "--db_type",
        nargs="?",
        default="sqlite",
        help="database type (either sqlite or mysql)",
    )

    parser.add_argument(
        "--db_psswd", nargs="+", default="", help="set db password (default none)"
    )

    parser.add_argument(
        "--ko_dir",
        nargs="?",
        default=DEFAULT_KO_DIR,
        help=f"ko directory, defaults to {DEFAULT_KO_DIR}",
    )

    parser.add_argument(
        "--skip_biodb",
        action="store_true",
        default=False,
        help="skip setting up biodb, might be useful if using a pre-existing db",
    )

    parser.add_argument(
        "--biosql_schema", nargs="+", help="location of the biosql schema"
    )

    args = vars(parser.parse_args())

    db = None
    if not args.get("skip_biodb"):
        print("Setting up the biosql schema")
        db = setup_biodb(args)
        create_data_table(db)
        create_versions_table(db)

    if db is None:
        db_name = args["db_name"]
        db_type = args["db_type"]
        db_psswd = args["db_psswd"]
        chlamdb_args = {
            "zdb.db_type": db_type,
            "zdb.db_name": db_name,
            "zdb.db_psswd": db_psswd,
        }
        db = db_utils.DB.load_db(db_name, chlamdb_args)

    if args.get("load_cog", False):
        print("Downloading COG tables")
        cog_dir = args.get("cog_dir", ".")
        download_cog_files(cog_dir)
        print("Loading COG tables")
        setup_cog(db, cog_dir)

    if args.get("load_kegg", False):
        print("Loading KEGG tables")
        ko_dir = args.get("ko_dir", DEFAULT_KO_DIR)
        genes = load_KO_references(db, args, ko_dir)
