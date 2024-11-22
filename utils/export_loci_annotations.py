import argparse
import os
import sys

import pandas as pd

zdbdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
print(zdbdir)
sys.path.append(os.path.join(zdbdir, "webapp"))
from lib import db_utils  # noqa


def format_cog(cog_id):
    return f"COG{int(cog_id):04d}"


def format_ko(ko_id):
    return f"K{int(ko_id):05d}"


def format_module(module_id):
    return f"M{module_id:05d}"


def format_pathway(pathway_id):
    return f"map{pathway_id:05d}"


def get_tables(db_file, taxonid):
    zdb_params = {'zdb.db_type': 'sqlite',
                  'zdb.db_name': 'George',
                  'zdb.psswd': ''}
    db = db_utils.DB.load_db(db_file, zdb_params)

    locus_term_id = db.server.adaptor.execute_one(
        "SELECT term_id FROM term WHERE name='locus_tag'")[0]

    query = (f"SELECT seqfeature_qualifier_value.seqfeature_id, seqfeature_qualifier_value.value FROM seqfeature_qualifier_value "
             f"INNER JOIN seqfeature ON seqfeature_qualifier_value.seqfeature_id=seqfeature.seqfeature_id "
             f"INNER JOIN bioentry ON bioentry.bioentry_id=seqfeature.bioentry_id "
             f"WHERE seqfeature_qualifier_value.term_id = {locus_term_id} AND bioentry.taxon_id = {taxonid}")
    loci = pd.DataFrame(db.server.adaptor.execute_and_fetchall(query),
                        columns=["seqid", "locus"]).set_index("seqid")

    query = ("SELECT sequence_hash_dictionnary.seqid, cog_hits.cog_id "
             "FROM cog_hits INNER JOIN sequence_hash_dictionnary "
             "ON sequence_hash_dictionnary.hsh=cog_hits.hsh")

    cog_hits = pd.DataFrame(db.server.adaptor.execute_and_fetchall(query),
                            columns=["seqid", "cog"]).set_index("seqid")

    cogs = loci.merge(cog_hits, left_index=True, right_index=True)

    query = ("SELECT sequence_hash_dictionnary.seqid, ko_hits.ko_id "
             "FROM ko_hits INNER JOIN sequence_hash_dictionnary "
             "ON sequence_hash_dictionnary.hsh=ko_hits.hsh")

    ko_hits = pd.DataFrame(db.server.adaptor.execute_and_fetchall(query),
                           columns=["seqid", "ko"]).set_index("seqid")

    kos = loci.merge(ko_hits, left_index=True, right_index=True)

    query = "SELECT ko_id, module_id from ko_to_module"
    ko_to_modules = pd.DataFrame(db.server.adaptor.execute_and_fetchall(query),
                                 columns=["ko", "module"])
    modules = kos.merge(ko_to_modules, left_on="ko", right_on="ko")
    modules = modules.drop("ko", axis=1)
    modules = modules.drop_duplicates()

    query = "SELECT ko_id, pathway_id from ko_to_pathway"
    ko_to_pathways = pd.DataFrame(db.server.adaptor.execute_and_fetchall(query),
                                  columns=["ko", "pathway"])
    pathways = kos.merge(ko_to_pathways, left_on="ko", right_on="ko")
    pathways = pathways.drop("ko", axis=1)
    pathways = pathways.drop_duplicates()

    cogs["cog"] = cogs["cog"].map(format_cog)
    kos["ko"] = kos["ko"].map(format_ko)
    modules["module"] = modules["module"].map(format_module)
    pathways["pathway"] = pathways["pathway"].map(format_pathway)
    return cogs, kos, modules, pathways


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Exports COG and KO table (ids vs loci)")

    parser.add_argument("db", help="Path to DB file")
    parser.add_argument("taxid", help="taxid")

    parser.add_argument("--outdir", nargs="?", default=".",
                        help="directory to which the tables should be written to.")

    args = parser.parse_args()
    cogs, kos, modules, pathways = get_tables(args.db, int(args.taxid))
    cogs.to_csv(os.path.join(args.outdir, "cogs.tsv"), index=False, sep="\t")
    kos.to_csv(os.path.join(args.outdir, "kos.tsv"), index=False, sep="\t")
    modules.to_csv(os.path.join(args.outdir, "modules.tsv"), index=False, sep="\t")
    pathways.to_csv(os.path.join(args.outdir, "pathways.tsv"), index=False, sep="\t")
