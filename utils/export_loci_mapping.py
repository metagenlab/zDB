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


def get_table(db_file, ref_taxid, taxids):
    zdb_params = {"zdb.db_type": "sqlite", "zdb.db_name": "George", "zdb.psswd": ""}
    db = db_utils.DB.load_db(db_file, zdb_params)

    locus_term_id = db.server.adaptor.execute_one(
        "SELECT term_id FROM term WHERE name='locus_tag'"
    )[0]

    cds_term_id = db.server.adaptor.execute_one(
        "SELECT term_id FROM term WHERE name='CDS'"
    )[0]

    query = (
        f"SELECT seqfeature_qualifier_value.seqfeature_id, seqfeature_qualifier_value.value, og_hits.orthogroup FROM seqfeature_qualifier_value "
        f"INNER JOIN seqfeature ON seqfeature_qualifier_value.seqfeature_id=seqfeature.seqfeature_id "
        f"INNER JOIN bioentry ON bioentry.bioentry_id=seqfeature.bioentry_id "
        f"INNER JOIN og_hits ON og_hits.seqid=seqfeature.seqfeature_id "
        f"WHERE seqfeature.type_term_id = {cds_term_id} AND seqfeature_qualifier_value.term_id = {locus_term_id} AND bioentry.taxon_id = {ref_taxid}"
    )
    loci = pd.DataFrame(
        db.server.adaptor.execute_and_fetchall(query),
        columns=["seqid", "locus", "orthogroup"],
    ).set_index("seqid")

    og_annot = db.get_genes_from_og(
        orthogroups=[int(el) for el in loci["orthogroup"].unique()],
        terms=["locus_tag", "gene"],
    )

    all_org = db.get_organism(og_annot.index.tolist(), as_df=True)

    # We add the gene names to the loci table
    loci = loci.merge(og_annot["gene"], left_index=True, right_index=True, how="left")

    # Now we add empty columns for each organism
    organisms = db.get_genomes_description(taxids)["description"]
    for org in organisms:
        loci[org] = ""

    # Now for each locus we will get the homologs with identities and gene name
    for seqid, row in loci.iterrows():
        identities = db.get_og_identity(og=row.orthogroup, ref_seqid=seqid)
        identities["identity"] = identities["identity"].apply(lambda x: f"{x:.1f}")
        og_loci = identities.merge(og_annot, left_index=True, right_index=True).merge(
            all_org, left_index=True, right_index=True
        )
        if og_loci.empty:
            continue

        og_loci.insert(
            len(og_loci.columns),
            "label",
            og_loci[["locus_tag", "gene", "identity"]].map(str).agg("|".join, axis=1),
        )

        matching_loci = (
            og_loci[["organism", "label"]]
            .groupby("organism")
            .agg(lambda x: ",".join(x))
            .label
        )
        for org in organisms:
            loci.loc[seqid, org] = matching_loci.get(org, "")

    return loci


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Exports a mapping between the loci of a reference genome and the loci of other genomes."
    )

    parser.add_argument("db", help="Path to DB file")
    parser.add_argument("ref_taxid", help="reference taxid")
    parser.add_argument(
        "taxids", help="Taxids in which to find matching loci", nargs="*"
    )

    parser.add_argument(
        "--outdir",
        nargs="?",
        default=".",
        help="directory to which the table should be written to.",
    )

    args = parser.parse_args()
    table = get_table(args.db, int(args.ref_taxid), [int(t) for t in args.taxids])
    table.to_csv(os.path.join(args.outdir, "loci_map.tsv"), index=False, sep="\t")
