import json

import numpy as np
from chlamdb.forms import make_circos_form
from django.conf import settings
from django.shortcuts import render
from lib.db_utils import DB
from matplotlib.cm import tab10
from matplotlib.cm import tab20
from views.object_type_metadata import my_locals
from views.utils import page2title


class CircosData:
    """
    Prepare json for CGView
    """

    def __init__(self, figure_div_id="heatmapChart"):
        self.features = []
        self.contigs = []
        self.plots = []
        self.settings = {}
        self.legend_items = []
        self.tracks = []

    def to_json(self):
        return {
            "cgview": {
                "version": "1.7.0",
                "features": self.features,
                "sequence": {"contigs": self.contigs},
                "plots": self.plots,
                "legend": {
                    "items": self.legend_items,
                    "position": "top-left",
                    "backgroundColor": "rgba(255,255,255,0.75)",
                },
                "tracks": self.tracks,
                "dividers": {
                    "track": {
                        "color": "rgba(0,0,0,0.7)",
                        "thickness": 3,
                        "spacing": 1.5,
                    },
                    "slot": {
                        "color": "rgba(128,128,128,0.5)",
                        "thickness": 1,
                        "spacing": 1,
                    },
                },
            }
        }

    def add_contigs_data(self, df_bioentry):
        """
        Input df with the following columns:
        - length
        - accession

        The index is used as id
        """
        for i, (bioentry, row) in enumerate(df_bioentry.iterrows()):
            self.contigs.append(
                {
                    "length": row.length,
                    "color": "#8dd3c7" if i % 2 == 0 else "#fb8072",
                    "name": bioentry,
                    "meta": {"name": row.accession},
                    "seq": row.seq,
                }
            )

    def add_gene_track(self, df):
        """
        Input df columns:
        - bioentry_id
        - start_pos
        - end_pos
        - locus_ref
        """
        loci = [
            {
                "name": row.locus_ref,
                "source": "reference",
                "contig": str(row.bioentry_id),
                "start": row.start_pos,
                "stop": row.end_pos,
                "strand": row.strand,
                "type": row.term_name,
                "meta": {"gene": f"{row.gene}", "product": f"{row.gene_product}"},
                "legend": row.legend,
            }
            for n, row in df.iterrows()
        ]
        self.features.extend(loci)
        self.tracks.append(
            {
                "name": "reference",
                "separateFeaturesBy": "strand",
                "position": "both",
                "thicknessRatio": 1,
                "dataType": "feature",
                "dataMethod": "source",
                "dataKeys": "reference",
            }
        )

    def add_heatmap_data(self, df, label, color):
        """
        Input df columns:
        - bioentry_id
        - start_pos
        - end_pos
        - identity
        - locus_tag

        Locus without homologs in the target genome should have a locus_tag with a value of None
        """
        heatmap_data = [
            {
                "contig": row.bioentry_id,
                "start": row.start_pos,
                "stop": row.end_pos,
                "type": row.term_name,
                "meta": {
                    "gene": f"{row.gene}",
                    "product": f"{row.gene_product}",
                    "identity": f"{row.identity}%",
                },
                "score": row.identity / 100,
                "name": row.locus_tag,
                "source": "targets",
                "legend": label,
            }
            for n, row in df.iterrows()
            if row.locus_tag is not None
        ]
        self.features.extend(heatmap_data)
        self.legend_items.append(
            {
                "name": label,
                "decoration": "score",
                "swatchColor": color,
            }
        )

    def add_heatmap_track(self):
        self.tracks.append(
            {
                "name": "targets",
                "position": "inside",
                "thicknessRatio": 1,
                "dataType": "feature",
                "dataMethod": "source",
                "dataKeys": "targets",
                "separateFeaturesBy": "legend",
            }
        )

    def add_histogram_track(self, df, label, color):
        """
        Minimal Input df columns:
        - bioentry_id
        - start_pos
        - end_pos
        - value
        """
        histogram_data = [
            {
                "contig": str(int(row.bioentry_id)),
                "start": row.start_pos,
                "stop": row.end_pos,
                "score": row.value,
                "source": label,
                "legend": label,
                "name": f"group_{n}",
                "type": "OG",
                "meta": {"prevalence": f"{int(100 * row.value)}%"},
            }
            for n, row in df.iterrows()
        ]
        self.features.extend(histogram_data)
        self.tracks.append(
            {
                "name": label,
                "position": "outside",
                "thicknessRatio": 1,
                "dataType": "feature",
                "dataMethod": "source",
                "dataKeys": label,
                "separateFeaturesBy": "none",
            }
        )
        self.legend_items.append(
            {
                "name": label,
                "decoration": "score",
                "swatchColor": color,
            }
        )


def get_circos_data(reference_taxon, target_taxons, highlight_og=None):
    biodb_path = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb_path)

    # "bioentry_id", "accession" ,"length"
    df_bioentry = db.get_bioentry_list(reference_taxon, min_bioentry_length=1000)

    # "bioentry_id", "seqfeature_id", "start_pos", "end_pos", "strand"
    df_feature_location = db.get_features_location(
        reference_taxon, search_on="taxon_id", seq_term_names=["CDS", "rRNA", "tRNA"]
    ).set_index(["seqfeature_id"])
    # df of target genomes
    df_targets = db.get_proteins_info(
        ids=target_taxons,
        search_on="taxid",
        as_df=True,
        to_return=["locus_tag"],
        inc_non_CDS=False,
        inc_pseudo=False,
    )

    # retrieve n_orthologs of list of seqids
    seq_og = db.get_og_count(df_feature_location.index.to_list(), search_on="seqid")
    count_all_genomes = db.get_og_count(
        seq_og["orthogroup"].to_list(), search_on="orthogroup"
    )
    n_genomes = db.get_number_of_genomes()
    orthogroup2frac_all = (
        count_all_genomes[count_all_genomes > 0].count(axis=1) / n_genomes
    )
    homologs_count = (
        df_feature_location.loc[df_feature_location.term_name == "CDS"]
        .join(seq_og)
        .reset_index()
        .set_index("orthogroup")
        .merge(orthogroup2frac_all.rename("value"), left_index=True, right_index=True)[
            ["bioentry_id", "start_pos", "end_pos", "value"]
        ]
    )

    # this query can be pretty slow
    df_identity = db.get_identity_closest_homolog(
        reference_taxon, target_taxons
    ).set_index(["target_taxid"])

    c = CircosData()
    c.add_contigs_data(df_bioentry)

    # sort taxons by number of homologs (from mot similar to most dissmilar)
    target_taxon_n_homologs = (
        df_identity.groupby(["target_taxid"])
        .count()["seqfeature_id_1"]
        .sort_values(ascending=False)
    )
    # "bioentry_id", "seqfeature_id", "start_pos", "end_pos", "strand"
    # "seqfeature_id_1", "seqfeature_id_2", "identity", "target_taxid"
    # join on seqfeature id
    df_feature_location["gene"] = df_feature_location["qualifier_value_gene"].fillna(
        "-"
    )
    df_feature_location["gene_product"] = df_feature_location[
        "qualifier_value_product"
    ].fillna("-")
    df_feature_location["locus_tag"] = df_feature_location[
        "qualifier_value_locus_tag"
    ].fillna("-")

    if highlight_og is not None:
        df_feature_location["legend"] = "locus"
        df_genes = db.get_genes_from_og(
            highlight_og, taxon_ids=[reference_taxon], terms=["locus_tag"]
        )
        df_feature_location.loc[df_genes["locus_tag"].index, "legend"] = (
            "extracted locus"
        )
    else:
        df_feature_location["legend"] = df_feature_location["term_name"]

    df_feature_location = df_feature_location.rename(columns={"locus_tag": "locus_ref"})
    c.add_gene_track(df_feature_location)

    import matplotlib as mpl

    c.add_histogram_track(
        homologs_count,
        "Prevalence",
        mpl.colors.get_named_colors_mapping()["aquamarine"],
    )

    # iterate ordered list of target taxids, add track to circos
    if len(target_taxon_n_homologs.index) <= 10:
        taxon_colors = {
            el: tab10(i) for i, el in enumerate(target_taxon_n_homologs.index)
        }
    else:
        taxon_colors = {
            el: tab20(i % 20) for i, el in enumerate(target_taxon_n_homologs.index)
        }
    for target_taxon in target_taxon_n_homologs.index:
        df_combined = df_feature_location.join(
            df_identity.loc[target_taxon]
            .reset_index()
            .set_index("seqfeature_id_1")
            .rename_axis("seqfeature_id")
        ).reset_index()
        df_combined.identity = df_combined.identity.fillna(0).astype(int)
        df_combined.bioentry_id = df_combined.bioentry_id.astype(str)

        # only keep the highest identity for each seqfeature id
        df_combined = (
            df_combined.sort_values("identity", ascending=False)
            .drop_duplicates(subset=["seqfeature_id", "start_pos"])
            .sort_index()
        )

        df_combined = df_combined.join(
            df_targets, on="seqfeature_id_2", how="left"
        ).reset_index()
        df_combined["locus_tag"] = (
            df_combined["locus_tag"].fillna(np.nan).replace([np.nan], [None])
        )
        c.add_heatmap_data(
            df_combined,
            f"target_{target_taxon}",
            mpl.colors.rgb2hex(taxon_colors[target_taxon]),
        )
    c.add_heatmap_track()
    return c.to_json()


def circos(request):
    biodb_path = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb_path)
    page_title = page2title["circos"]

    circos_form_class = make_circos_form(db)

    if request.method == "POST":
        form = circos_form_class(request.POST)

        if form.is_valid():
            target_taxons = form.get_target_taxids()
            reference_taxon = form.get_ref_taxid()

            form_display = "inherit"
            highlighted_ogs = None
            if "highlighted_ogs" in form.data:
                # This is only set when coming from the OG extraction view.
                # As this is not a field supported in the form, we hide the form
                highlighted_ogs = form.data.getlist("highlighted_ogs")
                form_display = None

            circos_json = json.dumps(
                get_circos_data(reference_taxon, target_taxons, highlighted_ogs)
            )
            envoi = True
    else:
        form = circos_form_class()

    local_vars = my_locals(locals())
    return render(request, "chlamdb/circos.html", local_vars)
