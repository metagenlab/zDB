import json
from operator import attrgetter

import matplotlib as mpl
import numpy as np
from chlamdb.forms import make_circos_form
from django.shortcuts import render
from django.views import View
from matplotlib.cm import tab10
from matplotlib.cm import tab20
from views.mixins import BaseViewMixin
from views.utils import optional2status


class CircosData:
    """
    Prepare json for CGView
    """

    def __init__(
        self, with_amr=False, with_vf=False, with_highlighted_ogs=False, with_gi=False
    ):
        self.with_amr = with_amr
        self.with_vf = with_vf
        self.with_highlighted_ogs = with_highlighted_ogs
        self.with_gi = with_gi
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
                        "color": "rgba(0,0,0,0.0)",
                        "thickness": 1,
                        "spacing": 8,
                    },
                    "slot": {
                        "color": "rgba(128,128,128,0.0)",
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

    def _gene_meta_extractor(self):
        meta_map = {
            "gene": attrgetter("gene"),
            "product": attrgetter("gene_product"),
        }
        if self.with_highlighted_ogs:
            meta_map["extracted"] = attrgetter("extracted")
        if self.with_amr:
            meta_map["AMR"] = attrgetter("amr")
        if self.with_vf:
            meta_map["VF"] = attrgetter("vf")

        return lambda row: {key: fun(row) for key, fun in meta_map.items()}

    def _legend_extractor(self):
        if self.with_highlighted_ogs:
            return lambda x: f"extracted {x.extracted}"
        return attrgetter("term_name")

    def add_gene_track(self, df):
        """
        Input df columns:
        - bioentry_id
        - start_pos
        - end_pos
        - locus_ref
        """

        meta_extractor = self._gene_meta_extractor()
        legend_extractor = self._legend_extractor()
        loci = [
            {
                "name": row.locus_ref,
                "source": "reference",
                "contig": str(row.bioentry_id),
                "start": row.start_pos,
                "stop": row.end_pos,
                "strand": row.strand,
                "type": row.term_name,
                "meta": meta_extractor(row),
                "legend": legend_extractor(row),
            }
            for n, row in df.iterrows()
        ]
        self.features.extend(loci)
        self.tracks.append(
            {
                "name": "reference",
                "separateFeaturesBy": "strand",
                "position": "inside",
                "thicknessRatio": 1,
                "dataType": "feature",
                "dataMethod": "source",
                "dataKeys": "reference",
            }
        )
        legend_item_names = [("CDS", "green"), ("tRNA", "magenta"), ("mRNA", "purple")]
        if self.with_highlighted_ogs:
            legend_item_names.extend(
                [("extracted yes", "magenta"), ("extracted no", "green")]
            )
        if self.with_amr:
            legend_item_names.extend([("AMR yes", "magenta"), ("AMR no", "green")])
        if self.with_vf:
            legend_item_names.extend([("VF yes", "magenta"), ("VF no", "green")])
        visible_items = {el["legend"] for el in loci}
        self.legend_items.extend(
            [
                {
                    "name": name,
                    "decoration": "arrow",
                    "swatchColor": color,
                    "visible": True if name in visible_items else False,
                }
                for name, color in legend_item_names
            ]
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
                "thicknessRatio": 2,
                "dataType": "feature",
                "dataMethod": "source",
                "dataKeys": label,
                "separateFeaturesBy": "none",
            },
        )
        self.legend_items.insert(
            0,
            {
                "name": label,
                "decoration": "score",
                "swatchColor": color,
            },
        )

    def add_gi_track(self, df):
        """
        Input df columns: gis_id, taxon_id, bioentry_id, cluster_id, start_pos, end_pos, length
        """
        loci = [
            {
                "name": f"GI{row.gis_id}",
                "source": "gi",
                "contig": str(row.bioentry_id),
                "start": int(row.start_pos),
                "stop": int(row.end_pos),
                "type": "genomic island",
                "legend": "Genomic island",
            }
            for n, row in df.iterrows()
        ]
        self.features.extend(loci)
        self.tracks.append(
            {
                "name": "Genomic island",
                "separateFeaturesBy": "none",
                "position": "outside",
                "dataKeys": "gi",
            }
        )
        self.legend_items.append(
            {
                "name": "Genomic island",
                "swatchColor": "gold",
            }
        )


class CircosView(BaseViewMixin, View):
    view_name = "circos"
    template = "chlamdb/circos.html"

    def dispatch(self, request, *args, **kwargs):
        self.form_class = make_circos_form(self.db)
        return super(CircosView, self).dispatch(request, *args, **kwargs)

    def get(self, request, *args, **kwargs):
        self.form = self.form_class()
        return render(request, self.template, self.get_context())

    def post(self, request, *args, **kwargs):
        self.form = self.form_class(request.POST)

        if self.form.is_valid():
            self.target_taxons = self.form.get_target_taxids()
            self.reference_taxon = self.form.get_ref_taxid()

            form_display = "inherit"
            self.highlighted_ogs = None
            if "highlighted_ogs" in self.form.data:
                # This is only set when coming from the OG extraction view.
                # As this is not a field supported in the form, we hide the form
                self.highlighted_ogs = self.form.data.getlist("highlighted_ogs")
                form_display = None

            self.prepare_circos_data()
            context = self.get_context(
                show_results=True,
                circos_json=json.dumps(self.data.to_json()),
                form_display=form_display,
                with_highlighted_ogs=self.data.with_highlighted_ogs,
                with_amr=self.data.with_amr,
                with_vf=self.data.with_vf,
                with_gi=self.data.with_gi,
            )
            return render(request, self.template, context)
        return render(request, self.template, self.get_context())

    def prepare_circos_data(self):
        self.data = CircosData(
            with_highlighted_ogs=self.highlighted_ogs is not None,
            with_amr=optional2status.get("amr", False),
            with_vf=optional2status.get("vf", False),
            with_gi=optional2status.get("gi", False),
        )
        # "bioentry_id", "accession" ,"length"
        df_bioentry = self.db.get_bioentry_list(
            self.reference_taxon, min_bioentry_length=1000
        )

        # "bioentry_id", "seqfeature_id", "start_pos", "end_pos", "strand"
        df_feature_location = self.db.get_features_location(
            self.reference_taxon,
            search_on="taxon_id",
            seq_term_names=["CDS", "rRNA", "tRNA"],
        ).set_index(["seqfeature_id"])
        # df of target genomes
        df_targets = self.db.get_proteins_info(
            ids=self.target_taxons,
            search_on="taxid",
            as_df=True,
            to_return=["locus_tag"],
            inc_non_CDS=False,
            inc_pseudo=False,
        )

        # retrieve n_orthologs of list of seqids
        seq_og = self.db.get_og_count(
            df_feature_location.index.to_list(), search_on="seqid"
        )
        count_all_genomes = self.db.get_og_count(
            seq_og["orthogroup"].to_list(), search_on="orthogroup"
        )
        n_genomes = self.db.get_number_of_genomes()
        orthogroup2frac_all = (
            count_all_genomes[count_all_genomes > 0].count(axis=1) / n_genomes
        )
        homologs_count = (
            df_feature_location.loc[df_feature_location.term_name == "CDS"]
            .join(seq_og)
            .reset_index()
            .set_index("orthogroup")
            .merge(
                orthogroup2frac_all.rename("value"), left_index=True, right_index=True
            )[["bioentry_id", "start_pos", "end_pos", "value"]]
        )

        # this query can be pretty slow
        df_identity = self.db.get_identity_closest_homolog(
            self.reference_taxon, self.target_taxons
        ).set_index(["target_taxid"])

        self.data.add_contigs_data(df_bioentry)

        # sort taxons by number of homologs (from mot similar to most dissmilar)
        target_taxon_n_homologs = (
            df_identity.groupby(["target_taxid"])
            .count()["seqfeature_id_1"]
            .sort_values(ascending=False)
        )
        # "bioentry_id", "seqfeature_id", "start_pos", "end_pos", "strand"
        # "seqfeature_id_1", "seqfeature_id_2", "identity", "target_taxid"
        # join on seqfeature id
        df_feature_location["gene"] = df_feature_location[
            "qualifier_value_gene"
        ].fillna("-")
        df_feature_location["gene_product"] = df_feature_location[
            "qualifier_value_product"
        ].fillna("-")
        df_feature_location["locus_tag"] = df_feature_location[
            "qualifier_value_locus_tag"
        ].fillna("-")

        if self.data.with_highlighted_ogs:
            df_feature_location["extracted"] = "no"
            df_genes = self.db.get_genes_from_og(
                self.highlighted_ogs,
                taxon_ids=[self.reference_taxon],
                terms=["locus_tag"],
            )
            df_feature_location.loc[df_genes["locus_tag"].index, "extracted"] = "yes"

        if self.data.with_amr:
            amrs = self.db.get_amr_hits_from_seqids(
                df_feature_location.index, columns=("seqid",)
            ).set_index("seqid")
            amrs["amr"] = "yes"
            df_feature_location = df_feature_location.join(amrs)
            df_feature_location.fillna({"amr": "no"}, inplace=True)

        if self.data.with_vf:
            vfs = self.db.vf.get_hits_from_seqids(
                df_feature_location.index, columns=("seqid",)
            ).set_index("seqid")
            vfs["vf"] = "yes"
            df_feature_location = df_feature_location.join(vfs)
            df_feature_location.fillna({"vf": "no"}, inplace=True)

        df_feature_location = df_feature_location.rename(
            columns={"locus_tag": "locus_ref"}
        )

        if self.data.with_gi:
            gis = self.db.gi.get_hits([self.reference_taxon])
            self.data.add_gi_track(gis)

        self.data.add_histogram_track(
            homologs_count,
            "Prevalence",
            "PaleVioletRed",
        )

        self.data.add_gene_track(df_feature_location)

        # iterate ordered list of target taxids, add track to circos
        if len(target_taxon_n_homologs.index) <= 10:
            taxon_colors = {
                el: tab10(i) for i, el in enumerate(target_taxon_n_homologs.index)
            }
        else:
            taxon_colors = {
                el: tab20(i % 20) for i, el in enumerate(target_taxon_n_homologs.index)
            }

        genome_descriptions = self.db.get_genomes_description(self.target_taxons)
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
            self.data.add_heatmap_data(
                df_combined,
                genome_descriptions.loc[target_taxon].description,
                mpl.colors.rgb2hex(taxon_colors[target_taxon]),
            )
        self.data.add_heatmap_track()
