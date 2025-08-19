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
                "legend": {"items": self.legend_items},
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
                "source": "genes",
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
                "name": "genes",
                "separateFeaturesBy": "strand",
                "position": "both",
                "thicknessRatio": 1,
                "dataType": "feature",
                "dataMethod": "source",
                "dataKeys": "genes",
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
                },
                "score": row.identity / 100,
                "name": row.locus_tag,
                "source": "orthogroups",
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
                "name": "orthogroups",
                "position": "inside",
                "thicknessRatio": 1,
                "dataType": "feature",
                "dataMethod": "source",
                "dataKeys": "orthogroups",
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
