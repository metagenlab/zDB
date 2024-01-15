from views.utils import page2title


class ComparisonViewMixin():

    type2objname = {
        "cog": "COGs",
        "pfam": "Pfam domains",
        "ko": "Kegg Orthologs",
        "orthogroup": "Orthologous groups",
        "amr": "AMR"
    }

    @property
    def compared_obj_name(self):
        return self.type2objname[self.comp_type]

    @property
    def page_title(self):
        return page2title[self.view_name]


class AmrAnnotationsMixin():

    def aggregate_amr_annotations(self, amr_annotations):
        gene_annot_counts = amr_annotations.gene.value_counts()
        for gene in gene_annot_counts[gene_annot_counts > 1].keys():
            amr_annotations[amr_annotations["gene"] == gene] = \
                self.aggregate_annotations_for_gene(gene, amr_annotations)

    def aggregate_annotations_for_gene(self, gene, annotations):
        """
        entries for gene symbols or descriptions are not unique
        and there can therefore be more than one entry per gene.
        We concatenante them, it's the best we can do
        """
        rows = annotations[annotations.gene == gene]
        aggregated = []
        for col in rows.columns:
            aggregated.append(" || ".join(rows[col].unique()))
        return aggregated
