from django.conf import settings
from lib.db_utils import DB

from views.utils import format_amr, format_hmm_url, page2title


class BaseViewMixin():

    _db = None

    @property
    def db(self):
        if self._db is None:
            biodb_path = settings.BIODB_DB_PATH
            self._db = DB.load_db_from_name(biodb_path)
        return self._db

    def page_title(self):
        return page2title[self.view_name]

    @property
    def object_name_plural(self):
        return f"{self.object_name}s"

    @property
    def object_name_singular_or_plural(self):
        return f"{self.object_name}(s)"


class ComparisonViewMixin(BaseViewMixin):

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


class AmrViewMixin(BaseViewMixin):

    object_type = "amr"
    object_name = "AMR gene"

    colname_to_header = {
        "gene": "Gene",
        "seq_name": "Description",
        "scope": "Scope",
        "type": "Type",
        "class": "Class",
        "subclass": "Subclass",
        "hmm_id": "HMM"
    }

    transforms = {
        "gene": (format_amr, {"to_url": True}),
        "hmm_id": (format_hmm_url, {}),
        }

    @property
    def get_hit_counts(self):
        return self.db.get_amr_hit_counts

    def get_hit_descriptions(self, ids, transformed=True):
        descriptions = self.db.get_amr_descriptions(ids)
        self.aggregate_amr_annotations(descriptions)
        if transformed:
            for colname, (transform, kwargs) in self.transforms.items():
                descriptions[colname] = descriptions[colname].apply(transform,
                                                                    **kwargs)
        return descriptions

    @staticmethod
    def format_entry(gene, to_url=False):
        return format_amr(gene, to_url=to_url)

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
