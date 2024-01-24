from django.conf import settings
from lib.db_utils import DB

from views.utils import (format_amr, format_cog, format_hmm_url, format_ko,
                         format_orthogroup, format_pfam, page2title)


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

    @property
    def object_column(self):
        return self.object_type


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
    object_column = "gene"

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
        descriptions = descriptions.drop_duplicates(subset=["gene"])
        descriptions.set_index("gene", drop=False)
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


class CogViewMixin(BaseViewMixin):

    object_type = "cog"
    object_name = "COG entry"
    object_name_plural = "COG entries"
    object_name_singular_or_plural = "COG entry(ies)"

    colname_to_header = {
        "cog": "ID",
        "function": "Function(s)",
        "description": "Description",
    }

    @property
    def get_hit_counts(self):
        return self.db.get_cog_hits

    def get_hit_descriptions(self, ids, transformed=True):
        descriptions = self.db.get_cog_summaries(
            ids, only_cog_desc=True, as_df=True)
        if transformed:
            cog_func = self.db.get_cog_code_description()
            descriptions["function"] = descriptions["function"].apply(lambda func: "<br>".join((cog_func[code] for code in func)))
        return descriptions

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_cog(entry, as_url=to_url)


class KoViewMixin(BaseViewMixin):

    object_type = "ko"
    object_name = "Kegg Ortholog"

    colname_to_header = {
        "ko": "KO",
        "description": "Description",
    }

    @property
    def get_hit_counts(self):
        return self.db.get_ko_hits

    def get_hit_descriptions(self, ids, transformed=True):
        descriptions = self.db.get_ko_desc(ids, as_df=True)
        descriptions.set_index(["ko"], drop=False)
        if transformed:
            descriptions["ko"] = descriptions["ko"].apply(self.format_entry)
        return descriptions

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_ko(entry, as_url=to_url)


class PfamViewMixin(BaseViewMixin):

    object_type = "pfam"
    object_name = "Pfam domain"

    colname_to_header = {
        "pfam": "PFAM",
        "def": "Description",
    }

    @property
    def get_hit_counts(self):
        return self.db.get_pfam_hits

    def get_hit_descriptions(self, ids, transformed=True):
        descriptions = self.db.get_pfam_def(ids)
        return descriptions

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_pfam(entry, to_url=to_url)


class OrthogroupViewMixin(BaseViewMixin):

    object_type = "orthogroup"
    object_name = "Orthologous group"

    @property
    def get_hit_counts(self):
        return self.db.get_og_count

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_orthogroup(entry, to_url=to_url)
