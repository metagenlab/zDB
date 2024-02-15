from django.conf import settings
from lib.db_utils import DB

from views.utils import (format_amr, format_cog, format_hmm_url, format_ko,
                         format_lst_to_html, format_orthogroup, format_pfam,
                         page2title, safe_replace)


class BaseViewMixin():

    _db = None
    _base_colname_to_header_mapping = {}
    _specific_colname_to_header_mapping = {}

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

    def colname_to_header(self, colname):
        return self._specific_colname_to_header_mapping.get(
            colname, self._base_colname_to_header_mapping.get(
                colname, colname.capitalize()))

    @property
    def table_headers(self):
        return [self.colname_to_header(col)
                for col in self.table_data_accessors]


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

    _base_colname_to_header_mapping = {
        "seq_name": "Description",
        "hmm_id": "HMM"
    }

    transforms = {
        "gene": (format_amr, {"to_url": True}),
        "hmm_id": (format_hmm_url, {}),
        "class": (safe_replace, {"args": ["/", " / "]}),
        "subclass": (safe_replace, {"args": ["/", " / "]})
        }

    @property
    def get_hit_counts(self):
        return self.db.get_amr_hit_counts

    def transform_data(self, descriptions):
        for colname, (transform, kwargs) in self.transforms.items():
            descriptions[colname] = descriptions[colname].apply(transform,
                                                                **kwargs)
        return descriptions

    def get_hit_descriptions(self, ids, transformed=True):
        descriptions = self.db.get_amr_descriptions(ids)
        self.aggregate_amr_annotations(descriptions)
        descriptions = descriptions.drop_duplicates(subset=["gene"])
        descriptions = descriptions.set_index("gene", drop=False)
        if transformed:
            descriptions = self.transform_data(descriptions)
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

    _base_colname_to_header_mapping = {
        "cog": "ID",
        "function": "Function(s)",
        "function_descr": "Function(s)",
    }

    _cog_code_descriptions = None

    @property
    def get_hit_counts(self):
        return self.db.get_cog_hits

    @property
    def cog_code_descriptions(self):
        if not self._cog_code_descriptions:
            self._cog_code_descriptions = self.db.get_cog_code_description()
        return self._cog_code_descriptions

    def format_function_descr(self, func):
        descr = [f"{self.cog_code_descriptions[abbr]} ({abbr})" for abbr in func]
        return format_lst_to_html(descr, False)

    def get_hit_descriptions(self, ids, transformed=True):
        descriptions = self.db.get_cog_summaries(
            ids, only_cog_desc=True, as_df=True)
        if transformed:
            descriptions["cog"] = descriptions["cog"].apply(
                self.format_entry, to_url=True)
            descriptions["function_descr"] = descriptions.function.apply(
                self.format_function_descr)
        return descriptions

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_cog(entry, as_url=to_url)


class KoViewMixin(BaseViewMixin):

    object_type = "ko"
    object_name = "Kegg Ortholog"

    _base_colname_to_header_mapping = {
        "ko": "KO",
    }

    @property
    def get_hit_counts(self):
        return self.db.get_ko_hits

    def get_hit_descriptions(self, ids, transformed=True):
        descriptions = self.db.get_ko_desc(ids, as_df=True)
        descriptions = descriptions.set_index(["ko"], drop=False)
        if transformed:
            descriptions["ko"] = descriptions["ko"].apply(self.format_entry,
                                                          to_url=True)
        return descriptions

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_ko(entry, as_url=to_url)


class PfamViewMixin(BaseViewMixin):

    object_type = "pfam"
    object_name = "Pfam domain"

    _base_colname_to_header_mapping = {
        "pfam": "Domain ID",
        "def": "Description",
        "ttl_cnt": "nDomains"
    }

    @property
    def get_hit_counts(self):
        return self.db.get_pfam_hits

    def get_hit_descriptions(self, ids, transformed=True, **kwargs):
        descriptions = self.db.get_pfam_def(ids, **kwargs)
        if transformed:
            descriptions["pfam"] = descriptions["pfam"].apply(
                self.format_entry, to_url=True)
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

    def get_hit_descriptions(self, ids, transformed=True):
        descriptions = self.db.get_genes_from_og(ids)
        if transformed:
            descriptions["orthogroup"] = descriptions["orthogroup"].apply(
                self.format_entry)
        return descriptions

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_orthogroup(entry, to_url=to_url)


class VfViewMixin(BaseViewMixin):

    object_type = "vf"
    object_name = "Virulence factor"

    @property
    def get_hit_counts(self):
        return self.db.vf.get_hit_counts

    def get_hit_descriptions(self, ids, transformed=True):
        descriptions = self.db.vf.get_hit_descriptions(ids)
        if transformed:
            descriptions["vf_id"] = descriptions["vf_id"].apply(self.format_entry)
        return descriptions

    @staticmethod
    def format_entry(entry, to_url=False):
        """Will point to the details page as soon as it exists
        """
        return entry
