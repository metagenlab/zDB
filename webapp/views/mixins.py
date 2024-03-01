import pandas as pd
from django.conf import settings
from lib.db_utils import DB

from views.analysis_view_metadata import analysis_views_metadata
from views.utils import (format_amr, format_cog, format_hmm_url, format_ko,
                         format_lst_to_html, format_orthogroup, format_pfam,
                         format_refseqid_to_ncbi, my_locals, page2title,
                         safe_replace)


class Transform():

    def __init__(self, colname, transform, kwargs=None, outcolname=None):
        self.colname = colname
        self.transform = transform
        self.kwargs = kwargs or {}
        self.outcolname = outcolname or colname

    def should_apply(self, df):
        return self.colname in df

    def apply(self, df):
        df[self.outcolname] = df[self.colname].apply(self.transform, **self.kwargs)

    def maybe_apply(self, df):
        if self.should_apply(df):
            self.apply(df)


class TransformWithAccessoryColumn(Transform):

    def __init__(self, colname, accessory_colname, transform):
        self.colname = colname
        self.accessory_colname = accessory_colname
        self.transform = transform

    def should_apply(self, df):
        return self.colname in df and self.accessory_colname in df

    def apply(self, df):
        df[self.colname] = df[[self.colname, self.accessory_colname]].apply(
            self.transform, axis=1)

    def maybe_add_accessory_col_for_query(self, columns, transformed):
        if not transformed or not columns or self.colname not in columns:
            return
        if self.accessory_colname not in columns:
            columns.append(self.accessory_colname)


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
        return page2title.get(getattr(self, "view_name", None),
                              page2title[f"{self.object_type}_comparison"])

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

    def transform_data(self, descriptions):
        for transform in self.transforms:
            transform.maybe_apply(descriptions)
        return descriptions

    @property
    def transforms(self):
        return [
            Transform(self.object_column, self.format_entry, {"to_url": True})]

    @property
    def available_views(self):
        return [view_metadata(self.object_type, self.object_name_plural)
                for view_metadata in analysis_views_metadata
                if view_metadata.available_for(self.object_type)]

    def get_context(self, **kwargs):
        context = {
            "page_title": self.page_title(),
            "object_type": self.object_type,
            "object_name": self.object_name_plural,
            "object_name_plural": self.object_name_plural,
            "available_views": self.available_views,
        }
        if hasattr(self, "form"):
            context["form"] = self.form
        context.update(kwargs)
        return my_locals(context)


class AmrViewMixin(BaseViewMixin):

    object_type = "amr"
    object_name = "AMR gene"
    object_column = "gene"

    _base_colname_to_header_mapping = {
        "seq_name": "Description",
        "hmm_id": "HMM"
    }

    table_data_accessors = ["gene", "seq_name", "scope", "type", "class",
                            "subclass", "hmm_id"]

    @property
    def transforms(self):
        return [
            Transform(self.object_column, self.format_entry, {"to_url": True}),
            Transform("hmm_id", format_hmm_url),
            Transform("class", safe_replace, {"args": ["/", " / "]}),
            Transform("subclass", safe_replace, {"args": ["/", " / "]})
        ]

    @property
    def get_hit_counts(self):
        return self.db.get_amr_hit_counts

    def get_hit_descriptions(self, ids, transformed=True, **kwargs):
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

    table_data_accessors = ["cog", "function_descr", "description"]

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

    def get_hit_descriptions(self, ids, transformed=True, **kwargs):
        descriptions = self.db.get_cog_summaries(
            ids, only_cog_desc=True, as_df=True)
        if transformed:
            descriptions = self.transform_data(descriptions)
        return descriptions

    @property
    def transforms(self):
        return [
            Transform(self.object_column, self.format_entry, {"to_url": True}),
            Transform("function", self.format_function_descr,
                      outcolname="function_descr"),
        ]

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_cog(entry, as_url=to_url)


class KoViewMixin(BaseViewMixin):

    object_type = "ko"
    object_name = "Kegg Ortholog"

    _base_colname_to_header_mapping = {
        "ko": "KO",
    }

    table_data_accessors = ["ko", "description", "modules", "pathways"]

    @property
    def get_hit_counts(self):
        return self.db.get_ko_hits

    def get_hit_descriptions(self, ids, transformed=True, **kwargs):
        descriptions = self.db.get_ko_desc(ids, as_df=True)
        descriptions = descriptions.set_index(["ko"], drop=False)
        if transformed:
            descriptions = self.transform_data(descriptions)
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

    table_data_accessors = ["pfam", "def"]

    @property
    def get_hit_counts(self):
        return self.db.get_pfam_hits

    def get_hit_descriptions(self, ids, transformed=True, **kwargs):
        descriptions = self.db.get_pfam_def(ids, **kwargs)
        if transformed:
            descriptions = self.transform_data(descriptions)
        return descriptions

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_pfam(entry, to_url=to_url)


class OrthogroupViewMixin(BaseViewMixin):

    object_type = "orthogroup"
    object_name = "Orthologous group"

    _base_colname_to_header_mapping = {
        "product": "Products",
        "gene": "Genes"
    }

    table_data_accessors = ["orthogroup", "product", "gene"]

    @property
    def get_hit_counts(self):
        return self.db.get_og_count

    def get_hit_descriptions(self, ids, transformed=True, **kwargs):
        genes = self.db.get_genes_from_og(ids)
        grouped = genes.groupby("orthogroup")
        genes = grouped["gene"].apply(list).apply(format_lst_to_html)
        products = grouped["product"].apply(list).apply(format_lst_to_html)
        descriptions = pd.DataFrame({"orthogroup": ids},
                                    index=pd.Index(ids, name="orthogroup"))

        descriptions = descriptions.merge(
            pd.DataFrame({"gene": genes, "product": products}),
            "left", left_index=True, right_index=True)
        if transformed:
            descriptions = self.transform_data(descriptions)

        return descriptions

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_orthogroup(entry, to_url=to_url)


class VfViewMixin(BaseViewMixin):

    object_type = "vf"
    object_column = "vf_gene_id"
    object_name = "Virulence factor"

    _base_colname_to_header_mapping = {
        "vf_gene_id": "VF gene ID",
        "prot_name": "Protein",
        "vfid": "VF ID",
        "gb_accession": "Protein ID",

    }

    table_data_accessors = ["vf_gene_id", "prot_name", "vfid", "category"]

    @property
    def get_hit_counts(self):
        return self.db.vf.get_hit_counts

    @property
    def get_hits(self):
        return self.db.vf.get_hits

    def get_hit_descriptions(self, ids, transformed=True, **kwargs):
        cols = kwargs.get("columns")
        self.transform_category.maybe_add_accessory_col_for_query(
            cols, transformed)
        descriptions = self.db.vf.get_hit_descriptions(ids, columns=cols)
        if "vf_gene_id" in descriptions:
            descriptions = descriptions.set_index("vf_gene_id", drop=False)
        if transformed:
            descriptions = self.transform_data(descriptions)
        return descriptions

    @staticmethod
    def format_entry(entry, to_url=False):
        """Will point to the details page as soon as it exists
        """
        if to_url:
            return f"<a href=\"/fam_vf/{entry}\">{entry}</a>"
        return entry

    @staticmethod
    def format_vfid(vfid):
        return f'<a href="http://www.mgc.ac.cn/cgi-bin/VFs/vfs.cgi?VFID={vfid}"'\
               f'target="_blank">{vfid}</a>'

    @staticmethod
    def format_vf_category(category_and_id):
        category, vfcid = category_and_id
        return f'<a href="http://www.mgc.ac.cn/cgi-bin/VFs/VFcategory.cgi?{vfcid}"'\
               f'target="_blank">{category}</a>'

    @property
    def transform_category(self):
        return TransformWithAccessoryColumn(
            "category", "vf_category_id", self.format_vf_category)

    @property
    def transforms(self):
        return [
            Transform(self.object_column, self.format_entry, {"to_url": True}),
            Transform("gb_accession", format_refseqid_to_ncbi),
            Transform("vfid", self.format_vfid),
            self.transform_category
        ]


class ComparisonViewMixin():
    """This class is somewhat of a hack to get pseudo inheritance
    from the correct mixin for views that get the object_type as
    parameter.
    """
    type2mixin = {
        "cog": CogViewMixin,
        "pfam": PfamViewMixin,
        "ko": KoViewMixin,
        "orthogroup": OrthogroupViewMixin,
        "amr": AmrViewMixin,
        "vf": VfViewMixin
    }

    mixin = None

    def __getattr__(self, attrname):
        return getattr(self.mixin, attrname)

    def dispatch(self, request, comp_type, *args, **kwargs):
        self.object_type = comp_type
        self.mixin = self.type2mixin[self.object_type]()
        return super(ComparisonViewMixin, self).dispatch(request, *args, **kwargs)
