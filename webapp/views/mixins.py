import pandas as pd
from django.conf import settings
from lib.db_utils import DB
from views.analysis_view_metadata import analysis_views_metadata
from views.object_type_metadata import AmrMetadata
from views.object_type_metadata import CogMetadata
from views.object_type_metadata import GiMetadata
from views.object_type_metadata import KoMetadata
from views.object_type_metadata import OrthogroupMetadata
from views.object_type_metadata import PfamMetadata
from views.object_type_metadata import VfMetadata
from views.object_type_metadata import my_locals
from views.utils import DataTableConfig
from views.utils import format_genome
from views.utils import format_genomic_island
from views.utils import format_hmm_url
from views.utils import format_ko_module
from views.utils import format_lst_to_html
from views.utils import format_refseqid_to_ncbi
from views.utils import get_genomes_data
from views.utils import page2title
from views.utils import safe_replace


class Transform:
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
            self.transform, axis=1
        )

    def maybe_add_accessory_col_for_query(self, columns, transformed):
        if not transformed or not columns or self.colname not in columns:
            return
        if self.accessory_colname not in columns:
            columns.append(self.accessory_colname)


class BaseViewMixin:
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
        return (
            page2title.get(getattr(self, "view_name", None))
            or page2title[f"{self.object_type}_comparison"]
        )

    @property
    def object_column(self):
        return self.object_type

    def colname_to_header(self, colname):
        return self._specific_colname_to_header_mapping.get(
            colname,
            self._base_colname_to_header_mapping.get(colname, colname.capitalize()),
        )

    @property
    def table_headers(self):
        return [self.colname_to_header(col) for col in self.table_data_accessors]

    def transform_data(self, descriptions):
        for transform in self.transforms:
            transform.maybe_apply(descriptions)
        return descriptions

    @property
    def transforms(self):
        return [Transform(self.object_column, self.format_entry, {"to_url": True})]

    @property
    def available_views(self):
        return [
            view_metadata(self.object_type, self.object_name_plural)
            for view_metadata in analysis_views_metadata
            if view_metadata.available_for(self.object_type)
        ]

    @property
    def metadata(self):
        if getattr(self, "_metadata_cls", None):
            return self._metadata_cls(self.object_type, self.object_name_plural)

    def get_context(self, **kwargs):
        context = {
            "page_title": self.page_title(),
        }
        if hasattr(self, "metadata"):
            context.update(
                {
                    "description": getattr(self.metadata, "description", ""),
                    "tab_name": getattr(self.metadata, "name", ""),
                }
            )
        if hasattr(self, "object_type"):
            context.update(
                {
                    "object_type": self.object_type,
                    "object_name": self.object_name_plural,
                    "object_name_plural": self.object_name_plural,
                    "available_views": self.available_views,
                }
            )
        if hasattr(self, "form"):
            context["form"] = self.form
        context.update(kwargs)
        return my_locals(context)

    def _format_column_headers_to_str(self, with_object_name=True):
        if with_object_name:
            first_col = 1
        else:
            first_col = 0
        return (
            ", ".join(self.table_headers[first_col:-1])
            + " and "
            + self.table_headers[-1]
        )


class AmrViewMixin(BaseViewMixin, AmrMetadata):
    object_column = "gene"

    _base_colname_to_header_mapping = {"seq_name": "Description", "hmm_id": "HMM"}

    table_data_accessors = [
        "gene",
        "seq_name",
        "scope",
        "type",
        "class",
        "subclass",
        "hmm_id",
    ]

    @property
    def transforms(self):
        return [
            Transform(self.object_column, self.format_entry, {"to_url": True}),
            Transform("hmm_id", format_hmm_url),
            Transform("class", safe_replace, {"args": ["/", " / "]}),
            Transform("subclass", safe_replace, {"args": ["/", " / "]}),
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

    def aggregate_amr_annotations(self, amr_annotations):
        gene_annot_counts = amr_annotations.gene.value_counts()
        for gene in gene_annot_counts[gene_annot_counts > 1].keys():
            amr_annotations[amr_annotations["gene"] == gene] = (
                self.aggregate_annotations_for_gene(gene, amr_annotations)
            )

    def aggregate_annotations_for_gene(self, gene, annotations):
        """
        entries for gene symbols or descriptions are not unique
        and there can therefore be more than one entry per gene.
        We concatenante them, it's the best we can do
        """
        rows = annotations[annotations.gene == gene]
        aggregated = []
        for col in rows.columns:
            aggregated.append(" || ".join(filter(None, rows[col].unique())))
        return aggregated


class CogViewMixin(BaseViewMixin, CogMetadata):
    _base_colname_to_header_mapping = {
        "cog": "ID",
        "function": "Function(s) cat.",
        "function_descr": "Function(s) descr.",
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

    def get_hit_descriptions(self, ids, transformed=True, only_cog_desc=True, **kwargs):
        descriptions = self.db.get_cog_summaries(
            ids, only_cog_desc=only_cog_desc, as_df=True
        )
        if transformed:
            descriptions = self.transform_data(descriptions)
        return descriptions

    @property
    def transforms(self):
        return [
            Transform(self.object_column, self.format_entry, {"to_url": True}),
            Transform(
                "function", self.format_function_descr, outcolname="function_descr"
            ),
        ]


class KoViewMixin(BaseViewMixin, KoMetadata):
    _base_colname_to_header_mapping = {
        "ko": "KO",
    }

    table_data_accessors = ["ko", "description", "pathways", "modules"]

    @property
    def get_hit_counts(self):
        return self.db.get_ko_hits

    def get_hit_descriptions(
        self,
        ids,
        transformed=True,
        extended_data=True,
        taxid_for_pathway_formatting=None,
        **kwargs,
    ):
        descriptions = self.db.get_ko_desc(ids, as_df=True)
        descriptions = descriptions.set_index(["ko"], drop=False)
        if extended_data:
            if not transformed:
                raise NotImplementedError(
                    "Can only use extended_data with transformed=True"
                )
            modules = self.db.get_ko_modules(ids, as_pandas=True)
            if not modules.empty:
                modules = modules.groupby("ko_id").apply(self.format_modules)
                descriptions = descriptions.merge(
                    modules.rename("modules"),
                    how="left",
                    left_index=True,
                    right_index=True,
                )
            else:
                descriptions["modules"] = None
            pathways = self.db.get_ko_pathways(ids, as_df=True)
            if not pathways.empty:
                pathways = pathways.groupby("ko").apply(
                    self.format_pathways, with_taxid=taxid_for_pathway_formatting
                )
                descriptions = descriptions.merge(
                    pathways.rename("pathways"),
                    how="left",
                    left_index=True,
                    right_index=True,
                )
            else:
                descriptions["pathways"] = None

        if transformed:
            descriptions = self.transform_data(descriptions)
            return descriptions.where(descriptions.notna(), "-")
        return descriptions

    @staticmethod
    def format_modules(modules, as_list=False):
        gen = (
            format_ko_module(row.module_id, row.desc) for i, row in modules.iterrows()
        )
        if as_list:
            return list(gen)
        return "<br>".join(gen)

    @staticmethod
    def format_pathways(pathways, with_taxid=None):
        if with_taxid is None:
            fmt_str = '<a href="/KEGG_mapp_ko/map{id:05d}">{descr}</a>'
        else:
            fmt_str = '<a href="/KEGG_mapp_ko/map{id:05d}/{taxid}">{descr}</a>'
        gen = (
            fmt_str.format(id=row.pathway, descr=row.description, taxid=with_taxid)
            for i, row in pathways.iterrows()
        )
        return "<br>".join(gen)


class PfamViewMixin(BaseViewMixin, PfamMetadata):
    _base_colname_to_header_mapping = {
        "pfam": "Domain ID",
        "def": "Description",
        "ttl_cnt": "nDomains",
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


class OrthogroupViewMixin(BaseViewMixin, OrthogroupMetadata):
    _base_colname_to_header_mapping = {"product": "Products", "gene": "Genes"}

    table_data_accessors = ["orthogroup", "products", "gene"]

    @property
    def get_hit_counts(self):
        return self.db.get_og_count

    def get_hit_descriptions(self, ids, transformed=True, **kwargs):
        genes = self.db.get_genes_from_og(ids)
        grouped = genes.groupby("orthogroup")
        genes = grouped["gene"].apply(list).apply(format_lst_to_html)
        products = grouped["product"].apply(list).apply(format_lst_to_html)
        descriptions = pd.DataFrame(
            {"orthogroup": ids}, index=pd.Index(ids, name="orthogroup")
        )

        descriptions = descriptions.merge(
            pd.DataFrame({"gene": genes, "products": products}),
            "left",
            left_index=True,
            right_index=True,
        )
        if transformed:
            descriptions = self.transform_data(descriptions)

        return descriptions


class VfViewMixin(BaseViewMixin, VfMetadata):
    object_column = "vf_gene_id"

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
        self.transform_category.maybe_add_accessory_col_for_query(cols, transformed)
        descriptions = self.db.vf.get_hit_descriptions(ids, columns=cols)
        if "vf_gene_id" in descriptions:
            descriptions = descriptions.set_index("vf_gene_id", drop=False)
        if transformed:
            descriptions = self.transform_data(descriptions)
        return descriptions

    @staticmethod
    def format_vfid(vfid):
        return (
            f'<a href="http://www.mgc.ac.cn/cgi-bin/VFs/vfs.cgi?VFID={vfid}"'
            f'target="_blank">{vfid}</a>'
        )

    @staticmethod
    def format_vf_category(category_and_id):
        category, vfcid = category_and_id
        return (
            f'<a href="http://www.mgc.ac.cn/cgi-bin/VFs/VFcategory.cgi?{vfcid}"'
            f'target="_blank">{category}</a>'
        )

    @property
    def transform_category(self):
        return TransformWithAccessoryColumn(
            "category", "vf_category_id", self.format_vf_category
        )

    @property
    def transforms(self):
        return [
            Transform(self.object_column, self.format_entry, {"to_url": True}),
            Transform("gb_accession", format_refseqid_to_ncbi),
            Transform("vfid", self.format_vfid),
            self.transform_category,
        ]


class GiViewMixin(BaseViewMixin, GiMetadata):
    object_column = "cluster_id"

    _base_colname_to_header_mapping = {
        "cluster_id": "Cluster ID",
        "gis_id": "ID",
        "bioentry_id": "Bioentry",
        "start_pos": "Start",
        "end_pos": "End",
    }

    table_data_accessors = ["cluster_id", "length"]

    def get_hit_descriptions(self, ids, transformed=True, **kwargs):
        descriptions = self.db.gi.get_hit_descriptions(ids)
        descriptions = descriptions.set_index("cluster_id", drop=False)
        if transformed:
            descriptions = self.transform_data(descriptions)
        return descriptions

    @property
    def get_hit_counts(self):
        return self.db.gi.get_hit_counts

    @property
    def get_hits(self):
        return self.db.gi.get_hits

    def get_gi_descriptions(self, ids, transformed=True, **kwargs):
        descriptions = self.db.gi.get_gi_descriptions(ids)
        descriptions["length"] = descriptions["end_pos"] - descriptions["start_pos"]
        descriptions = descriptions.set_index("gis_id", drop=False)
        if transformed:
            descriptions = self.transform_data(descriptions)
            descriptions["organism"] = descriptions["taxon_id"]
            descriptions.drop(columns=["taxon_id"], inplace=True)
        return descriptions

    @property
    def transforms(self):
        return [
            Transform(self.object_column, self.format_entry, {"to_url": True}),
            TransformWithAccessoryColumn("taxon_id", "organism", format_genome),
            Transform("gis_id", format_genomic_island),
        ]


class ComparisonViewMixin:
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
        "vf": VfViewMixin,
        "gic": GiViewMixin,
    }

    mixin = None

    def __getattr__(self, attrname):
        return getattr(self.mixin, attrname)

    def dispatch(self, request, comp_type, *args, **kwargs):
        self.object_type = comp_type
        self.mixin = self.type2mixin[self.object_type]()
        # Fix the metadata property
        self.mixin._metadata_cls = getattr(self, "_metadata_cls", None)
        return super(ComparisonViewMixin, self).dispatch(request, *args, **kwargs)


class GenomesTableMixin:
    _genome_table_help = (
        "This table contains the list of genomes included in the {} and a "
        "summary of their content. <br> Clicking on the genome name "
        "(first column) a second table with the protein content is displayed. "
        "It shows to which contig each protein belongs and provides the link "
        "to the locus tags. <br> Fasta and gbk files can be downloaded."
    )

    @property
    def genome_table_help(self):
        return self._genome_table_help.format(self.genome_source_object)

    def get_genomes_table(self, taxids=None):
        genomes_data = get_genomes_data(self.db, taxids=taxids)

        filenames_tax_id = self.db.get_filenames_to_taxon_id()
        filenames_tax_id_db = pd.DataFrame.from_dict(list(filenames_tax_id.items()))
        filenames_tax_id_db.columns = ["filename", "taxon_id"]
        filenames_tax_id_db.index = list(filenames_tax_id_db["taxon_id"])
        filenames_list = list(filenames_tax_id_db["filename"])

        path_template = settings.BLAST_DB_PATH + "/{ext}/{filename}.{ext}"
        link_template = '<a href="{}"> .{{ext}} </a>'.format(path_template)
        for ext in ["faa", "fna", "ffn", "gbk"]:
            filenames_tax_id_db[f"path_to_{ext}"] = [
                link_template.format(filename=filename, ext=ext)
                for filename in filenames_list
            ]

        filenames_tax_id_db = filenames_tax_id_db[
            ["path_to_faa", "path_to_fna", "path_to_ffn", "path_to_gbk"]
        ]
        genomes_data = genomes_data.join(filenames_tax_id_db, on="taxon_id")

        genomes_data["accession"] = genomes_data[["id", "description"]].apply(
            format_genome, axis=1
        )

        data_table_header = [
            "Name",
            "GC %",
            "N proteins",
            "N contigs",
            "Size (Mbp)",
            "Coding %",
            "Has plasmid(s)",
            "faa seq",
            "fna seq",
            "ffn seq",
            "gbk file",
        ]

        table_data_accessors = [
            "accession",
            "gc",
            "n_prot",
            "n_contigs",
            "length",
            "coding_density",
            "has_plasmid",
            "path_to_faa",
            "path_to_fna",
            "path_to_ffn",
            "path_to_gbk",
        ]

        table_data = genomes_data[table_data_accessors]

        return {
            "table_data": table_data,
            "table_headers": data_table_header,
            "data_table_config": DataTableConfig(),
            "table_data_accessors": table_data_accessors,
            "table_help": self.genome_table_help,
        }
