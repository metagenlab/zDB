from chlamdb.forms import make_metabo_from
from django.shortcuts import render
from django.views import View

from views.mixins import (AmrViewMixin, CogViewMixin, KoViewMixin,
                          OrthogroupViewMixin, PfamViewMixin, VfViewMixin)
from views.utils import (format_cog, format_ko, format_lst_to_html,
                         format_orthogroup, my_locals)


class TabularComparisonViewBase(View):

    template = 'chlamdb/tabular_comparison.html'
    hist_colour_index_shift = 0
    tab_name = "comp"
    table_headers = None

    def dispatch(self, request, *args, **kwargs):
        self.comp_metabo_form = self.make_metabo_from()
        self.show_comparison_table = False
        self._hash_to_taxon_dict = None
        return super(TabularComparisonViewBase, self).dispatch(request, *args, **kwargs)

    def make_metabo_from(self):
        return make_metabo_from(self.db)

    def get(self, request, *args, **kwargs):
        self.form = self.comp_metabo_form()
        return render(request, self.template, self.context)

    def post(self, request, *args, **kwargs):
        self.form = self.comp_metabo_form(request.POST)
        if self.form.is_valid():
            self.show_comparison_table = True
            self.targets = self.form.get_choices()
            if "comp_type" in self.form.cleaned_data:
                self.comp_type = self.form.cleaned_data["comp_type"]
            self.set_table_data()

        return render(request, self.template, self.context)

    @property
    def view_name(self):
        return f"{self.object_type}_comparison"

    @property
    def context(self):
        context = {
            "page_title": self.page_title,
            "form_title": self.form_title,
            "form_help": self.form_help,
            "form": self.form,
            "show_comparison_table": self.show_comparison_table,
            "tab_name": self.tab_name,
            "view_name": self.view_name,
            "object_type": self.object_type,
            }
        if self.show_comparison_table:
            context["table_headers"] = self.table_headers
            context["table_rows"] = self.table_rows
            context["table_title"] = self.table_title
            context["table_help"] = self.table_help
            context["first_coloured_row"] = self.first_coloured_row
            context["n_data_columns"] = len(self.base_info_headers) + \
                len(self.targets)
            context["hist_colour_index_shift"] = self.hist_colour_index_shift
        return my_locals(context)

    @property
    def hash_to_taxon_dict(self):
        if self._hash_to_taxon_dict is None:
            self._hash_to_taxon_dict = self.db.get_genomes_description().description.to_dict()
        return self._hash_to_taxon_dict

    def set_table_data(self):
        self.n_selected = len(self.targets)

        taxon_list = [self.hash_to_taxon_dict[target_id]
                      for target_id in self.targets]
        self.table_headers = self.base_info_headers + taxon_list

        self.table_rows = self.get_table_rows()
        self.n_rows = len(self.table_rows)

    def get_table_rows(self):
        """This method should return an iteratable object with each iteration
        yielding a row of the table. Apart from the data, rows can contain
        values used to color the cells, these are stored in the rows after
        the data and will be matched to the corresponding data value with a
        shift in index (row[i] colored by row[i + hist_colour_index_shift]).
        """
        raise NotImplementedError

    @property
    def table_title(self):
        return "Number of {} present at least once in 1 of the {} selected "\
               "genomes: <strong>{}</strong>".format(self.object_name_plural,
                                                     self.n_selected,
                                                     self.n_rows)

    @property
    def form_title(self):
        return "Compare the distribution of shared {}.".format(
            self.object_name_plural)

    @property
    def form_help(self):
        return "Compare the size of the {} shared by selected genomes (targets).".format(
            self.object_name_plural)

    @property
    def first_coloured_row(self):
        return len(self.base_info_headers)

    @property
    def base_info_headers(self):
        return [self.colname_to_header(colname)
                for colname in self.base_info_accessors]


class PfamComparisonView(TabularComparisonViewBase, PfamViewMixin):

    base_info_accessors = ["pfam", "def", "ttl_cnt"]

    table_help = """
    The ouput table contains the list of shared Pfam domains and the number of
    times each of them was identified in the selected genomes.
    <br>nDomain: total number of occurence of this domain in the
    complete database.
    <br>Click on Pfam accession to get detailed phylogenetic profile of the
    corresponding Pfam entry.
    """

    def get_table_rows(self):
        pfam_hits = self.get_hit_counts(ids=self.targets)
        pfam_defs = self.get_hit_descriptions(pfam_hits.index.tolist(),
                                              add_ttl_count=True)

        table_rows = []
        for key, values in pfam_hits.iterrows():
            entry_infos = pfam_defs.loc[key]
            base_infos = [entry_infos[accessor]
                          for accessor in self.base_info_accessors]
            table_rows.append(base_infos + values.values.tolist())

        return table_rows


class CogComparisonView(TabularComparisonViewBase, CogViewMixin):

    base_info_headers = ["COG accession", "Description", "# complete DB", "# genomes"]

    table_help = """
    The ouput table contains the list of COG annotated in selected genomes and
    the number of times each of them was identified in each genome.
    <br>Click on COG accession to get detailed phylogenetic profile of the
    corresponding COG entry.
    """

    def get_table_rows(self):
        cog_hits = self.db.get_cog_hits(
            ids=self.targets, search_on="taxid")
        # retrieve entry list
        cog_all = self.db.get_cog_hits(
            ids=list(self.hash_to_taxon_dict.keys()),
            search_on="taxid")

        # count frequency and n genomes
        cog_count = cog_all.sum(axis=1)
        cog_freq = cog_all[cog_all > 0].count(axis=1)
        cog_freq = cog_freq.loc[cog_hits.index]
        cog_count = cog_count.loc[cog_hits.index]
        # retrieve annotations
        cogs_summaries = self.db.get_cog_summaries(
            cog_hits.index.tolist(), as_df=True, only_cog_desc=True)

        cog2description = cogs_summaries.description.to_dict()
        cog_hits["accession"] = [format_cog(cog, as_url=True)
                                 for cog in cog_hits.index]
        cog_hits["description"] = [cog2description[cog] if cog in cog2description else '-'
                                   for cog in cog_hits.index]

        # combine into df
        combined_df = cog_hits.merge(
            cog_count.rename('count'),
            left_index=True,
            right_index=True).merge(
                cog_freq.rename('freq'),
                left_index=True,
                right_index=True).sort_values(["count"], ascending=False)

        cols = combined_df.columns.to_list()
        ordered_cols = cols[self.n_selected:] + cols[:self.n_selected]
        return combined_df[ordered_cols].values

    @property
    def first_coloured_row(self):
        return 4


class OrthogroupComparisonView(TabularComparisonViewBase, OrthogroupViewMixin):

    base_info_headers = ["Orthogroup", "Annotaion"]

    table_help = """
    The ouput table contains the number of homologs in the shared orthogroups
    of the selected genomes. Interesting for comparing the size of orthogroups
    within genomes.
    <br> Homolog counts can be reordrered by clicking on column headers.<br>
    <br>Click on Orthologous group to get all the homologs identified in the
    database and the phylogenetic profile.
    """

    @property
    def view_name(self):
        return "orthogroup_comparison"

    def get_table_rows(self):
        og_count = self.db.get_og_count(self.targets)
        annotations = self.db.get_genes_from_og(
            orthogroups=og_count.index.tolist(),
            taxon_ids=self.targets, terms=["product"])

        products = annotations.groupby("orthogroup")["product"].apply(list)

        og_data = []
        for og, items in og_count.iterrows():
            row = [format_orthogroup(og, to_url=True)]
            if og in products.index:
                row.append(format_lst_to_html(products.loc[og]))
            else:
                row.append("-")
            row.extend(items)

            og_data.append(row)

        return og_data


class KoComparisonView(TabularComparisonViewBase, KoViewMixin):

    table_help = """
    The ouput table contains the number of homologs in the shared Kegg
    Orthologs of the selected genomes and the total number of homologs
    of each Kegg Orthologs identified in the whole collection. Interesting
    for comparing the size of Kegg Orthologs within genomes.
    <br> Homolog counts can be reordrered by clicking on column headers.<br>
    <br>Click on the Ko entry and list the Ko modules and pathways of which it
    is part.
    """

    base_info_headers = ["KO", "Annot", "tot"]

    def get_table_rows(self):
        hits = self.db.get_ko_count(self.targets).unstack(level=0, fill_value=0)
        hits.columns = [col for col in hits["count"].columns.values]

        ko2annot = self.db.get_ko_desc(hits.index.tolist())
        df_ttl = self.db.get_ko_count(hits.index.tolist(), search_on="ko_id")
        ko2total_count = df_ttl.groupby("KO").sum()["count"].to_dict()
        table_rows = []
        for key, values in hits.iterrows():
            row = [format_ko(key, as_url=True), ko2annot[key], ko2total_count[key]]
            row.extend(values.values)
            table_rows.append(row)
        return table_rows


class AmrComparisonView(TabularComparisonViewBase, AmrViewMixin):

    _table_help = """
    The ouput table contains the number of times a given AMR {} appears
    in the selected genomes, color coded according to the quality
    (coverage*identity) of the best hit for that genome.<br>
    <br> Counts can be reordrered by clicking on column headers.<br>
    """

    _scope_hint = """
    <br> Note that genes are split into "core" and "plus" scopes, where
    "core" proteins are expected to have an effect on resistance while
    "plus" proteins are included with a less stringent criteria.<br>
    """

    type_choices = (("gene", "Gene"),
                    ("class", "Class"),
                    ("subclass", "Subclass"))

    def make_metabo_from(self):
        return make_metabo_from(self.db, type_choices=self.type_choices)

    @property
    def base_info_accessors(self):
        if self.comp_type == "gene":
            return ["gene", "scope", "class", "subclass", "seq_name", "hmm_id"]
        elif self.comp_type == "subclass":
            return ["subclass", "class"]
        elif self.comp_type == "class":
            return ["class"]

    @property
    def table_help(self):
        if self.comp_type == "gene":
            return self._table_help.format(self.comp_type) + self._scope_hint
        else:
            return self._table_help.format(self.comp_type)

    def get_table_rows(self):
        hits = self.db.get_amr_hits_from_taxonids(self.targets)
        hits = self.transform_data(hits)

        table_rows = []
        hits["quality"] = hits["coverage"] * hits["identity"] / 10000
        for groupid, data in hits.groupby(self.comp_type):
            row_data = data.iloc[0]
            row = [row_data[key] or "-" for key in self.base_info_accessors]
            taxonids = data["bioentry.taxon_id"]
            values = [len(taxonids[taxonids == target_id])
                      for target_id in self.targets]
            colours = [data[taxonids == target_id]["quality"].max() if value
                       else 0 for value, target_id in zip(values, self.targets)]
            row.extend(values)
            row.extend(colours)
            table_rows.append(row)
        return table_rows

    @property
    def hist_colour_index_shift(self):
        return len(self.targets)


class VfComparisonView(TabularComparisonViewBase, VfViewMixin):

    _table_help = """
    The ouput table contains the number of times a given {} appears
    in the selected genomes, color coded according to the e-value
    of the best hit for that genome.<br>
    <br> Counts can be reordrered by clicking on column headers.<br>
    """

    type_choices = (("vf_gene_id", "VF Gene"),
                    ("vfid", "VF"),
                    ("category", "VF Category"))

    def make_metabo_from(self):
        return make_metabo_from(self.db, type_choices=self.type_choices)

    @property
    def base_info_accessors(self):
        if self.comp_type == "vf_gene_id":
            return ["vf_gene_id", "prot_name", "vfid", "category"]
        elif self.comp_type == "vfid":
            return ["vfid", "category"]
        elif self.comp_type == "category":
            return ["category"]

    @property
    def table_help(self):
        return self._table_help.format(self.comp_type)

    def get_table_rows(self):
        hits = self.get_hits(self.targets)
        hits = self.transform_data(hits)
        hits = hits.where(hits.notna(), "-")
        table_rows = []
        for groupid, data in hits.groupby(self.comp_type):
            row = data.iloc[0][self.base_info_accessors].tolist()
            taxonids = data["taxon_id"]
            values = [len(taxonids[taxonids == target_id])
                      for target_id in self.targets]
            colours = [data[taxonids == target_id]["evalue"].max() if value
                       else 0 for value, target_id in zip(values, self.targets)]
            row.extend(values)
            row.extend(colours)
            table_rows.append(row)
        return table_rows

    @property
    def hist_colour_index_shift(self):
        return len(self.targets)
