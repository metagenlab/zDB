
from chlamdb.forms import make_venn_from
from django.conf import settings
from django.shortcuts import render
from django.views import View
from lib.db_utils import DB

from views.mixins import ComparisonViewMixin
from views.utils import (format_cog, format_ko, format_lst_to_html,
                         format_orthogroup, format_pfam, my_locals, page2title,
                         to_s)


def escape_quotes(unsafe):
    return unsafe.replace("\"", "\\\"")


class VennBaseView(View, ComparisonViewMixin):

    template = 'chlamdb/venn_generic.html'

    @property
    def view_name(self):
        return f"venn_{self.comp_type}"

    def dispatch(self, request, *args, **kwargs):
        biodb_path = settings.BIODB_DB_PATH
        self.db = DB.load_db_from_name(biodb_path)
        self.form_class = make_venn_from(self.db, label=self.compared_obj_name,
                                         limit=6, action=self.view_name)
        return super(VennBaseView, self).dispatch(request, *args, **kwargs)

    def get_context(self, **kwargs):
        context = {
            "page_title": self.page_title,
            "comp_type": self.comp_type,
            "compared_obj_name": self.compared_obj_name,
            "form_venn": self.form,
        }
        if getattr(self, "show_results", False):
            context.update({
                "show_results": True,
                "table_headers": self.table_headers,
                "table_data_descr": self.table_data_descr,
                "series": self.series,
                "data_dict": self.data_dict,
                })
        context.update(kwargs)
        return my_locals(context)

    def get(self, request, *args, **kwargs):
        self.form = self.form_class()
        return render(request, self.template, self.get_context())

    @staticmethod
    def _to_percent(count, tot):
        return (count / tot * 100).astype(int)

    def post(self, request, *args, **kwargs):
        self.form = self.form_class(request.POST)
        if not self.form.is_valid():
            self.form = self.form_class()
            # add error message in web page
            return render(request, self.template, self.get_context())

        self.targets = self.form.get_taxids()
        return self.render_venn(request)

    def render_venn(self, request):
        genomes = self.db.get_genomes_description()
        counts = self.get_counts(self.targets)

        self.series = []
        for taxon, taxon_counts in counts.items():
            self.series.append({
                "name": genomes.loc[int(taxon)].description,
                "data": [self.format_entry(key)
                         for key, cnt in taxon_counts.items() if cnt > 0]
                })

        context = self.prepare_data(counts, genomes)
        return render(request, self.template, context)


class VennOrthogroupView(VennBaseView):

    comp_type = "orthogroup"
    table_headers = ["Orthogroup", "Gene", "Description"]
    table_data_descr = "The table contains a list of the genes annotated "\
                       "in each Orthogroup and their description."

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_orthogroup(entry, to_url=to_url)

    @property
    def get_counts(self):
        return self.db.get_og_count

    def prepare_data(self, counts, genomes):
        og_list = counts.index.tolist()
        annotations = self.db.get_genes_from_og(
            orthogroups=og_list, taxon_ids=genomes.index.tolist())
        grouped = annotations.groupby("orthogroup")
        genes = grouped["gene"].apply(list)
        products = grouped["product"].apply(list)

        self.data_dict = {}
        for og in og_list:
            gene_data = "-"
            if og in genes.index:
                g = genes.loc[og]
                gene_data = format_lst_to_html(g, add_count=False)
            prod_data = "-"
            if og in products.index:
                p = products.loc[og]
                prod_data = format_lst_to_html(p, add_count=False)
            self.data_dict[self.format_entry(og)] = [
                self.format_entry(og, to_url=True), gene_data, prod_data]
        self.show_results = True
        return self.get_context()


class VennPfamView(VennBaseView):

    comp_type = "pfam"
    table_headers = ["PFAM", "Description"]
    table_data_descr = "The table contains a list of the Pfam entries and "\
                       "their description."

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_pfam(entry, to_url=to_url)

    @property
    def get_counts(self):
        return self.db.get_pfam_hits

    def prepare_data(self, counts, genomes):
        data = self.db.get_pfam_def(counts.index.tolist())
        self.data_dict = {}
        for pfam, pfam_info in data.iterrows():
            self.data_dict[self.format_entry(pfam)] = [
                self.format_entry(pfam, to_url=True), pfam_info["def"]]
        self.show_results = True
        return self.get_context()


class VennKoView(VennBaseView):

    comp_type = "ko"
    table_headers = ["KO", "Description"]
    table_data_descr = "The table contains a list of the Kegg Orthologs and "\
                       "their description."

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_ko(entry, as_url=to_url)

    def get_counts(self, targets):
        counts = self.db.get_ko_count(targets)["count"].unstack(
            level=0, fill_value=0)
        return counts

    def prepare_data(self, counts, genomes):
        ko_list = counts.index.get_level_values("KO").unique().to_list()
        data = self.db.get_ko_desc(ko_list)
        self.data_dict = {}
        for ko, ko_desc in data.items():
            self.data_dict[self.format_entry(ko)] = [
                self.format_entry(ko, to_url=True), ko_desc]
        self.show_results = True
        return self.get_context()


class VennKoSubsetView(VennKoView):

    def dispatch(self, request, category, *args, **kwargs):
        self.category = category.replace("+", " ")
        self.request = request
        return super(VennKoSubsetView, self).dispatch(request, *args, **kwargs)

    def get_counts(self, targets):
        counts = self.db.get_ko_count_cat(
            taxon_ids=targets, subcategory_name=self.category, index=False)
        counts = counts.drop("module_id", axis=1)
        counts = counts.drop_duplicates(subset=["taxon_id", "KO"])
        counts = counts.set_index(["taxon_id", "KO"])
        counts = counts.unstack(level=0, fill_value=0)["count"]
        return counts

    def get(self, request, *args, **kwargs):
        self.form = self.form_class()
        try:
            self.targets = [int(i) for i in self.request.GET.getlist('h')]
        except Exception:
            return render(request, self.template, self.get_context())

        if len(self.targets) > 5:
            self.targets = self.targets[0:6]

        return self.render_venn(request)


class VennCogView(VennBaseView):

    comp_type = "cog"
    table_headers = ["ID", "Category", "Description"]
    table_data_descr = "The table contains a list of COG definitions, their "\
                       "description and the category to which they belong."

    @staticmethod
    def format_entry(entry, to_url=False):
        return format_cog(entry, as_url=to_url)

    def get_counts(self, targets):
        return self.db.get_cog_hits(targets, search_on="taxid")

    def prepare_data(self, counts, genomes):
        data = self.db.get_cog_summaries(
            counts.index.tolist(), only_cog_desc=True, as_df=True)
        cog_codes = self.db.get_cog_code_description()
        self.data_dict = {}
        for cog, cog_data in data.iterrows():
            functions = [f"{cog_codes[abbr]} ({abbr})"
                         for abbr in cog_data.function]
            functions = format_lst_to_html(functions, False)
            self.data_dict[self.format_entry(cog)] = [
                self.format_entry(cog, to_url=True),
                functions,
                cog_data.description]
        self.show_results = True
        return self.get_context()


def cog_venn_subset(request, category):
    # Note: add error handling

    biodb = settings.BIODB_DB_PATH
    db = DB.load_db(biodb, settings.BIODB_CONF)

    targets = [int(i) for i in request.GET.getlist('h')]
    if len(targets) > 5:
        targets = targets[0:6]

    cog_hits = db.get_cog_hits(targets, search_on="taxid")
    genome_desc = db.get_genomes_description().description.to_dict()
    cog_description = db.get_cog_summaries(cog_hits.index.tolist(), only_cog_desc=True, as_df=True)
    selected_cogs = cog_description[cog_description.function.str.contains(category)]
    cog_codes = db.get_cog_code_description()

    cog2description_l = []
    for cog, data in selected_cogs.iterrows():
        name = data.description
        func = data.function
        functions = ",".join(f"\"{abbr}\"" for abbr in func)
        cog2description_l.append(f"h[\"{format_cog(cog)}\"] = [[{functions}], \"{name}\"]")
    cog2description = ";".join(cog2description_l)

    sel_cog_ids = selected_cogs.index
    cog_hits = cog_hits.reindex(sel_cog_ids)

    series_tab = []
    for target in targets:
        cogs = cog_hits[target]
        non_zero_cogs = cogs[cogs > 0]
        data = ",".join(f"\"{format_cog(cog)}\"" for cog, count in non_zero_cogs.items())
        series_tab.append(f"{{name: \"{genome_desc[target]}\", data: [{data}]}}")
    series = "[" + ",".join(series_tab) + "]"

    cog_func_dict = (f"\"{func}\": \"{descr}\"" for func, descr in cog_codes.items())
    cog_func_dict = "{" + ",".join(cog_func_dict) + "}"
    display_form = False
    envoi_venn = True
    return render(request, 'chlamdb/venn_cogs.html', my_locals(locals()))
