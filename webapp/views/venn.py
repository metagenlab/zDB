
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
                "envoi_venn": True,
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

        targets = self.form.get_taxids()
        genomes = self.db.get_genomes_description()
        counts = self.db.get_og_count(targets)
        context = self.prepare_data(counts, genomes)
        return render(request, self.template, context)


class VennOrthogroupView(VennBaseView):

    template = 'chlamdb/venn_orthogroup.html'
    comp_type = "orthogroup"

    @property
    def get_counts(self):
        return self.db.get_og_count

    def prepare_data(self, counts, genomes):
        fmt_data = []
        for taxon in counts:
            ogs = counts[taxon]
            ogs_str = ",".join(f"{to_s(format_orthogroup(og))}"
                               for og, cnt in ogs.items() if cnt > 0)
            genome = genomes.loc[int(taxon)].description
            fmt_data.append(f"{{name: {to_s(genome)}, data: [{ogs_str}]}}")
        series = "[" + ",".join(fmt_data) + "]"

        og_list = counts.index.tolist()
        annotations = self.db.get_genes_from_og(
            orthogroups=og_list, taxon_ids=genomes.index.tolist())
        grouped = annotations.groupby("orthogroup")
        genes = grouped["gene"].apply(list)
        products = grouped["product"].apply(list)

        orthogroup2description = []
        for og in og_list:
            forbidden = "\""
            gene_data = "-"
            if og in genes.index:
                g = genes.loc[og]
                gene_data = format_lst_to_html(g, add_count=False)
            prod_data = "-"
            if og in products.index:
                p = products.loc[og]
                prod_data = format_lst_to_html(p, add_count=False)
            og_info = "[\"" + gene_data + "\",\"" + prod_data + "\"]"
            og_item = f"h[{to_s(format_orthogroup(og))}] = {og_info}"
            orthogroup2description.append(og_item)
        orthogroup2description = "\n".join(orthogroup2description)
        self.show_results = True
        return self.get_context(series=series,
                                orthogroup2description=orthogroup2description)


def venn_pfam(request):
    biodb = settings.BIODB_DB_PATH
    db = DB.load_db(biodb, settings.BIODB_CONF)
    page_title = page2title["venn_pfam"]

    venn_form_class = make_venn_from(db, label="PFAM domain", limit=6, action="venn_pfam")
    if request.method != "POST":
        form_venn = venn_form_class()
        return render(request, 'chlamdb/venn_Pfam.html', my_locals(locals()))

    form_venn = venn_form_class(request.POST)
    if not form_venn.is_valid():
        # add error message
        form_venn = venn_form_class()
        return render(request, 'chlamdb/venn_Pfam.html', my_locals(locals()))

    targets = form_venn.get_taxids()
    pfam_hits = db.get_pfam_hits(targets, search_on="taxid", indexing="taxid")
    data = db.get_pfam_def(pfam_hits.index.tolist())
    genomes_desc = db.get_genomes_description().description.to_dict()

    series_tab = []
    for target in targets:
        pfams = pfam_hits[target]
        non_zero = pfams[pfams > 0]
        str_fmt = ",".join(f"\"{format_pfam(pfam)}\"" for pfam, _ in non_zero.items())
        series_tab.append(f"{{name: \"{genomes_desc[target]}\", data: [{str_fmt}]}}")
    series = "[" + ",".join(series_tab) + "]"

    descriptions = []
    for pfam, pfam_info in data.iterrows():
        pfam_def = escape_quotes(pfam_info["def"])
        descriptions.append(f"h[\"{format_pfam(pfam)}\"] = \"{pfam_def}\"")

    ctx = {"envoi_venn": True,
           "series": series,
           "pfam2description": ";".join(descriptions),
           "form_venn": form_venn}
    return render(request, 'chlamdb/venn_Pfam.html', my_locals(ctx))


def venn_ko(request):
    biodb = settings.BIODB_DB_PATH
    db = DB.load_db(biodb, settings.BIODB_CONF)
    page_title = page2title["venn_ko"]

    venn_form_class = make_venn_from(db, limit=6)
    display_form = True
    if request.method != "POST":
        form_venn = venn_form_class()
        return render(request, 'chlamdb/venn_ko.html', my_locals(locals()))

    form_venn = venn_form_class(request.POST)
    if not form_venn.is_valid():
        # add error message
        form_venn = venn_form_class()
        return render(request, 'chlamdb/venn_ko.html', my_locals(locals()))

    taxids = form_venn.get_taxids()
    genomes = db.get_genomes_description().description.to_dict()
    ko_counts = db.get_ko_count(taxids)

    fmt_data = []
    ko_list = ko_counts.index.get_level_values("KO").unique().to_list()
    for taxid in taxids:
        kos = ko_counts.loc[taxid].index.values
        kos_str = ",".join(f"{to_s(format_ko(ko))}" for ko in kos)
        genome = genomes[taxid]
        fmt_data.append(f"{{name: {to_s(genome)}, data: [{kos_str}]}}")
    series = "[" + ",".join(fmt_data) + "]"

    ko_descriptions = db.get_ko_desc(ko_list)
    ko2description = []
    for ko, ko_desc in ko_descriptions.items():
        cleaned_desc = escape_quotes(ko_desc)
        ko_item = f"h[{to_s(format_ko(ko))}] = [\"{cleaned_desc}\"];"
        ko2description.append(ko_item)
    ko2description = "\n".join(ko2description)
    envoi_venn = True
    return render(request, 'chlamdb/venn_ko.html', my_locals(locals()))


def venn_cog(request, sep_plasmids=False):
    """
    Will need to modify the signature of the method to remove the sep_plasmid
    parameter as it is not taken into account. Or put back the differentiate
    plasmid parameter in the web page.
    """

    biodb = settings.BIODB_DB_PATH
    db = DB.load_db(biodb, settings.BIODB_CONF)
    page_title = page2title["venn_cog"]

    display_form = True
    venn_form_class = make_venn_from(db, limit=6, label="COG")
    if request.method != "POST":
        form_venn = venn_form_class()
        return render(request, 'chlamdb/venn_cogs.html', my_locals(locals()))

    form_venn = venn_form_class(request.POST)
    if not form_venn.is_valid():
        # TODO: add error message
        return render(request, 'chlamdb/venn_cogs.html', my_locals(locals()))

    targets = form_venn.get_taxids()
    cog_hits = db.get_cog_hits(targets, search_on="taxid")
    data = db.get_cog_summaries(cog_hits.index.tolist(), only_cog_desc=True, as_df=True)
    genome_desc = db.get_genomes_description().description.to_dict()

    # necessary as some COG do not have a description
    # --> filter them out
    cog_hits = cog_hits.reindex(data.index)

    series_tab = []
    for target in targets:
        cogs = cog_hits[target]
        non_zero_cogs = cogs[cogs > 0]
        str_fmt = ",".join(f"\"{format_cog(cog)}\"" for cog, count in non_zero_cogs.items())
        series_tab.append(f"{{name: \"{genome_desc[target]}\", data: [{str_fmt}]}}")
    series = "[" + ",".join(series_tab) + "]"

    cog2description_l = []
    cog_codes = db.get_cog_code_description()
    for cog, data in data.iterrows():
        name = escape_quotes(data.description)
        func = data.function
        functions = ",".join(f"\"{abbr}\"" for abbr in func)
        cog2description_l.append(f"h[\"{format_cog(cog)}\"] = [[{functions}], \"{name}\"]")

    cog_func_dict = (f"\"{func}\": \"{descr}\"" for func, descr in cog_codes.items())
    cog_func_dict = "{" + ",".join(cog_func_dict) + "}"
    cog2description = ";".join(cog2description_l)
    envoi_venn = True
    return render(request, 'chlamdb/venn_cogs.html', my_locals(locals()))


def ko_venn_subset(request, category):
    biodb = settings.BIODB_DB_PATH
    db = DB.load_db(biodb, settings.BIODB_CONF)

    category = category.replace("+", " ")
    try:
        targets = [int(i) for i in request.GET.getlist('h')]
    except Exception:
        # add an error message
        return render(request, 'chlamdb/venn_ko.html', my_locals(locals()))

    if len(targets) > 5:
        targets = targets[0:6]

    genomes = db.get_genomes_description(targets).description.to_dict()
    ko_counts = db.get_ko_count_cat(taxon_ids=targets, subcategory_name=category, index=False)
    ko_count = ko_counts.groupby("taxon_id")["KO"].apply(list)

    # shameful copy/cape from venn_ko
    fmt_data = []
    ko_set = set()
    for taxid in targets:
        if taxid not in ko_count.index:
            continue
        kos = ko_count.loc[taxid]
        kos_str = ",".join(f"{to_s(format_ko(ko))}" for ko in kos)
        genome = genomes[taxid]
        fmt_data.append(f"{{name: {to_s(genome)}, data: [{kos_str}]}}")
        ko_set = ko_set.union({ko for ko in kos})
    series = "[" + ",".join(fmt_data) + "]"

    ko_list = list(ko_set)
    ko_descriptions = db.get_ko_desc(ko_list)
    ko2description = []
    for ko, ko_desc in ko_descriptions.items():
        forbidden = "\""
        safe_desc = escape_quotes(ko_desc)
        ko_item = f"h[{to_s(format_ko(ko))}] = {forbidden}{safe_desc}{forbidden}"
        ko2description.append(ko_item)

    ko2description = "".join(ko2description)
    display_form = False
    envoi_venn = True
    return render(request, 'chlamdb/venn_ko.html', my_locals(locals()))
