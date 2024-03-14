
from chlamdb.forms import make_venn_from
from django.shortcuts import render
from django.views import View

from views.errors import errors
from views.mixins import (AmrViewMixin, CogViewMixin, KoViewMixin,
                          OrthogroupViewMixin, PfamViewMixin, VfViewMixin)


def escape_quotes(unsafe):
    return unsafe.replace("\"", "\\\"")


class VennBaseView(View):

    template = 'chlamdb/venn_generic.html'

    _table_data_descr = "The table contains a list of the {0}, their {1}."

    @property
    def view_name(self):
        return f"venn_{self.object_type}"

    def dispatch(self, request, *args, **kwargs):
        self.form_class = make_venn_from(self.db, label=self.object_name_plural,
                                         limit=6, action=self.view_name)
        return super(VennBaseView, self).dispatch(request, *args, **kwargs)

    def get_context(self, **kwargs):
        context = super(VennBaseView, self).get_context(**kwargs)
        if getattr(self, "show_results", False):
            context.update({
                "show_results": True,
                "table_headers": self.table_headers,
                "table_data_descr": self.table_data_descr,
                "series": self.series,
                "data_dict": self.data_dict,
                })
        return context

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
            return render(request, self.template,
                          self.get_context(**errors["invalid_form"]))

        self.targets = self.form.get_taxids()
        if len(self.targets) > 6:
            self.form = self.form_class()
            return render(request, self.template,
                          self.get_context(**errors["too_many_targets"]))
        return self.render_venn(request)

    def render_venn(self, request):
        genomes = self.db.get_genomes_description()
        counts = self.get_hit_counts(self.targets)
        if counts.empty:
            return render(request, self.template,
                          self.get_context(**errors["no_hits"]))

        counts = self.prepare_data(counts, genomes)

        self.series = []
        for taxon, taxon_counts in counts.items():
            self.series.append({
                "name": genomes.loc[int(taxon)].description,
                "data": [self.format_entry(key)
                         for key, cnt in taxon_counts.items() if cnt > 0]
                })
        context = self.get_context()
        return render(request, self.template, context)

    def filter_data(self, data, counts):
        return data, counts

    def prepare_data(self, counts, genomes):
        data = self.get_hit_descriptions(counts.index.tolist())
        data, counts = self.filter_data(data, counts)
        self.data_dict = {}
        for key, info in data.iterrows():
            row = [info[accessor]
                   for accessor in self.table_data_accessors]
            row = [el if el is not None else "-" for el in row]
            self.data_dict[self.format_entry(key)] = row

        self.show_results = True
        return counts

    @property
    def table_data_descr(self):
        return self._table_data_descr.format(
            self.object_name_plural, self._format_column_headers_to_str())


class VennOrthogroupView(VennBaseView, OrthogroupViewMixin):

    pass


class VennPfamView(VennBaseView, PfamViewMixin):

    pass


class VennKoView(VennBaseView, KoViewMixin):

    pass


class VennSubsetBaseView(VennBaseView):

    def dispatch(self, request, category, *args, **kwargs):
        self.category = category.replace("+", " ")
        self.request = request
        return super(VennSubsetBaseView, self).dispatch(request, *args, **kwargs)

    def get(self, request, *args, **kwargs):
        self.form = self.form_class()
        try:
            self.targets = [int(i) for i in self.request.GET.getlist('h')]
        except Exception:
            return render(request, self.template, self.get_context())

        if len(self.targets) > 5:
            self.targets = self.targets[0:6]

        return self.render_venn(request)


class VennKoSubsetView(VennSubsetBaseView, KoViewMixin):

    table_data_accessors = ["ko", "description"]

    def get_hit_counts(self, targets):
        counts = self.db.get_ko_count_cat(
            taxon_ids=targets, subcategory_name=self.category, index=False)
        counts = counts.drop("module_id", axis=1)
        counts = counts.drop_duplicates(subset=["taxon_id", "KO"])
        counts = counts.set_index(["taxon_id", "KO"])
        counts = counts.unstack(level=0, fill_value=0)["count"]
        return counts


class VennCogView(VennBaseView, CogViewMixin):

    def filter_data(self, data, counts):
        return data, counts.reindex(data.index)


class VennCogSubsetView(VennSubsetBaseView, VennCogView):

    def filter_data(self, data, counts):
        filtered_data = data[data.function.str.contains(self.category)]
        return (filtered_data, counts.reindex(filtered_data.index))


class VennAmrView(VennBaseView, AmrViewMixin):

    pass


class VennVfView(VennBaseView, VfViewMixin):

    pass
