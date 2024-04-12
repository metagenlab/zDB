from chlamdb.forms import CustomPlotsForm
from django.conf import settings
from django.shortcuts import render
from django.views import View
from ete3 import Tree
from lib.db_utils import DB
from lib.ete_phylo import EteTree, SimpleColorColumn

from views.errors import errors
from views.mixins import ComparisonViewMixin
from views.object_type_metadata import my_locals
from views.utils import EntryIdIdentifier, ResultTab, TabularResultTab


class CusomPlotsView(View):

    title = "Custom plots"
    description = "Produce phylogenetic trees including annotations of your choice."
    template = 'chlamdb/custom_plots.html'
    _db = None

    form_class = CustomPlotsForm

    @property
    def view_name(self):
        return "custom_plots"

    def get_result_tabs(self, table):
        return [
            ResultTab("phylogenetic_tree", "Phylogenetic tree",
                      "chlamdb/result_asset.html",
                      asset_path=getattr(self, "tree_path", None)),
            TabularResultTab(
                "gwas_table", "Table",
                table_headers=table["headers"],
                table_data=table["data"],
                table_data_accessors=table["accessors"],
                )
            ]

    @property
    def db(self):
        if self._db is None:
            biodb_path = settings.BIODB_DB_PATH
            self._db = DB.load_db_from_name(biodb_path)
        return self._db

    def get_context(self, **kwargs):
        context = {
            "page_title": self.title,
            "description": self.description,
            "form": self.form
        }
        context.update(kwargs)
        return my_locals(context)

    def get(self, request, *args, **kwargs):
        self.form = self.form_class()
        return render(request, self.template, self.get_context())

    def post(self, request, *args, **kwargs):
        self.form = self.form_class(request.POST)
        if not self.form.is_valid():
            self.form = self.form_class()
            return render(request, self.template,
                          self.get_context(**errors["invalid_form"]))

        entries = self.form.get_entries()

        # We make 1 query for each entry, although we could of course make
        # a single query for each object type, but I don't expect any
        # performance issues here, so I'd rather keep it simple (and maintain
        # the order of the entries as defined by the user).
        entry_id_identifier = EntryIdIdentifier()
        counts = []
        for entry in entries:
            object_type, entry_id = entry_id_identifier.id_to_object_type(entry)
            print(entry, entry_id)
            mixin = ComparisonViewMixin.type2mixin[object_type]
            hits = mixin().get_hit_counts([entry_id], search_on=object_type)
            hits = hits.rename(str).rename({entry_id: entry})
            counts.append(hits)

        genome_descriptions = self.db.get_genomes_description()
        self.prepare_tree(counts, genome_descriptions)
        table = self.prepare_table(counts, genome_descriptions)

        context = self.get_context(
            show_results=True,
            result_tabs=self.get_result_tabs(table)
            )
        return render(request, self.template, context)

    def prepare_tree(self, counts_list, genome_descriptions):
        ref_tree = self.db.get_reference_phylogeny()
        ref_names = genome_descriptions.description.to_dict()

        tree = Tree(ref_tree)
        R = tree.get_midpoint_outgroup()
        if R is not None:
            tree.set_outgroup(R)
        tree.ladderize()
        e_tree = EteTree(tree)
        e_tree.rename_leaves(ref_names)
        for counts in counts_list:
            for label, count in counts.iterrows():
                col = SimpleColorColumn.fromSeries(count, header=label, color_gradient=True)
                e_tree.add_column(col)
        self.tree_path = "/temp/custom_tree.svg"
        path = settings.BASE_DIR + "/assets/" + self.tree_path
        e_tree.render(path, dpi=500)
        return e_tree

    def prepare_table(self, counts_list, genome_descriptions):
        accessors = ["description"]
        headers = ["Name"]
        data = genome_descriptions
        for counts in counts_list:
            data = data.merge(counts.T, how="left", left_on="taxon_id", right_index=True)
            accessors.extend(counts.index)
            headers.extend(counts.index)
        data.where(data.notna(), 0, inplace=True)
        data[accessors[1:]] = data[accessors[1:]].astype(int)
        return {"data": data, "headers": headers, "accessors": accessors}
