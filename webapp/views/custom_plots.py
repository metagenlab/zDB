from chlamdb.forms import make_custom_plots_form
from django.conf import settings
from django.shortcuts import render
from django.views import View
from ete3 import Tree
from lib.db_utils import DB
from lib.ete_phylo import EteTree, SimpleColorColumn

from views.errors import errors
from views.mixins import ComparisonViewMixin
from views.object_type_metadata import my_locals
from views.utils import EntryIdIdentifier, ResultTab


class CusomPlotsView(View):

    title = "Custom plots"
    description = "Produce phylogenetic trees including annotations of your choice."
    template = 'chlamdb/custom_plots.html'
    _db = None

    @property
    def view_name(self):
        return "custom_plots"

    def get_result_tabs(self):
        return [
            ResultTab("phylogenetic_tree", "Phylogenetic tree",
                      "chlamdb/result_asset.html",
                      asset_path=getattr(self, "tree_path", None))
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

    @property
    def form_class(self):
        return make_custom_plots_form(self.db)

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
            mixin = ComparisonViewMixin.type2mixin[object_type]
            hits = mixin().get_hit_counts([entry_id], search_on=object_type)
            counts.append((entry, hits))

        e_tree = self.prepare_tree(counts)
        self.tree_path = "/temp/custom_tree.svg"
        path = settings.BASE_DIR + "/assets/" + self.tree_path
        e_tree.render(path, dpi=500)

        context = self.get_context(
            show_results=True,
            result_tabs=self.get_result_tabs()
            )
        return render(request, self.template, context)

    def prepare_tree(self, counts_list):
        ref_tree = self.db.get_reference_phylogeny()
        ref_names = self.db.get_genomes_description().description.to_dict()

        tree = Tree(ref_tree)
        R = tree.get_midpoint_outgroup()
        if R is not None:
            tree.set_outgroup(R)
        tree.ladderize()
        e_tree = EteTree(tree)
        e_tree.rename_leaves(ref_names)
        for label, counts in counts_list:
            for gene, count in counts.iterrows():
                col = SimpleColorColumn.fromSeries(count, header=label, color_gradient=True)
                e_tree.add_column(col)
        return e_tree
