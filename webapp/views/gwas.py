import os
from tempfile import TemporaryDirectory

import pandas as pd
from chlamdb.forms import make_gwas_form
from django.conf import settings
from django.shortcuts import render
from django.views import View
from ete3 import Tree
from lib.ete_phylo import EteTree, MatchingColorColumn, ValueColoredColumn
from scoary import scoary

from views.analysis_view_metadata import GwasMetadata
from views.errors import errors
from views.mixins import (AmrViewMixin, CogViewMixin, KoViewMixin,
                          OrthogroupViewMixin, PfamViewMixin, VfViewMixin)
from views.utils import ResultTab, TabularResultTab


class GWASBaseView(View):

    template = 'chlamdb/gwas.html'
    _gwas_data_accessors = ["sensitivity", "specificity", "fisher_p",
                            "fisher_q", "empirical_p", "best", "worst", "fq*ep"]
    _specific_colname_to_header_mapping = {
        "Gene": "Entry",
        "fisher_p": "p-value",
        "fisher_q": "q-value",
        "empirical_p": "Empirical<br>p-value",
        "best": "Best pairwise<br>p-value",
        "worst": "Worst pairwise<br>p-value",
        "fq*ep": "fq*ep",
    }

    @property
    def view_name(self):
        return f"gwas_{self.object_type}"

    @property
    def table_data_accessors(self):
        mixin_accessors = super(GWASBaseView, self).table_data_accessors
        return [mixin_accessors[0],
                *self._gwas_data_accessors,
                *mixin_accessors[1:]]

    def get_result_tabs(self, results):
        return [
            TabularResultTab(
                "gwas_table", "Table", "chlamdb/result_table.html",
                table_headers=self.table_headers,
                table_data=results,
                table_data_accessors=self.table_data_accessors,
                display_index=True,
                colvis_button=True
                ),
            ResultTab("gwas_tree", "Phylogenetic tree",
                      "chlamdb/result_asset.html",
                      asset_path=getattr(self, "tree_path", None))
            ]

    def get_context(self, **kwargs):
        context = super(GWASBaseView, self).get_context(**kwargs)
        context["description"] = self.metadata.description
        return context

    def dispatch(self, request, *args, **kwargs):
        self.form_class = make_gwas_form(self.db)
        self.metadata = GwasMetadata(self.object_type, self.object_name_plural)
        return super(GWASBaseView, self).dispatch(request, *args, **kwargs)

    def get(self, request, *args, **kwargs):
        self.form = self.form_class()
        return render(request, self.template, self.get_context())

    def post(self, request, *args, **kwargs):
        self.form = self.form_class(request.POST, request.FILES)
        if not self.form.is_valid():
            self.form = self.form_class()
            return render(request, self.template,
                          self.get_context(**errors["invalid_form"]))

        qval_threshold = self.form.cleaned_data["bonferroni_cutoff"]
        max_genes = self.form.cleaned_data["max_number_of_hits"]
        phenotype = self.get_phenotype()

        if phenotype is None:
            context = self.get_context(
                error=True, error_title="Invalid phenotype file",
                error_message="File could not be parsed and matched to genomes.")
            return render(request, self.template, context)

        results, hits = self.run_gwas(phenotype, qval_threshold, max_genes)
        if results is None:
            context = self.get_context(
                error=True, error_title="No significant association",
                error_message=f"No association of {self.object_name} with the"
                f" phenotype had q-value<{qval_threshold}")
            return render(request, self.template, context)

        descriptions = self.get_hit_descriptions(results["Gene"].to_list())
        results = results.merge(descriptions,
                                left_on="Gene",
                                right_on=descriptions.index)

        for key in self._gwas_data_accessors:
            results[key] = results[key].apply(self.format_float)

        e_tree = self.prepare_tree(phenotype, hits, results)
        self.tree_path = "/temp/gwas_tree.svg"
        path = settings.BASE_DIR + "/assets/" + self.tree_path
        e_tree.render(path, dpi=500)

        context = self.get_context(
            show_results=True,
            qval_threshold=qval_threshold,
            result_tabs=self.get_result_tabs(results)
            )
        return render(request, self.template, context)

    def prepare_tree(self, phenotype, hits, results):
        # unique_og = intersect.orthogroup.unique().tolist()
        # red_color = set(tuple(entry) for entry in intersect.to_numpy())
        # df_og_count = db.get_og_count(list(unique_og), search_on="orthogroup").T
        ref_tree = self.db.get_reference_phylogeny()
        ref_names = self.db.get_genomes_description().description.to_dict()

        tree = Tree(ref_tree)
        R = tree.get_midpoint_outgroup()
        if R is not None:
            tree.set_outgroup(R)
        tree.ladderize()
        e_tree = EteTree(tree)
        e_tree.rename_leaves(ref_names)

        phenotype.trait = phenotype.trait.astype(int)
        phenotype = phenotype.set_index("taxids")
        neg_phenotype = 1 - phenotype
        neg_phenotype = neg_phenotype.to_dict()["trait"]
        phenotype = phenotype.to_dict()["trait"]

        e_tree.add_column(ValueColoredColumn(
            phenotype,
            header="Phenotype",
            col_func=lambda x: EteTree.BLUE if x == 0 else EteTree.RED))

        positive = results["supporting"] > results["opposing"]
        for ((gene, hit), pos) in zip(hits.iterrows(), positive):
            col_column = MatchingColorColumn(
                hit.to_dict(), phenotype if pos else neg_phenotype, header=gene,
                col_func=lambda x: EteTree.BLUE if x == 0 else EteTree.RED)
            e_tree.add_column(col_column)
        return e_tree

    def get_phenotype(self):
        phenotype = pd.read_csv(self.form.cleaned_data["phenotype_file"],
                                header=None, names=["taxids", "trait"])
        phenotype["trait"] = phenotype["trait"].astype(bool)
        genomes = self.db.get_genomes_description()
        if all(phenotype.taxids.isin(genomes.index)):
            return phenotype
        elif all(phenotype.taxids.isin(genomes.description)):
            mapping = {genome.description: taxid
                       for taxid, genome in genomes.iterrows()}
            phenotype.taxids = phenotype.taxids.apply(lambda x: mapping[x])
        else:
            return None
        return phenotype

    def run_gwas(self, phenotype, qval_threshold, max_genes):
        genomes_data = self.db.get_genomes_infos()
        all_taxids_str = [str(i) for i in genomes_data.index.to_list()]
        all_taxids = [i for i in genomes_data.index.to_list()]
        all_hits = self.get_hit_counts(all_taxids, search_on="taxid")
        if all_hits.empty:
            return (None, None)

        with TemporaryDirectory() as tmp:
            genotype_path = os.path.join(tmp, "genotype")
            with open(genotype_path, "w") as fh:
                fh.write(",".join(["Gene"] + all_taxids_str)+"\n")
                for entry, counts in all_hits.iterrows():
                    row = [str(entry), *[str(i) for i in counts[all_taxids]]]
                    fh.write(",".join(row) + "\n")

            phenotype_path = os.path.join(tmp, "phenotype")
            phenotype.to_csv(phenotype_path, index=False)

            tree_path = os.path.join(tmp, "tree.nw")
            tree = Tree(self.db.get_reference_phylogeny())
            tree.write(format=8, outfile=tree_path)

            scoary(genotype_path,
                   phenotype_path,
                   os.path.join(tmp, "out"),
                   newicktree=tree_path,
                   multiple_testing=f'bonferroni:{qval_threshold}',
                   max_genes=max_genes)

            # Summary file is absent if no gene passed filtering
            if not os.path.isfile(os.path.join(tmp, "out", "summary.tsv")):
                return (None, None)
            res = pd.read_csv(
                os.path.join(tmp, "out", "traits", "trait", "result.tsv"),
                sep="\t")

        return res, all_hits.loc[res["Gene"]]

    @staticmethod
    def format_float(val, precision=2):
        return f"{val:.{precision}g}"


class AmrGwasView(GWASBaseView, AmrViewMixin):

    pass


class PfamGwasView(GWASBaseView, PfamViewMixin):

    pass


class KoGwasView(GWASBaseView, KoViewMixin):

    pass


class CogGwasView(GWASBaseView, CogViewMixin):

    pass


class OrthogroupGwasView(GWASBaseView, OrthogroupViewMixin):

    pass


class VfGwasView(GWASBaseView, VfViewMixin):

    pass
