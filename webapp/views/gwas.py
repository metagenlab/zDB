
import os
from tempfile import TemporaryDirectory

import pandas as pd
from chlamdb.forms import make_venn_from
from django.shortcuts import render
from django.views import View
from ete3 import Tree
from scoary import scoary

from views.analysis_view_metadata import GwasMetadata
from views.mixins import (AmrViewMixin, CogViewMixin, KoViewMixin,
                          OrthogroupViewMixin, PfamViewMixin, VfViewMixin)


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

    def get_context(self, **kwargs):
        context = super(GWASBaseView, self).get_context(**kwargs)
        context["description"] = self.metadata.description
        return context

    def dispatch(self, request, *args, **kwargs):
        self.form_class = make_venn_from(self.db, label=self.object_name_plural,
                                         action=self.view_name)
        self.metadata = GwasMetadata(self.object_type, self.object_name_plural)
        return super(GWASBaseView, self).dispatch(request, *args, **kwargs)

    def get(self, request, *args, **kwargs):
        self.form = self.form_class()
        return render(request, self.template, self.get_context())

    def post(self, request, *args, **kwargs):
        self.form = self.form_class(request.POST)
        if not self.form.is_valid():
            self.form = self.form_class()
            # add error message in web page
            return render(request, self.template, self.get_context())

        self.targets = self.form.get_taxids()
        results = self.run_gwas()
        if results is None:
            context = self.get_context(
                error=True, error_title="No significant association",
                error_message="No association of {} with the phenotype had p-value<0.05")
            return render(request, self.template, context)

        descriptions = self.get_hit_descriptions(results["Gene"].to_list())
        results = results.merge(descriptions,
                                left_on="Gene",
                                right_on=descriptions.index)

        for key in self._gwas_data_accessors:
            results[key] = results[key].apply(self.format_float)

        context = self.get_context(
            table_headers=self.table_headers,
            table_data_accessors=self.table_data_accessors,
            table_data=results,
            show_results=True,
            )
        return render(request, self.template, context)

    def run_gwas(self):
        genomes_data = self.db.get_genomes_infos()
        all_taxids_str = [str(i) for i in genomes_data.index.to_list()]
        all_taxids = [i for i in genomes_data.index.to_list()]
        all_hits = self.get_hit_counts(all_taxids_str, search_on="taxid")

        for i in all_taxids:
            if i not in all_hits:
                all_hits[i] = 0

        with TemporaryDirectory() as tmp:
            genotype_path = os.path.join(tmp, "genotype")
            with open(genotype_path, "w") as fh:
                fh.write(",".join(["Gene"] + all_taxids_str)+"\n")
                for entry, counts in all_hits.iterrows():
                    row = [str(entry), *[str(i) for i in counts[all_taxids]]]
                    fh.write(",".join(row) + "\n")

            phenotype_path = os.path.join(tmp, "phenotype")
            with open(phenotype_path, "w") as fh:
                fh.write("Trait,trait-1\n")
                for taxid in all_taxids:
                    fh.write(f"{taxid},{taxid in self.targets}\n")

            tree_path = os.path.join(tmp, "tree.nw")
            tree = Tree(self.db.get_reference_phylogeny())
            tree.write(format=8, outfile=tree_path)

            scoary(genotype_path, phenotype_path, os.path.join(tmp, "out"),
                   newicktree=tree_path)
            # Summary file is absent if no gene passed filtering
            if not os.path.isfile(os.path.join(tmp, "out", "summary.tsv")):
                return None
            res = pd.read_csv(
                os.path.join(tmp, "out", "traits", "trait-1", "result.tsv"),
                sep="\t")

        return res

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
