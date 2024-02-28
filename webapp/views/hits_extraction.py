import collections

import pandas as pd
from chlamdb.forms import make_extract_form
from django.shortcuts import render
from django.views import View

from views.mixins import (AmrViewMixin, CogViewMixin, KoViewMixin,
                          OrthogroupViewMixin, PfamViewMixin, VfViewMixin)
from views.utils import (format_cog, format_cog_url, format_ko,
                         format_ko_modules, format_ko_path, format_ko_url,
                         format_locus, format_lst_to_html, format_orthogroup)

ResultTab = collections.namedtuple("Tab", ["id", "title", "template"])


class ExtractHitsBaseView(View):

    template = 'chlamdb/extract_hits.html'

    results_table_help = """
            <br> <b>{0} table</b>: it contains the list of
            {0} shared among the selected (and absent
            from the excluded) genomes. For each entry, the table lists its
            {1} as well as the fraction of included genomes
            containing it {2} and the fraction of genomes in the whole database
            containing it {2}.
            {3}
            """

    _table_help_complement = ""
    _col_descriptions = {}

    @property
    def table_cols_description(self):
        return [self._col_descriptions.get(header, header.lower())
                for header in self.table_headers[:-2]]

    @property
    def table_help(self):
        once = "exactly once" if getattr(self, "single_copy", False)\
                else "at least once"
        col_descr = ", ".join(self.table_cols_description[:-1])
        col_descr += " and " + self.table_cols_description[-1]
        return self.results_table_help.format(self.object_name_plural,
                                              col_descr,
                                              once,
                                              self._table_help_complement)

    @property
    def view_name(self):
        return f"extract_{self.object_type}"

    @property
    def result_tabs(self):
        return [
            ResultTab(1, self.object_name_plural, "chlamdb/extract_hits_results_table.html"),
            ]

    @property
    def table_count_headers(self):
        return ["Presence in selection (%)",
                "Presence in database (%)"]

    @property
    def _table_headers(self):
        return super(ExtractHitsBaseView, self).table_headers

    @property
    def table_headers(self):
        return self._table_headers + self.table_count_headers

    def dispatch(self, request, *args, **kwargs):
        self.extract_form_class = make_extract_form(
            self.db, self.view_name, plasmid=True, label=self.object_name_plural)
        return super(ExtractHitsBaseView, self).dispatch(request, *args, **kwargs)

    def get_context(self, **kwargs):
        context = super(ExtractHitsBaseView, self).get_context(**kwargs)
        context["table_help"] = self.table_help,
        if getattr(self, "show_results", False):
            context.update({
                "show_results": True,
                "n_missing": self.n_missing,
                "n_hits": self.n_hits,
                "table_headers": self.table_headers,
                "table_data": self.table_data,
                "included_taxids": self.included_taxids,
                "excluded_taxids": self.excluded_taxids,
                "selection": self.selection,
                "result_tabs": self.result_tabs,
                })
        return context

    def get(self, request, *args, **kwargs):
        self.form = self.extract_form_class()
        return render(request, self.template, self.get_context())

    @staticmethod
    def _to_percent(count, tot):
        return (count / tot * 100).astype(int)

    def post(self, request, *args, **kwargs):
        self.form = self.extract_form_class(request.POST)
        if not self.form.is_valid():
            self.form = self.extract_form_class()
            # add error message in web page
            return render(request, self.template, self.get_context())

        # Extract form data
        self.included_taxids, self.included_plasmids = self.form.get_include_choices()
        self.excluded_taxids, self.excluded_plasmids = self.form.get_exclude_choices()
        self.n_missing = self.form.get_n_missing()
        self.single_copy = "checkbox_single_copy" in request.POST

        self.n_included = len(self.included_taxids)
        if self.included_plasmids is not None:
            self.n_included += len(self.included_plasmids)

        if self.n_missing >= self.n_included:
            context = self.get_context(wrong_n_missing=True)
            return render(request, self.template, context)

        self.min_fraq = int((self.n_included - self.n_missing) /
                            self.n_included * 100)

        self.n_excluded = len(self.excluded_taxids)
        if self.excluded_plasmids is not None:
            self.n_excluded += len(self.excluded_plasmids)

        # Count hits for selection
        hit_counts = self.get_hit_counts(self.included_taxids,
                                         plasmids=self.included_plasmids,
                                         search_on="taxid")
        if not self.single_copy:
            hit_counts["presence"] = self._to_percent(
                hit_counts[hit_counts > 0].count(axis=1), self.n_included)
            hit_counts["selection"] = hit_counts.presence >= self.min_fraq
        else:
            excluded_hits = hit_counts[hit_counts > 1].count(axis=1)
            hit_counts["presence"] = self._to_percent(
                hit_counts[hit_counts == 1].count(axis=1), self.n_included)
            hit_counts["selection"] = ((hit_counts.presence >= self.min_fraq)
                                       & (excluded_hits == 0))

        if self.n_excluded > 0:
            mat_exclude = self.get_hit_counts(
                self.excluded_taxids, plasmids=self.excluded_plasmids,
                search_on="taxid")
            mat_exclude["presence"] = mat_exclude[mat_exclude > 0].count(axis=1)
            mat_exclude["exclude"] = mat_exclude.presence > 0
            neg_index = mat_exclude[mat_exclude.exclude].index
        else:
            neg_index = pd.Index([])

        pos_index = hit_counts[hit_counts.selection].index
        self.selection = pos_index.difference(neg_index).tolist()
        if len(self.selection) == 0:
            context = self.get_context(no_match=True)
            return render(request, self.template, context)

        # Count hits for all genomes
        hit_counts_all = self.get_hit_counts(
            self.selection, search_on=self.object_column)
        self.n_hits = len(hit_counts_all.index)
        self.max_n = len(self.db.taxon_ids())

        if not self.single_copy:
            hit_counts_all = self._to_percent(
                hit_counts_all[hit_counts_all > 0].count(axis=1), self.max_n)
        else:
            hit_counts_all = self._to_percent(
                hit_counts_all[hit_counts_all == 1].count(axis=1), self.max_n)
        context = self.prepare_data(hit_counts, hit_counts_all)

        if context is None:
            context = self.get_context(no_match=True)
            return render(request, self.template, context)

        return render(request, self.template, context)

    def prepare_data(self, hit_counts, hit_counts_all):
        self.table_data = []
        # retrieve descriptions
        descriptions = self.get_hit_descriptions(self.selection)
        for entry in self.selection:
            amr_annot = descriptions.loc[entry]
            data = [amr_annot[key] for key in self.table_data_accessors]
            data.extend([hit_counts.presence.loc[entry], hit_counts_all.loc[entry]])
            data = [el if el is not None else "-" for el in data]
            self.table_data.append(data)

        self.show_results = True
        return self.get_context()


class ExtractOrthogroupView(ExtractHitsBaseView, OrthogroupViewMixin):

    _table_headers = ["Orthogroup", "Genes", "Products"]

    _table_help_complement = """
    The annotation(s) of orthologous groups is a
    consensus of the annotation of all members of the group, and only the two
    most frequent annotations are reported.
    <br>
    <b>Details tabke</b>: a complete list of the members of each orthologous
    group shared by the selected genome is displayed. Gene loci are reported
    and quickly linked to additional details about locus annotations.
    """

    _col_descriptions = {"COG": "COG category",
                         "KO": "KO assignment"}

    @property
    def result_tabs(self):
        return [
            ResultTab(1, self.object_name_plural, "chlamdb/extract_hits_results_table.html"),
            ResultTab(2, "Details table", "chlamdb/extract_hits_details_table.html"),
            ]

    @property
    def table_headers(self):
        return self._table_headers + getattr(self, "opt_header", [])\
               + self.table_count_headers

    @staticmethod
    def get_optional_annotations(db, seqids):
        header = []
        config_table = db.get_config_table()
        annotations = []
        if config_table.get("KEGG", False):
            header.append("KO")
            ko_hits = db.get_ko_hits(seqids, search_on="seqid", indexing="seqid")
            annotations.append(ko_hits)
        if config_table.get("COG", False):
            header.append("COG")
            cog_hits = db.get_cog_hits(seqids, indexing="seqid", search_on="seqid")
            annotations.append(cog_hits)

        if len(annotations) == 2:
            return header, annotations[0].join(annotations[1], how="outer")
        elif len(annotations) == 1:
            return header, annotations[0]
        return header, pd.DataFrame()

    @staticmethod
    def get_table_details(db, annotations):
        header = ["Orthogroup", "Organism", "Locus", "Gene", "Product"]
        hsh_organisms = db.get_organism(annotations.index.tolist())
        infos = []
        for seqid, data in annotations.iterrows():
            organism = hsh_organisms[seqid]
            og = format_orthogroup(data.orthogroup, to_url=True)
            gene = data.gene
            if pd.isna(data.gene):
                gene = "-"
            locus = format_locus(data.locus_tag, to_url=True)
            entry = [og, organism, locus, gene, data["product"]]
            infos.append(entry)
        return header, infos

    def prepare_data(self, hit_counts, hit_counts_all):
        self.table_data = []

        all_taxids = self.included_taxids
        if self.included_plasmids is not None:
            all_taxids += self.included_plasmids

        annotations = self.db.get_genes_from_og(
            orthogroups=self.selection,
            taxon_ids=all_taxids,
            terms=["gene", "product", "locus_tag"])
        if annotations.empty:
            return None

        self.opt_header, optional_annotations = self.get_optional_annotations(
            self.db, seqids=annotations.index.tolist())
        annotations = annotations.join(optional_annotations)
        grouped = annotations.groupby("orthogroup")
        genes = grouped["gene"].apply(list)
        products = grouped["product"].apply(list)

        if "COG" in self.opt_header:
            cogs = grouped["cog"].apply(list)

        if "KO" in self.opt_header:
            kos = grouped["ko"].apply(list)

        for row, count in hit_counts_all.items():
            cnt_in = hit_counts.presence.loc[row]
            g = genes.loc[row]
            gene_data = format_lst_to_html(
                "-" if pd.isna(entry) else entry for entry in g)
            prod_data = format_lst_to_html(products.loc[row])
            column_header = format_orthogroup(row, to_url=True)
            optional = []
            if "KO" in self.opt_header and row in kos:
                optional.append(format_lst_to_html(
                    kos.loc[row], add_count=True, format_func=format_ko_url))
            if "COG" in self.opt_header:
                optional.append(format_lst_to_html(
                    cogs.loc[row], add_count=True, format_func=format_cog_url))
            entry = [column_header, gene_data, prod_data, *optional, cnt_in, count]
            self.table_data.append(entry)

        ref_genomes = self.db.get_genomes_description(
        ).loc[self.included_taxids].reset_index()

        details_header, details_data = self.get_table_details(self.db, annotations)

        self.show_results = True
        context = self.get_context(ref_genomes=ref_genomes,
                                   details_header=details_header,
                                   details_data=details_data,
                                   show_circos_form=True)
        return context


class ExtractPfamView(ExtractHitsBaseView, PfamViewMixin):

    table_data_accessors = ["pfam", "def"]


class ExtractAmrView(ExtractHitsBaseView, AmrViewMixin):

    table_data_accessors = ["gene", "seq_name", "scope", "type", "class",
                            "subclass", "hmm_id"]


class ExtractVfView(ExtractHitsBaseView, VfViewMixin):

    table_data_accessors = ["vf_gene_id", "prot_name", "vfid", "category"]


class ExtractKoView(ExtractHitsBaseView, KoViewMixin):

    _table_headers = ["KO", "Description", "Kegg Pathways", "Kegg Modules"]

    _col_descriptions = {"Description": "description including the corresponding"
                                        " EC numbers used in enzyme nomenclature",
                         "Kegg Modules": "Kegg modules to whihch it belongs"}

    def prepare_data(self, hit_counts, hit_counts_all):
        ko_desc = self.get_hit_descriptions(self.selection)
        ko_mod = self.db.get_ko_modules(self.selection)
        ko_path = self.db.get_ko_pathways(self.selection)
        self.table_data = []
        for ko in self.selection:
            kof = format_ko(ko, as_url=True)
            kod = ko_desc.get(ko, "-")
            kop = format_ko_path(ko_path, ko)
            kom = format_ko_modules(ko_mod, ko)
            kot = hit_counts_all.loc[ko]
            data = [kof, kod, kop, kom, hit_counts.presence.loc[ko], kot]
            self.table_data.append(data)

        self.show_results = True
        return self.get_context()


class ExtractCogView(ExtractHitsBaseView, CogViewMixin):

    _table_headers = ["COG", "Function", "Description"]

    _table_help_complement = """
    <br><b> COG categories barchart</b>: this plot displays for each COG
    category in
    <span style="color: rgb(135, 186, 245)"><b>light-blue</b></span> the number
    of genes shared by the selected genomes, while in
    <span style="color: rgb(16, 76, 145)"><b>blue</b></span> , the total number
    of genes annotated with that COG in the selected genomes
    (shared and unique to each selected genome). Next to the light-blue
    barcharts there is a percentage value as the result of the number of
    reported COGs for each category divided by the number of shared COG in all
    categories, while next to the blue barcharts the value represents COGs
    count are divided by the total COGs unique and shared in the selected
    genomes.
    <br>A longer light-blue bar can be interpreted as an enrichment of that
    shared COG category in the selected genomes compared to their complete COG
    profiles.
    <br> <b>Locus list reference</b>
    """

    @property
    def result_tabs(self):
        return [
            ResultTab(1, self.object_name_plural, "chlamdb/extract_hits_results_table.html"),
            ResultTab(2, "COG categories barchart", "chlamdb/extract_hits_cog_barcharts.html"),
            ]

    def prepare_data(self, hit_counts, hit_counts_all):
        cat_count = {}
        cogs_summaries = self.db.get_cog_summaries(hit_counts_all.index.tolist())
        cogs_funct = self.db.get_cog_code_description()
        self.table_data = []
        for cog_id in self.selection:

            # some cogs do not have a description, skip those
            if cog_id not in cogs_summaries:
                continue

            data = [format_cog(cog_id)]
            func_acc = []
            for func, func_descr, cog_descr in cogs_summaries[cog_id]:
                func_acc.append((func, func_descr))
                inc, not_incl = cat_count.get(func, (0, 0))
                cat_count[func] = (inc + hit_counts.presence.loc[cog_id], not_incl)
            funcs = "<br>".join(f"{func} ({func_desc})" for func, func_desc in func_acc)
            data = (format_cog(cog_id, as_url=True), funcs, cog_descr,
                    hit_counts.presence.loc[cog_id], hit_counts_all.loc[cog_id])
            self.table_data.append(data)

        # get the categories for all cogs
        for cog_id, details_lst in cogs_summaries.items():
            for func, func_descr, cog_descr in details_lst:
                inc, not_incl = cat_count.get(func, (0, 0))
                cat_count[func] = (inc, not_incl + hit_counts_all.loc[cog_id])

        # Code to generate the barchart diagrams
        cat_map_str = ",".join([f"\"{func}\": \"{descr}\"" for func, descr in cogs_funct.items()])
        category_map = f"var category_description = {{{cat_map_str}}};"
        ttl_sel = sum([c1 for func, (c1, c2) in cat_count.items()])
        ttl_all = sum([c2 for func, (c1, c2) in cat_count.items()])

        cat_count_comp = ",".join([f"\"{func}\": [\"{c2}\", \"{c1}\"]" for func, (c1, c2) in cat_count.items()])
        category_count_complete = f"var category_count_complete = {{{cat_count_comp}}};"

        serie_selection_val = [str(round(float(c1) / ttl_sel, 2)) for func, (c1, c2) in cat_count.items()]
        serie_all_val = [str(round(float(c2) / ttl_all, 2)) for func, (c1, c2) in cat_count.items()]
        serie_selection = f"{{labels: \"selection\", values: [{','.join(serie_selection_val)}]}}"
        serie_all = f"{{labels: \"complete genomes\", values: [{','.join(serie_all_val)}]}}"
        series_str = ",".join([serie_all, serie_selection])
        series = f"[{series_str}]"
        labels_str = ",".join([f"\"{funct}\"" for funct in cat_count.keys()])
        labels = f"[{labels_str}]"

        self.show_results = True
        return self.get_context(
            labels=labels,
            series=series,
            category_map=category_map,
            category_count_complete=category_count_complete,
            )
