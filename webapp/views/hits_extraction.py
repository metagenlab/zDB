import collections

import pandas as pd
from chlamdb.forms import make_extract_form
from django.conf import settings
from django.shortcuts import render
from django.views import View
from lib.db_utils import DB

from views.mixins import AmrAnnotationsMixin, ComparisonViewMixin
from views.utils import (format_cog, format_cog_url, format_gene_to_ncbi_hmm,
                         format_ko, format_ko_modules, format_ko_path,
                         format_ko_url, format_locus, format_lst_to_html,
                         format_orthogroup, format_pfam, my_locals)

ResultTab = collections.namedtuple("Tab", ["id", "title", "template"])


class ExtractHitsBaseView(View, ComparisonViewMixin):

    template = 'chlamdb/extract_hits.html'

    @property
    def view_name(self):
        return f"extract_{self.comp_type}"

    @property
    def result_tabs(self):
        return [
            ResultTab(1, self.compared_obj_name, "chlamdb/extract_hits_results_table.html"),
            ]

    @property
    def table_count_headers(self):
        return [f"Presence in selection",
                f"Presence in database"]

    @property
    def table_headers(self):
        return self._table_headers + self.table_count_headers

    def dispatch(self, request, *args, **kwargs):
        biodb_path = settings.BIODB_DB_PATH
        self.db = DB.load_db_from_name(biodb_path)
        self.extract_form_class = make_extract_form(
            self.db, self.view_name, plasmid=True, label=self.compared_obj_name)
        return super(ExtractHitsBaseView, self).dispatch(request, *args, **kwargs)

    def get_context(self, **kwargs):
        context = {
            "page_title": self.page_title,
            "comp_type": self.comp_type,
            "compared_obj_name": self.compared_obj_name,
            "table_help": self.table_help,
            "form": self.form,
        }
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
        context.update(kwargs)

        return my_locals(context)

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
        single_copy = "checkbox_single_copy" in request.POST

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
        if not single_copy:
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
            self.selection, search_on=self.comp_type)
        self.n_hits = len(hit_counts_all.index)
        self.max_n = len(hit_counts_all.columns)

        if not single_copy:
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


class ExtractOrthogroupView(ExtractHitsBaseView):

    comp_type = "orthogroup"

    _table_headers = ["Orthogroup", "Genes", "Products"]

    table_help = """
    Two tables have been generated:<br>
    <b>Orthogroups</b> table: it contains the list of orthologous groups shared
    among the selected genomes. The annotation(s) of orthologous groups is a
    consensus of the annotation of all members of the group, and only the two
    most frequent annotations are reported.  For each orthogroup the table
    displays the gene name, the product and the COG category. Additionally,
    the number of occurences of the annotation in the whole set database and
    in the selected genomes.<br>
    <b>Table detail</b>: a complete list of the members of each orthologous
    group shared by the selected genome is displayed. Gene loci are reported
    and quickly linked to additional details about locus annotations.
    """

    @property
    def result_tabs(self):
        return [
            ResultTab(1, self.compared_obj_name, "chlamdb/extract_hits_results_table.html"),
            ResultTab(2, "Details table", "chlamdb/extract_hits_details_table.html"),
            ]

    @property
    def table_headers(self):
        return self._table_headers + self.opt_header + self.table_count_headers

    @property
    def get_hit_counts(self):
        return self.db.get_og_count

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

        self.opt_header, optional_annotations = get_optional_annotations(
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

        details_header, details_data = get_table_details(self.db, annotations)

        self.show_results = True
        context = self.get_context(ref_genomes=ref_genomes,
                                   details_header=details_header,
                                   details_data=details_data)
        return context


class ExtractPfamView(ExtractHitsBaseView):

    comp_type = "pfam"

    _table_headers = ["Pfam entry", "Description"]

    table_help = """
    <br> <b>Pfam</b> table: it contains the list of Pfam domains shared among
    the selected genomes. Pfam entry, its description, its frequency in the
    selected genomes and in all genomes is reported.
    """

    @property
    def get_hit_counts(self):
        return self.db.get_pfam_hits

    def prepare_data(self, hit_counts, hit_counts_all):
        self.table_data = []
        pfam_defs = self.db.get_pfam_def(self.selection)
        for pfam in self.selection:
            pfam_def = pfam_defs["def"].loc[pfam]
            data = [format_pfam(pfam, to_url=True), pfam_def,
                    hit_counts.presence.loc[pfam], hit_counts_all.loc[pfam]]
            self.table_data.append(data)

        self.show_results = True
        return self.get_context()


class ExtractAmrView(ExtractHitsBaseView, AmrAnnotationsMixin):

    comp_type = "amr"

    _table_headers = ["Gene", "Description", "Scope", "Type", "Class",
                      "Subclass"]

    table_help = """
    <br> <b>Pfam</b> table: it contains the list of Pfam domains shared among
    the selected genomes. Pfam entry, its description, its frequency in the
    selected genomes and in all genomes is reported.
    """

    @property
    def get_hit_counts(self):
        return self.db.get_amr_hit_counts

    def prepare_data(self, hit_counts, hit_counts_all):
        self.table_data = []
        # retrieve annotations
        amr_annotations = self.db.get_amr_descriptions(self.selection)
        self.aggregate_amr_annotations(amr_annotations)

        for gene in self.selection:
            amr_annot = amr_annotations[amr_annotations.gene == gene].iloc[0]
            data = [format_gene_to_ncbi_hmm(amr_annot[["gene", "hmm_id"]]),
                    amr_annot.seq_name, amr_annot.scope, amr_annot.type,
                    amr_annot["class"], amr_annot.subclass,
                    hit_counts.presence.loc[gene], hit_counts_all.loc[gene]]
            self.table_data.append(data)

        self.show_results = True
        return self.get_context()


class ExtractKoView(ExtractHitsBaseView):

    comp_type = "ko"

    _table_headers = ["KO", "Description", "Kegg Pathways", "Kegg Modules"]

    table_help = """
    The output consists of a table containing the KO accession number of
    the Kegg Orthologs shared among the selected genomes. For each entry
    its description, the corresponding EC numbers used in Enzyme nomenclature,
    the Kegg Pathways and Kegg modules to whihch it belongs are reported.
    Additionally the frequency of the Ko entry in the selected genomes and in
    the whole database is computed.
    """

    @property
    def get_hit_counts(self):
        return self.db.get_ko_hits

    def prepare_data(self, hit_counts, hit_counts_all):
        ko_desc = self.db.get_ko_desc(self.selection)
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


class ExtractCogView(ExtractHitsBaseView):

    comp_type = "cog"

    _table_headers = ["COG", "Category", "Name"]

    table_help = """
    <br><b>COGs</b> table: it contains the list of COG categories shared among
    the selected genomes. For each category the most frequent annotations and
    the number of occurences in the whole set database and in the selected
    genomes are reported.
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
    def get_hit_counts(self):
        return self.db.get_cog_hits

    @property
    def result_tabs(self):
        return [
            ResultTab(1, self.compared_obj_name, "chlamdb/extract_hits_results_table.html"),
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