# -*- coding: utf-8 -*-

# todo circos gc file curently written in home directory, move it to other place
# todo save temp files in temp folder

import collections
import random
import re
import string
import time
from io import StringIO
from tempfile import NamedTemporaryFile

import bibtexparser
import chlamdb.circosjs as circosjs
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import (NcbiblastnCommandline,
                                    NcbiblastpCommandline,
                                    NcbiblastxCommandline,
                                    NcbitblastnCommandline)
from chlamdb.forms import (make_blast_form, make_circos_form, make_metabo_from,
                           make_module_overview_form,
                           make_pathway_overview_form, make_plot_form,
                           make_single_genome_form, make_venn_from)
from django.conf import settings
from django.http import JsonResponse
from django.shortcuts import render
from django.views import View
from ete3 import StackedBarFace, TextFace, Tree
from lib import search_bar as sb
from lib.db_utils import DB
from lib.ete_phylo import (Column, EteTree, KOAndCompleteness,
                           ModuleCompletenessColumn, SimpleColorColumn)
from lib.KO_module import ModuleParser
from reportlab.lib import colors

from views.analysis_view_metadata import (AccumulationRarefactionMetadata,
                                          CategoriesBarchartMetadata,
                                          CategoriesCountHeatmapMetadata,
                                          CategoriesFreqHeatmapMetadata,
                                          HeatmapMetadata)
from views.errors import errors
from views.mixins import (BaseViewMixin, CogViewMixin, ComparisonViewMixin,
                          GenomesTableMixin, KoViewMixin)
from views.object_type_metadata import MetadataGetter, my_locals
from views.utils import (TabularResultTab, format_cog, format_gene, format_ko,
                         format_locus, format_orthogroup,
                         genomic_region_df_to_js, get_genomes_data,
                         locusx_genomic_region, make_div, optional2status,
                         page2title, to_s)


def id_generator(size=6, chars=string.ascii_uppercase + string.ascii_lowercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def help(request):
    return render(request, 'chlamdb/help.html', my_locals(locals()))


def about(request):
    path = settings.ASSET_ROOT + '/bibliography/references.bib'
    with open(path) as bibtex_file:
        bib_database = bibtexparser.load(bibtex_file)

    entry_list = []

    for entry in bib_database.entries:
        ref = f"<b>{re.sub('[{}]', '', entry['title'])}</b><br>"\
              f"{entry['author']}, {entry['journal']},"
        if entry.get("volume", False):
            ref += f" {entry['volume']}({entry['number']}):{entry['pages']},"
        ref += f" {entry['year']}"
        url = entry["url"]
        entry_list.append([ref, url])

    return render(request, 'chlamdb/credits.html', my_locals(locals()))


class StackedBarColumn(Column):
    def __init__(self, values, header, colours=None, relative=False,
                 face_params=None, header_params=None):
        super().__init__(header, face_params, header_params)
        self.values = values
        self.colours = colours
        self.relative = relative
        if relative:
            self.max = max(values.values())
            self.min = min(values.values())

    def get_face(self, index):
        val = self.values[int(index)]

        if self.relative and self.max != self.min:
            val = 100 * float(val - self.min) / (self.max - self.min)
        elif self.relative:
            val = 100

        face = StackedBarFace([val, 100 - val], width=50,
                              height=9, colors=self.colours)
        self.set_default_params(face)
        face.inner_border.color = "black"
        face.inner_border.width = 0
        return face


def home(request):
    biodb_path = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb_path)

    hsh_files = db.get_filenames_to_taxon_id()
    number_of_files = len(hsh_files)

    number_ort = db.get_n_orthogroups()

    if optional2status.get("BLAST_swissprot", False):
        number_swissprot = db.get_number_of_swissprot_entries()
    if optional2status.get("vf", False):
        number_vf = db.vf.get_number_of_entries()
    if optional2status.get("amr", False):
        number_amr = db.get_number_of_amr_entries()
    if optional2status.get("cog", False):
        number_cog = db.get_number_of_cog_entries()
    if optional2status.get("ko", False):
        number_ko = db.get_number_of_ko_entries()
    if optional2status.get("pfam", False):
        number_pfam = db.get_number_of_pfam_entries()
    versions = db.get_versions_table()
    return render(request, 'chlamdb/home.html', my_locals(locals()))


class ComparisonIndexView(ComparisonViewMixin, View):

    def get(self, request):
        boxes = []
        for i, el in enumerate(self.available_views):
            if i % 3 == 0:
                boxes.append([])
            boxes[-1].append(el)
        context = self.get_context(boxes=boxes)
        return render(request, 'chlamdb/index_comp.html', context)


class Genomes(BaseViewMixin, View, GenomesTableMixin):

    template = 'chlamdb/genomes.html'
    view_name = "genomes"
    genome_source_object = "database"

    def get(self, request, *args, **kwargs):
        results = self.get_genomes_table()
        return render(request, self.template, self.get_context(results=results))


def extract_contigs(request, genome):
    page_title = page2title["extract_contigs"]

    taxid = int(genome)
    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)

    descr = db.get_genomes_description().description.to_dict()
    prot_infos = db.get_proteins_info(
        [taxid], search_on="taxid", as_df=True,
        to_return=["locus_tag", "product", "gene"])
    seqids = prot_infos.index.tolist()
    ogs = db.get_og_count(seqids, search_on="seqid")

    loc = db.get_gene_loc(seqids, as_hash=False).set_index("seqid")
    contigs = db.get_contigs_to_seqid(taxid)

    def lambda_format_og(og):
        return format_orthogroup(og, to_url=True, from_str=False)

    all_infos = prot_infos.join(ogs).join(loc).join(contigs)
    all_infos.gene = all_infos.gene.map(format_gene)
    all_infos.locus_tag = all_infos.locus_tag.map(format_locus)
    all_infos.orthogroup = all_infos.orthogroup.map(lambda_format_og)

    organism = descr[taxid]
    data_table_header = [
        "Gene",
        "Product",
        "Locus_tag",
        "Orthogroup",
        "Contig",
        "Strand",
        "Start",
        "Stop",
    ]
    data_table = all_infos[[
        "gene",
        "product",
        "locus_tag",
        "orthogroup",
        "contig",
        "strand",
        "start",
        "end"]].values.tolist()
    return render(request, 'chlamdb/extract_contigs.html', my_locals(locals()))


# to be moved somewhere else at some point
def to_color_code(c):
    red = int(256 * c.red)
    green = int(256 * c.green)
    blue = int(256 * c.blue)
    return f"#{red: x}{green: x}{blue: x}"


class LocusHeatmapColumn(SimpleColorColumn):
    def __init__(self, values, ref_taxon=None, header=None):
        super().__init__(values, header)
        self.ref_taxon = ref_taxon
        if len(values) > 0:
            self.min_val = min(v for k, v in values.items())
            self.max_val = max(v for k, v in values.items())

    def get_face(self, index):
        index = int(index)
        if index == self.ref_taxon:
            text_face = TextFace("-".center(11))
            text_face.inner_background.color = EteTree.GREEN
            self.set_default_params(text_face)
            return text_face

        val = self.values.get(index, None)
        if val is None:
            return TextFace("-")

        color = colors.linearlyInterpolatedColor(
            colors.gray, colors.firebrick, self.min_val, self.max_val, val)
        text_face = TextFace(str(int(val)).center(12 - len(str(int(val)))))
        text_face.inner_background.color = to_color_code(color)
        self.set_default_params(text_face)
        return text_face


def str_if_none(s):
    if s is None:
        return "-"
    return s


def search_suggest(request,):
    index = sb.ChlamdbIndex.use_index(settings.SEARCH_INDEX)
    params = request.GET
    user_query = params["term"]
    results = list(index.search(user_query, limit=None))
    data = [
        {"label": f"{i.get('name')}: {i.get('description')} ({i.get('entry_type')})",
         "value": f"{i.get('name')}: {i.get('description')}"}
        for i in results
    ]
    return JsonResponse(data, safe=False)

# NOTE: should refactor this code to avoid duplicated code


def search_bar(request):
    index = sb.ChlamdbIndex.use_index(settings.SEARCH_INDEX)
    user_query = request.GET.get("search2")

    results = index.search(user_query, limit=None)
    df = pd.DataFrame.from_records([res.fields() for res in results])
    if len(results) == 0:
        ctx = {"search_failed": True, "search_term": user_query}
        return render(request, "chlamdb/search.html", my_locals(ctx))

    # when iterating over the df, accessing the name attribute will
    # give back the name of the series, so we cannot use that attribute
    if sb.GeneEntry.entry_type in df.entry_type.values:
        df["gene"] = df.get("name")

    df = df.where(df.notna(), "-")

    tabs = []
    for entry_type_name in df.entry_type.unique():
        entry_type = sb.entry_type_to_cls[entry_type_name]()
        object_type = entry_type.object_type
        sel = df[df.entry_type == entry_type_name]
        if object_type == "locus":
            sel[object_type] = sel.locus_tag.apply(format_locus, to_url=True)
            tabs.append(TabularResultTab(
                object_type,
                "Genes",
                show_badge=True,
                table_headers=["Accession", "Gene", "Product", "Organism"],
                table_data=sel,
                table_data_accessors=[object_type, "gene", "description", "organism"]))
        else:
            metadata = MetadataGetter().object_type_to_metadata[object_type]
            if not metadata.is_enabled:
                continue
            sel[object_type] = sel.name.apply(entry_type.get_entry_id).apply(
                metadata.format_entry, to_url=True)
            tabs.append(TabularResultTab(
                object_type,
                metadata.object_name_plural,
                show_badge=True,
                table_headers=[metadata.object_name, "Description"],
                table_data=sel,
                table_data_accessors=[object_type, "description"]))

    ctx = {"search_term": user_query,
           "result_tabs": tabs,
           }
    return render(request, "chlamdb/search.html", my_locals(ctx))


def format_module(mod_id, base=None, to_url=False):
    if base is None:
        formated = f"M{mod_id:05d}"
    else:
        formated = base

    if to_url:
        return f"<a href=/KEGG_module_map/{formated}>{formated}</a>"
    return formated


class CogPhyloHeatmap(CogViewMixin, View):

    def post(self, request, *args, **kwargs):
        return render(request, 'chlamdb/cog_phylo_heatmap.html',
                      self.get_context())

    def get(self, request, frequency, *args, **kwargs):
        freq = frequency != "False"
        if freq:
            self._metadata_cls = CategoriesFreqHeatmapMetadata
        else:
            self._metadata_cls = CategoriesCountHeatmapMetadata

        tree = self.db.get_reference_phylogeny()
        descr = self.db.get_genomes_description()
        all_taxids = descr.index.tolist()

        all_cog_hits = self.db.get_cog_hits(all_taxids, search_on="taxid")
        all_cog_funcs = self.db.get_cog_summaries(
            all_cog_hits.index.unique().tolist(),
            only_cog_desc=True, as_df=True)
        all_cog_hits = all_cog_hits.join(all_cog_funcs.function)
        summed_entries = all_cog_hits.groupby("function").sum()

        # this is necessary as some cog are assigned several functions, in
        # which case, all functions are concatenated into a single string
        # e.g. MCT to say that a cog has the three functions
        for func, entries in summed_entries.iterrows():
            if len(func) == 1:
                continue
            for single_func in func:
                if single_func in summed_entries.index:
                    summed_entries.loc[single_func] += entries
                else:
                    summed_entries.loc[single_func] = entries
        grouped_by_functions = summed_entries[summed_entries.index.map(len) == 1].T

        t1 = Tree(tree)
        R = t1.get_midpoint_outgroup()
        t1.set_outgroup(R)
        t1.ladderize()

        funcs_descr = self.db.get_cog_code_description()
        e_tree = EteTree(t1)
        e_tree.rename_leaves(descr.description.to_dict())
        ttl_cnt = grouped_by_functions.sum(axis=1)
        for func in grouped_by_functions.columns:
            detailed_func = funcs_descr[func]
            func_count = grouped_by_functions[func]
            if freq:
                func_count /= ttl_cnt
                func_count *= 100
                func_count = func_count.round(2)
            col = SimpleColorColumn.fromSeries(
                func_count, header=detailed_func + "(" + func + ")",
                color_gradient=True)
            e_tree.add_column(col)

        freq = frequency
        path = settings.ASSET_ROOT + f"/temp/COG_tree_{freq}.svg"
        asset_path = f"/temp/COG_tree_{freq}.svg"
        e_tree.render(path, dpi=600)
        context = self.get_context(envoi=True, freq=freq, asset_path=asset_path)
        return render(request, 'chlamdb/cog_phylo_heatmap.html', context)


class KOModuleChooser:
    def __init__(self, hsh):
        self.hsh = hsh

    def get_color(self, index):
        return self.hsh[index]


def KEGG_module_map(request, module_name):
    page_title = page2title["KEGG_module_map"]

    if request.method != "GET":
        return render(request, 'chlamdb/KEGG_module_map.html',
                      my_locals(locals()))
    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    try:
        module_id = int(module_name[len("M"):])
    except Exception:
        context = my_locals(locals())
        context.update(errors["unknown_accession"])
        return render(request, 'chlamdb/KEGG_module_map.html',
                      my_locals(locals()))

    module_infos = db.get_modules_info([module_id])
    if len(module_infos) != 1:
        context = my_locals(locals())
        context.update(errors["no_module_info"])
        return render(request, 'chlamdb/KEGG_module_map.html', context)
    else:
        mod_id, module_descr, module_def, cat, sub_cat = module_infos[0]

    parser = ModuleParser(module_def)
    expr_tree = parser.parse()
    ko_ids = db.get_module_kos(module_id)
    mat = db.get_ko_hits(ko_ids, search_on="ko", indexing="taxid").T
    map_data = [(format_ko(ko_id), ko_desc)
                for ko_id, ko_desc in db.get_ko_desc(ko_ids).items()]
    if mat.empty:
        # should add an error message: no gene was associated for any
        # of the KO of the current module
        envoi = True
        menu = True
        valid_id = True
        return render(request, 'chlamdb/KEGG_module_map.html',
                      my_locals(locals()))

    seqids = db.get_ko_hits(ko_ids, search_on="ko", indexing="seqid")
    associated_ogs = db.get_og_count(seqids.index.tolist(), search_on="seqid")
    og_taxid = db.get_og_count(associated_ogs.orthogroup.unique().tolist(),
                               search_on="orthogroup").T
    ko_to_og_mapping = seqids.join(associated_ogs).groupby(
        ["ko"])["orthogroup"].unique()
    leaf_to_name = db.get_genomes_description().description.to_dict()

    tree = Tree(db.get_reference_phylogeny())
    R = tree.get_midpoint_outgroup()
    tree.set_outgroup(R)
    tree.ladderize()
    e_tree = EteTree(tree)

    hsh_pres = collections.defaultdict(dict)
    for ko in ko_ids:
        if ko not in mat.columns:
            e_tree.add_column(SimpleColorColumn({}, header=format_ko(ko),
                              default_val="-"))
            continue
        if ko not in ko_to_og_mapping.index:
            e_tree.add_column(SimpleColorColumn.fromSeries(mat[ko],
                              header=format_ko(ko)))
            continue

        associated_ogs = ko_to_og_mapping.loc[ko]
        n_homologs = og_taxid[associated_ogs].sum(axis=1)
        hsh_col, hsh_val = {}, {}
        for taxid, count in mat[ko].items():
            hsh_curr = hsh_pres[taxid]
            if count > 0:
                hsh_curr[ko] = 1
                hsh_val[taxid] = count
                hsh_col[taxid] = EteTree.RED
            elif taxid in n_homologs.index:
                cnt = n_homologs.loc[taxid]
                hsh_val[taxid] = cnt
                if cnt > 0:
                    hsh_col[taxid] = EteTree.GREEN
                    hsh_curr[ko] = 1
            else:
                hsh_val[taxid] = 0
        col_chooser = KOModuleChooser(hsh_col)
        e_tree.add_column(SimpleColorColumn(hsh_val, header=format_ko(ko),
                                            col_func=col_chooser.get_color))

    hsh_n_missing = {}
    for bioentry, _ in leaf_to_name.items():
        index = int(bioentry)
        if index not in mat.index:
            n_missing = expr_tree.get_n_missing({})
        else:
            n_missing = expr_tree.get_n_missing(hsh_pres[index])
        hsh_n_missing[index] = n_missing

    e_tree.add_column(SimpleColorColumn({}, default_val=" "))
    completeness = ModuleCompletenessColumn(
        hsh_n_missing, "", add_missing=False)
    e_tree.add_column(completeness)
    e_tree.rename_leaves(leaf_to_name)

    big = len(mat.columns) >= 40
    dpi = 800 if big else 1200
    path = settings.ASSET_ROOT + '/temp/KEGG_tree_%s.svg' % module_name
    asset_path = '/temp/KEGG_tree_%s.svg' % module_name
    e_tree.render(path, dpi=dpi)
    envoi = True
    menu = True
    valid_id = True
    return render(request, 'chlamdb/KEGG_module_map.html', my_locals(locals()))


def gen_pathway_profile(db, ko_ids):
    mat = db.get_ko_hits(ko_ids, search_on="ko", indexing="taxid").T
    seqids = db.get_ko_hits(ko_ids, search_on="ko", indexing="seqid")
    associated_ogs = db.get_og_count(seqids.index.tolist(), search_on="seqid")
    og_taxid = db.get_og_count(associated_ogs.orthogroup.unique().tolist(),
                               search_on="orthogroup").T
    ko_to_og_mapping = seqids.join(associated_ogs).groupby(
        ["ko"])["orthogroup"].unique()
    leaf_to_name = db.get_genomes_description().description.to_dict()

    tree = Tree(db.get_reference_phylogeny())
    R = tree.get_midpoint_outgroup()
    tree.set_outgroup(R)
    tree.ladderize()
    e_tree = EteTree(tree)
    for ko in ko_ids:
        if ko not in mat.columns:
            continue
        if ko not in ko_to_og_mapping.index:
            e_tree.add_column(SimpleColorColumn.fromSeries(mat[ko],
                              header=format_ko(ko)))
            continue

        associated_ogs = ko_to_og_mapping.loc[ko]
        n_homologs = og_taxid[associated_ogs].sum(axis=1)
        hsh_col, hsh_val = {}, {}
        for taxid, count in mat[ko].items():
            if count > 0:
                hsh_val[taxid] = count
                hsh_col[taxid] = EteTree.RED
            elif taxid in n_homologs.index:
                cnt = n_homologs.loc[taxid]
                hsh_val[taxid] = cnt
                if cnt > 0:
                    hsh_col[taxid] = EteTree.GREEN
            else:
                hsh_val[taxid] = 0
        col_chooser = KOModuleChooser(hsh_col)
        e_tree.add_column(SimpleColorColumn(hsh_val, header=format_ko(ko),
                                            col_func=col_chooser.get_color))
    e_tree.rename_leaves(leaf_to_name)
    return e_tree


def extract_map(db, request):
    if request.method != "POST":
        raise Exception("Wrong method")

    if "pathway" not in request.POST:
        raise Exception("Missing argument")

    return int(request.POST["pathway"])


def KEGG_mapp_ko(request, map_name=None, taxon_id=None):
    page_title = page2title["KEGG_mapp_ko"]

    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)

    if map_name is None:
        try:
            pathway = extract_map(db, request)
        except Exception:
            ctx = {"error": True, "error_message": "No pathway specified",
                   "error_title": "Error", "page_title": page_title}
            return render(request, 'chlamdb/KEGG_map_ko.html', my_locals(ctx))
        map_name = format_pathway(pathway)
    else:
        pathway = int(map_name[len("map"):])

    if taxon_id is not None:
        taxid = int(taxon_id)

    kos = db.get_ko_pathways([pathway], search_on="pathway", as_df=True)

    ko_list = kos["ko"].unique().tolist()
    ko_hits = db.get_ko_hits(ko_list, search_on="ko", indexing="taxid")

    if ko_hits.empty:
        ctx = {"error": True, "error_title": "No hits for this pathway",
               "page_title": page_title}
        return render(request, 'chlamdb/KEGG_map_ko.html', my_locals(ctx))

    ko_desc = db.get_ko_desc(ko_list)
    ko_ttl_count = ko_hits.sum(axis=1)
    hsh_organisms = db.get_genomes_description().description.to_dict()

    header = ["KO", "Description", "All occurrences"]
    data = []
    all_kos = []
    if taxon_id is not None:
        header.insert(2, "#in this genome")
    for ko_id, descr in ko_desc.items():
        if ko_id in ko_ttl_count.index:
            ttl = ko_ttl_count.loc[ko_id]
        else:
            ttl = 0
        if ko_id in ko_hits.index and taxon_id is not None:
            in_this_genome = ko_hits[taxid].loc[ko_id]
            if ko_hits[taxid].loc[ko_id] > 0:
                all_kos.append(format_ko(ko_id))
        else:
            in_this_genome = 0

        if taxon_id is None:
            entry = (format_ko(ko_id, as_url=True), descr, ttl)
        else:
            entry = (format_ko(ko_id, as_url=True), descr, in_this_genome, ttl)
        data.append(entry)
    e_tree = gen_pathway_profile(db, ko_list)
    path = settings.ASSET_ROOT + f"/temp/{map_name}.svg"
    asset_path = f"/temp/{map_name}.svg"
    e_tree.render(path, dpi=800)
    ctx = {"pathway_num": kos.iloc[0].pathway,
           "pathway": kos.iloc[0].description,
           "header": header,
           "data": data,
           "asset_path": asset_path,
           "url": map_name + "+" + "+".join(all_kos),
           "envoi": True,
           "page_title": page_title
           }

    if taxon_id is not None:
        ctx["organism"] = hsh_organisms[taxid]

    return render(request, 'chlamdb/KEGG_map_ko.html', my_locals(ctx))


def get_cog(request, taxon_id, category):
    biodb = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb)

    cog_hits = db.get_cog_hits([int(taxon_id)], indexing="seqid", search_on="taxid")
    cog_ids = cog_hits.cog.unique().tolist()
    cog_summaries = db.get_cog_summaries(cog_ids, only_cog_desc=True)
    prot_infos = db.get_proteins_info(cog_hits.index.unique().to_list())
    organisms = db.get_genomes_description().description.to_dict()
    functions = db.get_cog_code_description()

    data = []
    organism = organisms[int(taxon_id)]

    for seqid, cog_hit_data in cog_hits.iterrows():
        cog_id = cog_hit_data.cog
        if cog_id not in cog_summaries:
            continue
        cog_func, cog_desc = cog_summaries[cog_id]
        if category not in cog_func:
            continue
        locus_tag = prot_infos[seqid][0]
        product = prot_infos[seqid][3]
        data.append([organism, locus_tag, format_cog(cog_id), cog_desc, product])

    data_type = "cog"
    description = functions[category]
    return render(request, 'chlamdb/cog_info.html', my_locals(locals()))


def module_cat_info(request, taxid, category):
    # Not really efficient code: it would be better
    # to get the cat_id associated with category to select
    # in the category list to avoid multiple string comparison

    biodb_path = settings.BIODB_DB_PATH
    db = DB.load_db(biodb_path, settings.BIODB_CONF)

    organisms = db.get_genomes_description().description.to_dict()
    taxid = int(taxid)
    if len(organisms) == 0 or taxid not in organisms:
        return render(request, 'chlamdb/cog_info.html', my_locals(locals()))
    else:
        organism = organisms[taxid]

    category = category.replace("+", " ")
    ko_counts = db.get_ko_count([taxid], keep_seqids=True, as_multi=False)
    ko_modules = db.get_ko_modules(
        ko_counts["KO"].values.tolist(), as_pandas=True, compact=True)
    ko_modules_info = db.get_modules_info(ko_modules["module_id"].unique().tolist(),
                                          as_pandas=True)
    filtered_modules = ko_modules_info[ko_modules_info.subcat == category]
    selected_kos = filtered_modules.merge(ko_modules,
                                          left_on="module_id",
                                          right_on="module_id",
                                          how="inner")["ko_id"].unique()
    selected_seqids = ko_counts[ko_counts.KO.isin(selected_kos)]
    seqids = selected_seqids["seqid"].unique().tolist()
    hsh_to_prot = db.get_proteins_info(seqids)
    hsh_ko_desc = db.get_ko_desc(selected_kos.tolist())

    # description, locus, KO, KO name, KO description
    data = []
    for index, row in selected_seqids[["KO", "seqid"]].iterrows():
        seqid, ko_id = row.seqid, row.KO
        if seqid not in hsh_to_prot:
            continue
        locus, prot_id, gene, product = hsh_to_prot[seqid]
        ko_desc = hsh_ko_desc[ko_id]
        piece = [organism, locus, format_ko(ko_id), ko_desc, product]
        data.append(piece)
    description = category
    data_type = 'ko'
    return render(request, 'chlamdb/cog_info.html', my_locals(locals()))


def js_bioentries_to_description(hsh):
    taxon_map = 'var taxon2description = { '
    mid = ",".join(f"{to_s(bioentry)}: {to_s(description)}"
                   for bioentry, description in hsh.items())
    return taxon_map + mid + "};"


class KoBarchart(KoViewMixin, View):

    _metadata_cls = CategoriesBarchartMetadata

    @property
    def form_class(self):
        return make_venn_from(self.db, label=self.object_name_plural,
                              action="/ko_barchart/")

    def get(self, request, *args, **kwargs):
        self.form = self.form_class()
        return render(request, 'chlamdb/ko_barplot.html', self.get_context())

    def post(self, request, *args, **kwargs):
        self.form = self.form_class(request.POST)
        if not self.form.is_valid():
            return render(request, 'chlamdb/ko_barplot.html',
                          self.get_context())

        taxids = self.form.get_taxids()

        ko_counts = self.db.get_ko_count(
            taxids, keep_seqids=True, as_multi=False)
        ko_ids = ko_counts.KO.unique()
        ko_module_ids = self.db.get_ko_modules(
            ko_ids.tolist(), as_pandas=True, compact=True)
        ko_modules_info = self.db.get_modules_info(
            ko_module_ids["module_id"].unique().tolist(), as_pandas=True)

        merged = ko_counts.merge(ko_module_ids, left_on="KO",
                                 right_on="ko_id", how="inner")
        merged = merged.merge(ko_modules_info, left_on="module_id",
                              right_on="module_id", how="inner")
        cat_count = merged[["taxid", "subcat", "seqid"]].groupby(["taxid", "subcat"]).nunique()
        subcategories_list = cat_count.index.unique(level="subcat").to_list()
        subcategories = ",".join(f"{to_s(cat)}" for cat in subcategories_list)
        labels = f"[{subcategories}]"
        series_data = []

        # not ideal, but I'm really fed up with multi-indices. Be my guest
        # if you want to improve on this.
        cat_count_dict = cat_count["seqid"].to_dict()
        taxids = cat_count.index.unique(level="taxid")
        for taxid in taxids:
            entry_data = []
            for subcat in subcategories_list:
                if (taxid, subcat) in cat_count_dict:
                    entry_data.append(cat_count_dict[(taxid, subcat)])
                else:
                    entry_data.append(0)
            str_entry_data = (str(entry) for entry in entry_data)
            string = f"{{label: {to_s(taxid)}, values: [" + ",".join(str_entry_data) + "]}"
            series_data.append(string)

        taxon2description = self.db.get_genomes_description().description.to_dict()
        taxon_map = 'var taxon2description = {'
        taxon_map_lst = (f"\"{target}\": \"{taxon2description[target]}\""
                         for target in taxids)
        taxon_map = taxon_map + ",".join(taxon_map_lst) + '};'

        taxids = "?" + "&".join((f"h={i}" for i in taxids))
        series = "[" + ",".join(series_data) + "]"

        context = self.get_context(
            envoi=True, series=series, labels=labels, taxids=taxids,
            taxon_map=taxon_map)
        return render(request, 'chlamdb/ko_barplot.html', context)


class CogBarchart(CogViewMixin, View):

    _metadata_cls = CategoriesBarchartMetadata

    @property
    def form_class(self):
        return make_venn_from(self.db, label=self.object_name_plural,
                              action="/cog_barchart/")

    def get(self, request, *args, **kwargs):
        self.form = self.form_class()
        return render(request, 'chlamdb/cog_barplot.html', self.get_context())

    def post(self, request, *args, **kwargs):
        self.form = self.form_class(request.POST)
        if not self.form.is_valid():
            return render(request, 'chlamdb/cog_barplot.html',
                          self.get_context())

        target_bioentries = self.form.get_taxids()

        hsh_counts = self.db.get_cog_counts_per_category(target_bioentries)
        taxon2description = self.db.get_genomes_description().description.to_dict()
        category_dico = self.db.get_cog_code_description()

        # create a dictionnary to convert cog category description and one letter code
        category_map = 'var description2category = {'
        map_lst = (f"\"{func_descr}\": \"{func}\""
                   for func, func_descr in category_dico.items())
        category_map = category_map + ",".join(map_lst) + '};'

        taxon_map = 'var taxon2description = {'
        taxon_map_lst = (f"\"{target}\": \"{taxon2description[target]}\""
                         for target in target_bioentries)
        taxon_map = taxon_map + ",".join(taxon_map_lst) + '};'

        # Not too happy with this code and its level of indentation
        # Could also be made faster by avoiding string comparisons and list lookup
        taxon2category2count = {}
        all_categories = []
        for bioentry, hsh_cnt in hsh_counts.items():
            bioentry_str = str(bioentry)
            if bioentry_str not in taxon2category2count:
                taxon2category2count[bioentry_str] = {}

            for func, cnt in hsh_cnt.items():
                # a cog can have multiple functions
                for i in range(0, len(func)):
                    f = func[i]
                    category = category_dico[f]
                    if category in taxon2category2count[bioentry_str]:
                        taxon2category2count[bioentry_str][category] += cnt
                    else:
                        taxon2category2count[bioentry_str][category] = cnt
                        if category not in all_categories:
                            all_categories.append(category)

        labels_template = '[\n' \
                          '%s\n' \
                          ']\n'

        serie_template = '[%s\n' \
                         ']\n'

        one_serie_template = '{\n' \
                             'label: "%s",\n' \
                             'values: [%s]\n' \
                             '},\n'

        all_series_templates = []
        for taxon in taxon2category2count:
            one_category_list = []
            for category in all_categories:
                count = taxon2category2count[taxon].get(category, 0)
                one_category_list.append(count)
            one_category_list = [str(i) for i in one_category_list]

            all_series_templates.append(
                one_serie_template % (taxon, ','.join(one_category_list)))

        taxids = "?" + "&".join((f"h={i}" for i in target_bioentries))
        series = serie_template % ''.join(all_series_templates)
        labels = labels_template % ('"' + '","'.join(all_categories) + '"')
        context = self.get_context(
            envoi=True, series=series, labels=labels, taxids=taxids,
            category_map=category_map, taxon_map=taxon_map)
        return render(request, 'chlamdb/cog_barplot.html', context)


class PanGenome(ComparisonViewMixin, View):

    _metadata_cls = AccumulationRarefactionMetadata

    @property
    def form_class(self):
        return make_venn_from(self.db, label=self.object_name_plural,
                              limit=2, limit_type="lower",
                              action=f"/pan_genome/{self.object_type}")

    def get(self, request, *args, **kwargs):
        self.form = self.form_class()
        return render(request, 'chlamdb/pan_genome.html',
                      self.get_context(form=self.form))

    def post(self, request, *args, **kwargs):
        self.form = self.form_class(request.POST)
        if not self.form.is_valid():
            return render(request, 'chlamdb/pan_genome.html',
                          self.get_context(form=self.form))

        taxids = self.form.get_taxids()
        df_hits = self.get_hit_counts(taxids, search_on="taxid")
        if df_hits.empty:
            return render(request, 'chlamdb/pan_genome.html',
                          self.get_context(form=self.form,
                                           **errors["no_hits"]))
        unique, counts = np.unique(np.count_nonzero(df_hits, axis=1),
                                   return_counts=True)
        unique_to_count = dict(zip(unique, counts))

        data_count = []
        for i in range(1, len(taxids) + 1):
            count = unique_to_count.get(i, 0)
            data_count.append(count)

        acc_set = set()
        core_set = set(df_hits.index.tolist())
        sum_og = []
        core_og = []
        for col in df_hits:
            curr_col = df_hits[col]
            tmp_set = set(curr_col.index[curr_col != 0].unique())
            acc_set = acc_set.union(tmp_set)
            core_set = core_set.intersection(tmp_set)

            sum_og.append(len(acc_set))
            core_og.append(len(core_set))

        js_data_count = "[" + ",".join(str(i) for i in data_count) + "]"
        js_data_acc = "[" + ",".join(str(i) for i in sum_og) + "]"
        js_data_core = "[" + ",".join(str(i) for i in core_og) + "]"

        context = self.get_context(
            envoi=True, js_data_acc=js_data_acc,
            js_data_count=js_data_count, js_data_core=js_data_core,
            form=self.form
            )
        return render(request, 'chlamdb/pan_genome.html', context)


blast_input_dir = {"blastp": "faa",
                   "tblastn": "fna",
                   "blastn_fna": "fna",
                   "blastn_ffn": "ffn",
                   "blastx": "faa"}

blast_command = {"blastp": NcbiblastpCommandline,
                 "tblastn": NcbitblastnCommandline,
                 "blastn_fna": NcbiblastnCommandline,
                 "blastn_ffn": NcbiblastnCommandline,
                 "blastx": NcbiblastxCommandline}


def gen_blast_heatmap(db, blast_res, blast_type, no_query_name=False):
    parsed_results = NCBIXML.parse(StringIO(blast_res))

    # collects the bitscore and the query accession
    hits = collections.defaultdict(list)
    accessions = set()
    for record in parsed_results:
        if len(record.alignments) == 0:
            continue
        query = record.query
        if no_query_name:
            query = "query"
        for hit in record.alignments:
            # We cannot use hit.accession as this strips part of the id when
            # it ends with a version number (e.g. FLZO01000013.1). We recover
            # the accessions from the hit_id, which sometimes contains a source
            # in the form emb|FLZO01000013.1| and sometimes not.
            accession = hit.hit_id.rstrip("|").split("|")[-1]
            hsp = hit.hsps[0]
            scores = hits[query]
            scores.append((accession,
                           100.0 * hsp.identities / hsp.align_length))
            accessions.add(accession)

    if blast_type in ["blastp", "blastx", "blastn_ffn"]:
        acc_to_taxid = db.get_taxid_from_accession(list(accessions),
                                                   look_on="locus_tag")
    elif blast_type in ["tblastn", "blastn_fna"]:
        acc_to_taxid = db.get_taxid_from_accession(list(accessions),
                                                   look_on="contig")
    else:
        raise Exception("Unknown blast type: " + blast_type)

    all_infos = []
    for query, lst_vals in hits.items():
        hsh_taxid_to_score = {}
        for accession, score in lst_vals:
            taxid = acc_to_taxid.loc[accession].taxid
            if hsh_taxid_to_score.get(taxid, 0) < score:
                hsh_taxid_to_score[taxid] = int(round(score, 0))
        all_infos.append((query, hsh_taxid_to_score))

    min_val = min(min(hsh.values()) for _, hsh in all_infos)
    max_val = max(max(hsh.values()) for _, hsh in all_infos)
    tree = db.get_reference_phylogeny()
    descr = db.get_genomes_description()

    t1 = Tree(tree)
    R = t1.get_midpoint_outgroup()
    t1.set_outgroup(R)
    t1.ladderize()
    e_tree = EteTree(t1)
    e_tree.rename_leaves(descr.description.to_dict())

    for query, hsh_values in all_infos:
        col = SimpleColorColumn(
            hsh_values, header=query, color_gradient=True,
            default_val="-", gradient_value_range=[min_val, max_val])
        e_tree.add_column(col)

    base_file_name = time.strftime("blast_%d_%m_%y_%H_%M.svg", time.gmtime())
    path = settings.ASSET_ROOT + f"/temp/{base_file_name}"
    asset_path = f"/temp/{base_file_name}"
    e_tree.render(path, dpi=600)
    return asset_path


# for now, simplified the tblastn output to the same output as the
# other blast versions. May be re-introduced in future versions
def blast(request):
    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    page_title = page2title["blast"]
    blast_form_class = make_blast_form(db)

    if request.method != "POST":
        form = blast_form_class()
        return render(request, 'chlamdb/blast.html',
                      my_locals({"form": form, "page_title": page_title}))

    form = blast_form_class(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/blast.html',
                      my_locals({"form": form, "page_title": page_title}))

    number_blast_hits = form.cleaned_data['max_number_of_hits']
    blast_type = form.cleaned_data['blast']
    target = form.get_target()

    query_file = NamedTemporaryFile(mode='w')
    SeqIO.write(form.records, query_file, "fasta")
    query_file.flush()

    if target == 'all':
        my_db = 'merged'
    else:
        dictionary_acc_names = db.get_taxon_id_to_filenames()
        my_db = dictionary_acc_names[target]

    blast_args = {"query": query_file.name, "outfmt": 5}
    blast_args["db"] = settings.BLAST_DB_PATH + "/" + blast_input_dir[blast_type] + "/" + my_db
    if number_blast_hits != 'all':
        blast_args["max_target_seqs"] = number_blast_hits

    blast_cline = blast_command[blast_type](**blast_args)

    blastType = "locus"
    if blast_type == "tblastn" or blast_type == "blastn_fna":
        blastType = "genome"

    try:
        blast_stdout, blast_stderr = blast_cline()
    except Exception as e:
        context = {
            "error_message": "Error in blast " + str(e), "wrong_format": True,
            "error_title": "Blast error", "envoi": True, "form": form,
            "page_title": page_title}
        return render(request, 'chlamdb/blast.html', my_locals(context))

    if blast_stdout.find("No hits found") != -1:
        blast_no_hits = blast_stdout
    elif len(blast_stderr) != 0:
        blast_err = blast_stderr
    else:
        if target == "all":
            asset_path = gen_blast_heatmap(db, blast_stdout,
                                           blast_type, form.no_query_name)
        rand_id = id_generator(6)
        blast_file_l = settings.ASSET_ROOT + '/temp/%s.xml' % rand_id
        f = open(blast_file_l, 'w')
        f.write(blast_stdout)
        f.close()
        asset_blast_path = '/assets/temp/%s.xml' % rand_id
        js_out = True
        phylo_distrib = target == "all"
    envoi = True
    return render(request, 'chlamdb/blast.html', my_locals(locals()))


# Might be interesting to draw the regions in the same order as they appear
# in the phylogenetic tree (assuming single region per genome)
def coalesce_regions(genomic_regions, seqids):
    seqids_set = set(seqids)

    index = []
    for idx, (seqid, region, start, end) in enumerate(genomic_regions):
        region_seqids = set(region.index.unique())
        intersect = seqids_set & region_seqids
        index.append((idx, len(intersect), intersect))

    index.sort(key=lambda x: x[1], reverse=True)
    filtered_results = []
    for idx, _, intersect in index:
        seqid, region, start, end = genomic_regions[idx]
        if len(intersect & seqids_set) == 0:
            continue
        filtered_results.append(genomic_regions[idx])
        seqids_set = seqids_set - intersect
    return filtered_results


def optimal_region_order(regions):
    # We try to figure out an optimal order to display the sequences in.
    # For that we calculate the number of common orthogroups and whether
    # they appear in the same order. We then use that as score and calculate
    # a decent path using nearest-neighbour optimization.
    n_regions = len(regions)
    similarities = np.zeros((n_regions, n_regions))
    #  Make sure the genes are sorted, as we use the og order for the score:
    for seqid, region, start, end, _, _ in regions:
        region.sort_values("start_pos", inplace=True)

    for i, (seqid1, region1, start1, end1, _, _) in enumerate(regions):
        ogs1 = region1["orthogroup"].unique()
        neighboring_ogs1 = {(ogs1[i], ogs1[i+1]) for i in range(len(ogs1) - 1)}
        for j, (seqid2, region2, start2, end2, _, _) in enumerate(regions):
            if j <= i:
                continue
            ogs2 = region2["orthogroup"].unique()
            neighboring_ogs2 = {(ogs2[i], ogs2[i+1]) for i in range(len(ogs2) - 1)}
            score = (len(set(ogs1).intersection(ogs2)) +
                     len(neighboring_ogs1.intersection(neighboring_ogs2)))
            similarities[i, j] = similarities[j, i] = score

    region_indices = np.arange(n_regions)
    best_path = range(n_regions)
    max_score = sum(similarities[best_path[i], best_path[i + 1]] for i in range(n_regions - 1))
    for start in region_indices:
        path = [start]
        current = start
        while len(path) < n_regions:
            remaining_choices = [True if i not in path else False for i in region_indices]
            max_index = np.argmax(similarities[current, remaining_choices])
            next_choice = region_indices[remaining_choices][max_index]
            path.append(next_choice)
        score = sum(similarities[path[i], path[i + 1]] for i in range(n_regions - 1))
        if score > max_score:
            max_score = score
            best_path = path
    return best_path


def plot_region(request):
    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    page_title = page2title["plot_region"]
    form_class = make_plot_form(db)

    if request.method != "POST":
        form = form_class()
        return render(request, 'chlamdb/plot_region.html',
                      my_locals({"form": form, "page_title": page_title}))

    form = form_class(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/plot_region.html',
                      my_locals({"form": form, "page_title": page_title}))

    seqid = form.get_seqid()
    genomes = form.get_genomes()
    all_homologs = form.get_all_homologs()
    hsh_loc = db.get_gene_loc([seqid])
    ids = db.get_og_identity(ref_seqid=seqid)
    all_seqids = ids.index.unique().tolist()
    all_seqids.append(seqid)
    organisms = db.get_organism(all_seqids, as_df=True, as_taxid=True)
    all_infos = organisms.join(ids, how="inner")
    ref_strand, _, _ = hsh_loc[seqid]
    hsh_description = db.get_genomes_description().description.to_dict()

    if not all_homologs:
        best_matches = all_infos.groupby(["taxid"]).idxmax()
        ref_taxid = organisms.loc[seqid].taxid

        # redundant: do not include the genome of the locus tag, even if listed
        # by the user. The locus tag is its best match.
        if ref_taxid in genomes:
            genomes.remove(ref_taxid)
        best_matches = best_matches.reindex(genomes).dropna()
        if not best_matches.empty:
            seqids = best_matches["identity"].astype(int).tolist()
        else:
            seqids = []
        seqids.append(seqid)
    else:
        selection = all_infos[all_infos.taxid.isin(genomes)]
        seqids = selection.index.astype(int).tolist()
        seqids.append(seqid)

    if len(seqids) > 20:
        context = {
            "form": form,
            "error": True,
            "errors": ["Too many regions to display"],
            "page_title": page_title
        }
        return render(request, 'chlamdb/plot_region.html', my_locals(context))

    region_size = form.get_region_size()

    locus_tags = db.get_proteins_info(
        seqids, to_return=["locus_tag"], as_df=True)
    to_highlight = locus_tags["locus_tag"].tolist()
    all_regions = []
    connections = []
    prev_infos = None
    all_identities = []

    genomic_regions = []
    for seqid in seqids:
        region, start, end, contig_size, contig_topology = locusx_genomic_region(
            db, int(seqid), region_size / 2)
        genomic_regions.append([seqid, region, start, end, contig_size, contig_topology])

    # remove overlapping regions (e.g. if two matches are on the same
    # region, avoid displaying this region twice).
    # XXX : coalesce_regions does not do what it is intended to do, to fix
    filtered_regions = genomic_regions  # coalesce_regions(genomic_regions, seqids)

    # Let's add the ogs to the regions info
    for genomic_region in filtered_regions:
        seqid, region, start, end, _, _ = genomic_region
        if region["strand"].loc[seqid] * ref_strand == -1:
            mean_val = (end + start) / 2
            region["strand"] *= -1
            dist_vector_start = region["start_pos"]-mean_val
            len_vector = region["end_pos"]-region["start_pos"]
            region["end_pos"] = region["start_pos"]-2*dist_vector_start
            region["start_pos"] = region["end_pos"] - len_vector
        ogs = db.get_og_count(region.index.tolist(), search_on="seqid")
        # need to reset index to keep seqid in the next merge
        genomic_region[1] = region.join(ogs).reset_index()

    # determine the optimal prepresentation order and sort regions accordingly.
    best_path = optimal_region_order(filtered_regions)
    filtered_regions = [filtered_regions[i] for i in best_path]

    for seqid, region, start, end, contig_size, contig_topology in filtered_regions:
        if prev_infos is not None:
            # BM: horrible code (I wrote it, I should know).
            # would be nice to refactor it in a more efficient and clean way.

            # We need to drop na from orthogroup, as some genes, like rRNA
            # or tRNA are not assigned any orthogroup
            common_og = region.dropna(subset=["orthogroup"]).merge(
                    prev_infos, on="orthogroup")[
                        ["locus_tag_x",
                         "locus_tag_y",
                         "seqid_x",
                         "seqid_y",
                         "orthogroup"]]
            related = []
            ogs = common_og.orthogroup.astype(int).tolist()
            p1 = common_og.seqid_x.tolist()
            p2 = common_og.seqid_y.tolist()
            identities = db.get_og_identity(og=ogs, pairs=zip(p1, p2))
            identities = identities.set_index(["seqid_x", "seqid_y"]).identity.to_dict()
            hsh_agg = {}
            for i, v in common_og.iterrows():
                if v.seqid_x == v.seqid_y:
                    ident = 100
                elif (v.seqid_x, v.seqid_y) in identities:
                    ident = identities[(v.seqid_x, v.seqid_y)]
                else:
                    ident = identities[(v.seqid_y, v.seqid_x)]
                all_identities.append(ident)
                og_val = to_s(int(v.orthogroup))
                arr = hsh_agg.get(v.locus_tag_x, [])
                arr.append(f"[{to_s(v.locus_tag_y)}, {og_val}, {ident: .2f}]")
                hsh_agg[v.locus_tag_x] = arr
            related = (f"{to_s(loc)}: [" + ",".join(values) + "]"
                       for loc, values in hsh_agg.items())
            connections.append("{" + ",".join(related) + "}")

        taxid = organisms.loc[seqid].taxid
        genome_name = hsh_description[taxid]
        contig_name = db.get_bioentry_qualifiers(int(region["bioentry_id"][0]))\
                        .set_index("term").loc["accessions"].value
        region_name = f"{genome_name} - {contig_name} - {int(start)}:{int(end)}"
        js_val = genomic_region_df_to_js(region, start, end, contig_size,
                                         contig_topology, region_name)
        all_regions.append(js_val)
        prev_infos = region[["orthogroup", "locus_tag", "seqid"]]

    ctx = {"form": form,
           "genomic_regions": "[" + "\n,".join(all_regions) + "]",
           "to_highlight": to_highlight,
           "envoi": True,
           "connections": "[" + ",".join(connections) + "]",
           "page_title": page_title}
    if len(all_identities) > 0:
        ctx["max_ident"] = max(all_identities)
        ctx["min_ident"] = min(all_identities)
    return render(request, 'chlamdb/plot_region.html', my_locals(ctx))


def circos_main(request):
    biodb_path = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb_path)
    if request.method == 'POST':

        reference_taxon = request.POST["reference_taxid"]
        include_taxids = eval(request.POST["include_taxids"])
        exclude_taxids = eval(request.POST["exclude_taxids"])
        og_list = eval(request.POST["og_list"])

        target_taxons = include_taxids + exclude_taxids

        target_taxons.pop(target_taxons.index(int(reference_taxon)))

        js_code = get_circos_data(int(reference_taxon),
                                  [int(i) for i in target_taxons],
                                  highlight_og=og_list)

        envoi = True
        envoi_region = True

        return render(request, 'chlamdb/circos_main.html', my_locals(locals()))

    return render(request, 'chlamdb/circos_main.html', my_locals(locals()))


def get_circos_data(reference_taxon, target_taxons, highlight_og=False):
    biodb_path = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb_path)

    if highlight_og:
        df_genes = db.get_genes_from_og(
            highlight_og, taxon_ids=[reference_taxon], terms=["locus_tag"])
        locus_list = df_genes["locus_tag"].to_list()
    else:
        locus_list = []

    # "bioentry_id", "accession" ,"length"
    df_bioentry = db.get_bioentry_list(reference_taxon, min_bioentry_length=1000)

    # "bioentry_id", "seqfeature_id", "start_pos", "end_pos", "strand"
    df_feature_location = db.get_features_location(
        reference_taxon, search_on="taxon_id",
        seq_term_names=["CDS", "rRNA", "tRNA"]).set_index(["seqfeature_id"])

    # retrieve n_orthologs of list of seqids

    seq_og = db.get_og_count(
        df_feature_location.index.to_list(), search_on="seqid")
    count_all_genomes = db.get_og_count(
        seq_og["orthogroup"].to_list(), search_on="orthogroup")
    orthogroup2count_all = count_all_genomes[count_all_genomes > 0].count(axis=1)
    homologs_count = df_feature_location.loc[df_feature_location.term_name == 'CDS']\
        .join(seq_og)\
        .reset_index()\
        .set_index("orthogroup")\
        .merge(orthogroup2count_all.rename('value'),
               left_index=True,
               right_index=True)[["bioentry_id", "start_pos", "end_pos", "value"]]

    df_identity = db.get_identity_closest_homolog(
        reference_taxon, target_taxons).set_index(["target_taxid"])

    c = circosjs.CircosJs()

    c.add_contigs_data(df_bioentry)

    # sort taxons by number of homologs (from mot similar to most dissmilar)
    target_taxon_n_homologs = df_identity.groupby(["target_taxid"])\
                                         .count()["seqfeature_id_1"]\
                                         .sort_values(ascending=False)
    locus2seqfeature_id = db.get_hsh_locus_to_seqfeature_id()
    seqfeature_id2locus_tag = {
        value: key for key, value in locus2seqfeature_id.items()}
    # "bioentry_id", "seqfeature_id", "start_pos", "end_pos", "strand"
    # "seqfeature_id_1", "seqfeature_id_2", "identity", "target_taxid"
    # join on seqfeature id
    df_feature_location["gene"] = df_feature_location["qualifier_value_gene"].fillna("-")
    df_feature_location["gene_product"] = df_feature_location["qualifier_value_product"].fillna("-")
    df_feature_location["locus_tag"] = df_feature_location["qualifier_value_locus_tag"].fillna("-")

    df_feature_location["color"] = "grey"
    df_feature_location["color"][df_feature_location["term_name"] == 'tRNA'] = "magenta"
    df_feature_location["color"][df_feature_location["term_name"] == 'rRNA'] = "magenta"

    df_feature_location = df_feature_location.rename(columns={"locus_tag": "locus_ref"})
    # = [seqfeature_id2locus_tag[seqfeature_id] for seqfeature_id in df_feature_location.index]
    minus_strand = df_feature_location.set_index("strand").loc[-1]
    plus_strand = df_feature_location.set_index("strand").loc[1]
    c.add_histogram_track(homologs_count, "n_genomes",
                          radius_diff=0.1, sep=0.005, outer=True)
    c.add_gene_track(minus_strand, "minus", sep=0,
                     radius_diff=0.04, highlight_list=locus_list)
    c.add_gene_track(plus_strand, "plus", sep=0,
                     radius_diff=0.04, highlight_list=locus_list)

    # iterate ordered list of target taxids, add track to circos
    for n, target_taxon in enumerate(target_taxon_n_homologs.index):
        df_combined = df_feature_location.join(
            df_identity.loc[target_taxon].reset_index()
                                         .set_index("seqfeature_id_1")
            ).reset_index()
        df_combined.identity = df_combined.identity.fillna(0).astype(int)
        df_combined.bioentry_id = df_combined.bioentry_id.astype(str)

        # only keep the highest identity for each seqfeature id
        df_combined = df_combined.sort_values('identity', ascending=False)\
                                 .drop_duplicates('index')\
                                 .sort_index()
        loci = []
        for seqid in df_combined.seqfeature_id_2:
            # ugly hack... to be fixed
            if pd.isna(seqid) or int(seqid) not in seqfeature_id2locus_tag:
                loci.append(None)
                continue
            loci.append(seqfeature_id2locus_tag[seqid])
        df_combined["locus_tag"] = loci

        # comp is a custom scale with missing orthologs coloured in light blue
        if n == 0:
            sep = 0.03
        else:
            sep = 0.01
        c.add_heatmap_track(df_combined, f"target_{target_taxon}",
                            color="comp", sep=sep, radius_diff=0.04)

    c.add_line_track(df_bioentry, "GC_content", radius_diff=0.12,
                     fillcolor="green")

    js_code = c.get_js_code()
    return js_code


def circos(request):
    biodb_path = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb_path)
    page_title = page2title["circos"]

    circos_form_class = make_circos_form(db)

    if request.method == 'POST':

        form = circos_form_class(request.POST)

        if form.is_valid():

            target_taxons = form.get_target_taxids()
            reference_taxon = form.get_ref_taxid()

            js_code = get_circos_data(reference_taxon, target_taxons)

            envoi = True
    else:
        form = circos_form_class()

    local_vars = my_locals(locals())
    return render(request, 'chlamdb/circos.html', local_vars)


def alignment(request, input_fasta):

    handle = open(input_fasta, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        pass
    return render(request, 'chlamdb/alignment.html', my_locals(locals()))


class PlotHeatmap(ComparisonViewMixin, View):

    _metadata_cls = HeatmapMetadata

    @property
    def form_class(self):
        return make_venn_from(self.db, label=self.object_name_plural,
                              limit=2, limit_type="lower",
                              action=f"/plot_heatmap/{self.object_type}")

    def get(self, request, *args, **kwargs):
        self.form = self.form_class()
        return render(request, 'chlamdb/plot_heatmap.html',
                      self.get_context(form=self.form))

    def post(self, request, *args, **kwargs):
        self.form = self.form_class(request.POST)
        if not self.form.is_valid():
            return render(request, 'chlamdb/plot_heatmap.html',
                          self.get_context(form=self.form))

        import plotly.graph_objects as go
        import scipy.cluster.hierarchy as shc
        from scipy.cluster import hierarchy

        taxon_ids = self.form.get_taxids()
        mat = self.get_hit_counts(taxon_ids, search_on="taxid")
        if mat.empty:
            ctx = self.get_context(**errors["no_hits"],
                                   form=self.form)
            return render(request, 'chlamdb/plot_heatmap.html', ctx)

        mat.index = [self.format_entry(i, to_url=True) for i in mat.index]

        target2description = self.db.get_genomes_description().description.to_dict()
        mat.columns = (target2description[i] for i in mat.columns.values)
        # reorder row and columns based on clustering
        Z_rows = shc.linkage(mat.T, method='ward')
        order_rows = hierarchy.leaves_list(Z_rows)
        new_index = [mat.columns.values[i] for i in order_rows]

        Z_genomes = shc.linkage(mat, method='ward')
        order_genomes = hierarchy.leaves_list(Z_genomes)

        # set number of paralogs >1 as 2 to simplify the color code
        mat[mat > 1] = 2
        colors = ["#ffffff", "#2394d9", "#d923ce"]
        new_cols = [mat.index.tolist()[i] for i in order_genomes]

        new_mat = mat.T.reindex(new_index)[new_cols]
        fig = go.Figure(data=go.Heatmap(z=new_mat, colorscale=colors,
                                        y=new_mat.index, x=new_mat.columns))
        fig.update_traces(showlegend=False, showscale=False)

        html_plot = make_div(fig, div_id="heatmap")
        context = self.get_context(envoi_heatmap=True, html_plot=html_plot,
                                   form=self.form)
        return render(request, 'chlamdb/plot_heatmap.html', context)


def format_pathway(path_id, base=None, to_url=False, taxid=None):
    if base is None:
        base_string = f"map{path_id:05d}"
    else:
        base_string = base

    if not to_url:
        return base_string
    if taxid is None:
        to_page = f"/KEGG_mapp_ko/{base_string}"
    else:
        to_page = f"/KEGG_mapp_ko/{base_string}/{taxid}"
    return f"<a href=\"{to_page}\">{base_string}</a>"


def kegg_genomes(request):

    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    single_genome_form = make_single_genome_form(db)
    hsh_organisms = db.get_genomes_description().description.to_dict()
    if request.method != "POST":
        form = single_genome_form()
        return render(request, 'chlamdb/kegg_genomes.html',
                      my_locals(locals()))

    form = single_genome_form(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/kegg_genomes.html',
                      my_locals(locals()))

    taxid = form.get_genome()
    ko_hits = db.get_ko_hits([taxid], search_on="taxid")
    kos = ko_hits.index.tolist()
    ko_pathways = db.get_ko_pathways(kos)
    genomes = db.get_genomes_description().description.to_dict()

    header = ["Pathway", "Description", hsh_organisms[taxid]]

    page_title = page2title["kegg_genomes"] + f': {hsh_organisms[taxid]}'

    data = []
    pathway_count = collections.Counter()
    hsh_path_to_descr = {}
    for ko, pathways in ko_pathways.items():
        pathway_count.update(path_id for path_id, pathway in pathways)
        hsh_path_to_descr.update(pathways)

    for element, count in pathway_count.items():
        descr = hsh_path_to_descr[element]
        entry = (format_pathway(element, to_url=True, taxid=taxid),
                 descr, count)
        data.append(entry)

    ctx = {"envoi": True,
           "data": data,
           "header": header,
           "organism": hsh_organisms[taxid],
           "page_title": page_title}
    return render(request, 'chlamdb/kegg_genomes.html', my_locals(ctx))


def kegg_genomes_modules(request):
    page_title = page2title["kegg_genomes_modules"]

    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    single_genome_form = make_single_genome_form(db)
    hsh_organisms = db.get_genomes_description().description.to_dict()
    if request.method != "POST":
        form = single_genome_form()
        return render(request, 'chlamdb/kegg_genomes_modules.html',
                      my_locals(locals()))

    form = single_genome_form(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/kegg_genomes_modules.html',
                      my_locals(locals()))

    taxid = form.get_genome()
    list_taxid = [""]
    list_taxid.append(taxid)
    module_hits = db.get_ko_count_cat(taxon_ids=list_taxid)
    grouped_count = module_hits.groupby(["taxon_id", "module_id"]).sum()
    grouped_count = grouped_count.unstack(level="taxon_id")
    grouped_count.columns = grouped_count["count"].columns
    grouped_count.fillna(0, inplace=True, downcast="infer")
    taxids = grouped_count.columns
    genomes = db.get_genomes_description().description.to_dict()
    modules_info = db.get_modules_info(grouped_count.index.tolist(),
                                       as_pandas=True)
    all_infos = modules_info.set_index("module_id").join(grouped_count)

    header = ["Module", "Category", "Sub-category", "Description"]
    entries = []

    taxons = []
    for taxid in taxids:
        taxons.append(genomes[taxid])

    for module_id, data in all_infos.iterrows():
        line = [format_module(module_id, to_url=True), data["cat"],
                data.subcat, data.descr]
        for taxid in taxids:
            line.append(data[taxid])
        entries.append(line)
    envoi_comp = True
    return render(request, 'chlamdb/kegg_genomes_modules.html',
                  my_locals(locals()))


def kegg(request):
    page_title = page2title["kegg"]

    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    module_overview_form = make_module_overview_form(db)
    form_cat = module_overview_form()

    single_genome_form = make_single_genome_form(db, "kegg_genomes_modules")
    form_genome_modules = single_genome_form()

    single_genome_form = make_single_genome_form(db, "kegg_genomes")
    form_genome = single_genome_form()

    module_overview_form = make_module_overview_form(db, True)
    form_subcat = module_overview_form()

    comp_metabo_form = make_metabo_from(db, action="module_comparison")
    form_module = comp_metabo_form()

    pathway_overview_form = make_pathway_overview_form(db)
    form_pathway = pathway_overview_form()

    return render(request, 'chlamdb/kegg.html', my_locals(locals()))


def kegg_module_subcat(request):
    page_title = page2title["kegg_module_subcat"]

    biodb = settings.BIODB_DB_PATH
    db = DB.load_db(biodb, settings.BIODB_CONF)

    module_overview_form = make_module_overview_form(db, True)
    if request.method != "POST":
        form = module_overview_form()
        return render(request, 'chlamdb/module_subcat.html',
                      my_locals(locals()))

    form = module_overview_form(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/module_subcat.html',
                      my_locals(locals()))

    cat = form.cleaned_data["subcategory"]
    ko_count_subcat = db.get_ko_count_cat(subcategory=cat)

    modules_infos = db.get_modules_info(ids=[cat], search_on="subcategory")

    subcat_name = modules_infos[0][-1]
    cat_count = db.get_ko_count_cat(category=cat)
    leaf_to_name = db.get_genomes_description()
    grouped_count = ko_count_subcat.groupby(["taxon_id", "module_id"]).sum()
    if grouped_count.empty:
        n_occurences = pd.DataFrame()
    else:
        n_occurences = grouped_count.groupby(["module_id"]).count()

    data = []
    expression_tree = {}
    for module_id, descr, definition, *other in modules_infos:
        occurences = 0
        if module_id in n_occurences.index:
            occurences = n_occurences["count"].loc[module_id]
            parser = ModuleParser(definition)
            expression_tree[module_id] = parser.parse()
        data.append((format_module(module_id, to_url=True), descr, occurences))

    unique_module_ids = ko_count_subcat.index.get_level_values(
        "module_id").unique().tolist()
    if len(unique_module_ids) == 0:
        envoi = True
        context = my_locals(locals())
        context.update(errors["no_hits"])
        return render(request, 'chlamdb/module_subcat.html', context)

    grouped_count = grouped_count.unstack(level=1, fill_value=0)
    grouped_count.columns = [col for col in grouped_count["count"].columns]
    ref_tree = db.get_reference_phylogeny()
    e_tree = EteTree.default_tree(ref_tree)
    for module, counts in grouped_count.items():
        values = counts.to_dict()
        header = format_module(module)
        expr_tree = expression_tree[module]
        n_missing = {}
        for bioentry, count in counts.items():
            # hack for now
            s = ko_count_subcat.loc[bioentry].index.get_level_values("KO").tolist()
            n_missing[bioentry] = expr_tree.get_n_missing({ko: 1 for ko in s})
        new_col = KOAndCompleteness(values, n_missing, header)
        e_tree.add_column(new_col)
    e_tree.rename_leaves(leaf_to_name.description.to_dict())
    path = settings.ASSET_ROOT + '/temp/metabo_tree.svg'
    asset_path = '/temp/metabo_tree.svg'
    e_tree.render(path, dpi=500, w=800)
    envoi = True
    return render(request, 'chlamdb/module_subcat.html', my_locals(locals()))


def kegg_module(request):
    page_title = page2title["kegg_module"]

    biodb = settings.BIODB_DB_PATH
    db = DB.load_db(biodb, settings.BIODB_CONF)
    module_overview_form = make_module_overview_form(db)
    if request.method != "POST":
        form = module_overview_form()
        return render(request, 'chlamdb/module_overview.html',
                      my_locals(locals()))

    form = module_overview_form(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/module_overview.html',
                      my_locals(locals()))

    cat = form.cleaned_data["category"]
    modules_infos = db.get_modules_info(ids=[cat], search_on="category")
    cat_count = db.get_ko_count_cat(category=cat)
    leaf_to_name = db.get_genomes_description()
    cat_name_id = db.get_module_categories(module_ids=cat)
    cat_name = [name for id, name in cat_name_id][0]
    grouped_count = cat_count.groupby(["taxon_id", "module_id"]).sum()
    if grouped_count.empty:
        n_occurences = pd.DataFrame()
    else:
        n_occurences = grouped_count.groupby(["module_id"]).count()

    data = []
    expression_tree = {}
    for module_id, descr, definition, *other in modules_infos:
        occurences = 0
        if module_id in n_occurences.index:
            occurences = n_occurences["count"].loc[module_id]
            parser = ModuleParser(definition)
            expression_tree[module_id] = parser.parse()
        data.append((format_module(module_id, to_url=True), descr, occurences))

    if grouped_count.empty:
        envoi = True
        context = my_locals(locals())
        context.update(errors["no_hits"])
        return render(request, 'chlamdb/module_overview.html', context)

    grouped_count = grouped_count.unstack(level=1, fill_value=0)
    grouped_count.columns = [col for col in grouped_count["count"].columns]
    ref_tree = db.get_reference_phylogeny()
    e_tree = EteTree.default_tree(ref_tree)
    for module, counts in grouped_count.items():
        values = counts.to_dict()
        header = format_module(module)
        expr_tree = expression_tree[module]
        n_missing = {}
        for bioentry, count in counts.items():
            # hack for now
            s = cat_count.loc[bioentry].index.get_level_values("KO").tolist()
            n_missing[bioentry] = expr_tree.get_n_missing({ko: 1 for ko in s})
        new_col = KOAndCompleteness(values, n_missing, header)
        e_tree.add_column(new_col)
    e_tree.rename_leaves(leaf_to_name.description.to_dict())
    path = settings.ASSET_ROOT + '/temp/metabo_tree.svg'
    asset_path = '/temp/metabo_tree.svg'
    e_tree.render(path, dpi=500, w=800)

    envoi = True
    return render(request, 'chlamdb/module_overview.html', my_locals(locals()))


def module_comparison(request):
    page_title = page2title["module_comparison"]

    biodb = settings.BIODB_DB_PATH
    db = DB.load_db(biodb, settings.BIODB_CONF)
    comp_metabo_form = make_metabo_from(db)

    if request.method != "POST":
        form = comp_metabo_form()
        return render(request, 'chlamdb/module_comp.html', my_locals(locals()))

    form = comp_metabo_form(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/module_comp.html', my_locals(locals()))

    taxids = form.get_choices()

    module_hits = db.get_ko_count_cat(taxon_ids=taxids)
    tipo = type(taxids)
    grouped_count = module_hits.groupby(["taxon_id", "module_id"]).sum()
    grouped_count = grouped_count.unstack(level="taxon_id")
    grouped_count.columns = grouped_count["count"].columns
    grouped_count.fillna(0, inplace=True, downcast="infer")
    taxids = grouped_count.columns
    genomes = db.get_genomes_description().description.to_dict()
    modules_info = db.get_modules_info(grouped_count.index.tolist(),
                                       as_pandas=True)
    all_infos = modules_info.set_index("module_id").join(grouped_count)

    header = ["Module", "Category", "Sub-category", "Description"]
    entries = []

    taxons = []
    for taxid in taxids:
        taxons.append(genomes[taxid])

    for module_id, data in all_infos.iterrows():
        line = [format_module(module_id, to_url=True), data["cat"],
                data.subcat, data.descr]
        for taxid in taxids:
            line.append(data[taxid])
        entries.append(line)
    envoi_comp = True
    return render(request, 'chlamdb/module_comp.html', my_locals(locals()))


def faq(request):
    a = 2
    return render(request, 'chlamdb/FAQ.html', my_locals(locals()))


def phylogeny(request):
    page_title = page2title["phylogeny_intro"]

    biodb_path = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb_path)

    genomes_data = get_genomes_data(db)

    asset_path = "/temp/species_tree.svg"
    path = settings.ASSET_ROOT + '/temp/species_tree.svg'

    core = db.get_n_orthogroups(only_core=True)

    # plot phylo only of not already in assets
    tree = db.get_reference_phylogeny()
    t1 = Tree(tree)
    R = t1.get_midpoint_outgroup()
    if R is not None:
        t1.set_outgroup(R)
    t1.ladderize()

    e_tree = EteTree(t1)
    header_params = {"rotation": -30}
    stacked_face_params = {"margin_right": 8, "margin_left": 5}

    tree_params = [
        # serie_name, header, color and is_relative
        ["length", "Size (Mbp)", "#91bfdb", True],
        ["gc", "GC %", "#fc8d59", False],
        ["coding_density", "Coding density %", "#99d594", False],
        ["completeness", "Completeness %", "#d7191c", False],
        ["contamination", "Contamination %", "black", False]]

    for serie_name, header, col, is_relative in tree_params:
        data = genomes_data[serie_name]
        e_tree.add_column(SimpleColorColumn.fromSeries(
            data, header=None, use_col=False))
        stack = StackedBarColumn(data.to_dict(), colours=[col, "white"],
                                 relative=is_relative, header=header,
                                 header_params=header_params,
                                 face_params=stacked_face_params)
        e_tree.add_column(stack)

    e_tree.rename_leaves(genomes_data.description.to_dict())
    e_tree.render(path, dpi=500)
    return render(request, 'chlamdb/phylogeny_intro.html', my_locals(locals()))
