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
import matplotlib.colors as mpl_col
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import (NcbiblastnCommandline,
                                    NcbiblastpCommandline,
                                    NcbiblastxCommandline,
                                    NcbitblastnCommandline)
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from django.conf import settings
from django.http import JsonResponse
from django.shortcuts import render
from django.views import View
from ete3 import SeqMotifFace, StackedBarFace, TextFace, Tree, TreeStyle

import chlamdb.circosjs as circosjs

from lib import search_bar as sb
from lib.db_utils import (DB, NoPhylogenyException)
from lib.ete_phylo import (Column, EteTree, KOAndCompleteness,
                                       ModuleCompletenessColumn,
                                       SimpleColorColumn)
from lib.KO_module import ModuleParser

from reportlab.lib import colors

from chlamdb.forms import (make_blast_form, make_circos_form,
                           make_extract_form, make_metabo_from,
                           make_module_overview_form,
                           make_pathway_overview_form, make_plot_form,
                           make_single_genome_form, make_venn_from)
from chlamdb.utils import safe_replace

# could also be extended to cache the results of frequent queries
# (e.g. taxid -> organism name) to avoid db queries
with DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF) as db:
    hsh_config = db.get_config_table(ret_mandatory=True)
    optional2status = {name: value for name,
                       (mandatory, value) in hsh_config.items() if not mandatory}
    missing_mandatory = [name for name, (mandatory, value) in hsh_config.items()
                         if mandatory and not value]

page2title = {
    'extract_orthogroup': 'Comparisons: orthologous groups',
    'extract_pfam': 'Comparisons: PFAM domains',
    'extract_ko': 'Comparisons: Kegg Orthologs (KO)',
    'extract_cog': 'Comparisons: Clusters of Orthologous groups (COGs)',
    'venn_orthogroup': 'Comparisons: orthologous groups',
    'venn_ko': 'Comparisons: Kegg Orthologs (KO)',
    'venn_pfam': 'Comparisons: Pfam domains',
    'venn_cog': 'Comparisons: Clusters of Orthologous groups (COGs)',
    'heatmap_orthogroup': 'Comparisons: orthologous groups',
    'heatmap_pfam': 'Comparisons: PFAM domains',
    'heatmap_ko': 'Comparisons: Kegg Orthologs (KO)',
    'heatmap_COG': 'Comparisons: Clusters of Orthologous groups (COGs)',
    'pan_genome_orthology': 'Comparisons: orthologous groups',
    'pan_genome_Pfam': 'Comparisons: PFAM domains',
    'pan_genome_ko': 'Comparisons: Kegg Orthologs (KO)',
    'pan_genome_COG': 'Comparisons: Clusters of Orthologous groups (COGs)',
    'orthogroup_comparison': 'Comparisons: orthologous groups',
    'plot_heatmap_orthology': 'Comparisons: orthologous groups',
    'plot_heatmap_ko': 'Comparisons: KEGG orthologs',
    'plot_heatmap_Pfam': 'Comparisons: Pfam domains ',
    'plot_heatmap_COG': 'Comparisons: Clusters of Orthologous groups (COGs)',
    'index_comp_ko': 'Comparisons: KEGG orthologs',
    'index_comp_Pfam': 'Comparisons: Pfam domains ',
    'index_comp_COG': 'Comparisons: Clusters of Orthologous groups (COGs) ',
    'index_comp_orthology': 'Comparisons: orthologous groups',
    'cog_barchart': 'Comparisons: Clusters of Orthologous groups (COGs) ',
    'pan_genome_orthology': 'Comparisons: orthologous groups',
    'entry_list_pfam': 'Comparisons: PFAM domains',
    'entry_list_ko': 'Comparisons: Kegg Orthologs (KO)',
    'entry_list_cog': 'Comparisons: Clusters of Orthologous groups (COGs)',
    'ko_comparison': 'Comparisons: Kegg Orthologs (KO)',
    'pfam_comparison': 'Comparisons: PFAM domains',
    'amr_gene_comparison': 'Comparisons: Antimicrobial Resistance',
    'amr_class_comparison': 'Comparisons: Antimicrobial Resistance',
    'amr_subclass_comparison': 'Comparisons: Antimicrobial Resistance',
    'module_barchart': 'Comparisons: Kegg Orthologs (KO)',
    'blast': 'Homology search: Blast',
    'plot_region': 'Genome alignments: Plot region',
    'circos': 'Genome alignments: Circos plot',
    'cog_comparison': 'Comparisons: Clusters of Orthologous groups (COGs)',
    'kegg': 'Metabolism: kegg based',
    'kegg_genomes': 'Kegg metabolic pathways ',
    'kegg_genomes_modules': 'Metabolism: kegg based',
    'KEGG_mapp_ko': 'Metabolism: kegg based',
    'kegg_module': 'Metabolism: kegg based',
    'kegg_module_subcat': 'Metabolism: kegg based',
    'module_comparison': 'Metabolism: kegg based',
    'fam_pfam': 'Pfam domain ',
    'fam_ko': 'Kegg Ortholog ',
    'fam_cog': 'COG Ortholog ',
    'KEGG_module_map': 'Kegg module ',
    'phylogeny_intro': 'Phylogeny',
    'genomes_intro': 'Genomes: table of contents',
    'extract_contigs': 'Genomes: table of contents',
    'COG_phylo_heatmap': 'Comparisons: Clusters of Orthologous groups (COGs)',
    'genomes': 'Table of content: Genomes '

}


def my_locals(local_dico):
    local_dico["optional2status"] = optional2status
    local_dico["missing_mandatory"] = missing_mandatory
    return local_dico


def id_generator(size=6, chars=string.ascii_uppercase + string.ascii_lowercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))


def help(request):
    return render(request, 'chlamdb/help.html', my_locals(locals()))


def about(request):
    path = settings.BASE_DIR + '/assets/bibliography/references.bib'
    with open(path) as bibtex_file:
        bib_database = bibtexparser.load(bibtex_file)

    entry_list = []

    for entry in bib_database.entries:
        string = ("<b>%s</b><br> %s, %s, %s(%s):%s, %s" % (re.sub('[{}]', '', entry["title"]),
                                                           entry["author"],
                                                           entry["journal"],
                                                           entry["volume"],
                                                           entry["number"],
                                                           entry["pages"],
                                                           entry["year"],
                                                            ))
        url = entry["url"]
        entry_list.append([string, url])

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

    genomes_data = db.get_genomes_infos()
    genomes_descr = db.get_genomes_description()

    asset_path = "/temp/species_tree.svg"
    path = settings.BASE_DIR + '/assets/temp/species_tree.svg'

    genomes_data = genomes_data.join(genomes_descr)

    genomes_data.gc = genomes_data.gc.apply(round)
    genomes_data.coding_density = genomes_data.coding_density.apply(
        lambda x: round(100 * x))
    genomes_data.length = genomes_data.length.apply(
        lambda x: round(x / pow(10, 6), 2))

    data_table_header = ["Name", "%GC", "N proteins",
                         "N contigs", "Size (Mbp)", "Percent coding"]
    data_table = genomes_data[["description", "gc", "n_prot",
                               "n_contigs", "length", "coding_density"]].values.tolist()

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
        ["completeness", "Completeness", "#d7191c", False],
        ["contamination", "Contamination", "black", False]]

    for serie_name, header, col, is_relative in tree_params:
        data = genomes_data[serie_name]
        e_tree.add_column(SimpleColorColumn.fromSeries(
            data, header=None, use_col=False))
        stack = StackedBarColumn(data.to_dict(), colours=[col, "white"],
                                 relative=is_relative, header=header,
                                 header_params=header_params, face_params=stacked_face_params)
        e_tree.add_column(stack)

    e_tree.rename_leaves(genomes_descr.description.to_dict())
    e_tree.render(path, dpi=500)

    hsh_files = db.get_filenames_to_taxon_id()
    number_of_files = len(hsh_files)

    number_ort = db.get_n_orthogroups()
    taxids = list(genomes_descr.index)

    return render(request, 'chlamdb/home.html', my_locals(locals()))


def format_lst_to_html(lst_elem, add_count=True, format_func=lambda x: x):
    dict_elem = {}
    for elem in lst_elem:
        if pd.isna(elem):
            elem = "-"
        cnt = dict_elem.get(elem, 0)
        dict_elem[elem] = cnt + 1

    elems = []
    for k, v in dict_elem.items():
        if k != "-":
            token = format_func(k)
        else:
            token = k
        if add_count and k != "-":
            elems.append(f"{token} ({v})")
        else:
            elems.append(f"{token}")
    return "<br/>".join(elems)


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


def extract_orthogroup(request):
    biodb = settings.BIODB_DB_PATH
    db = DB.load_db(biodb, settings.BIODB_CONF)

    page_title = page2title["extract_orthogroup"]

    extract_form_class = make_extract_form(
        db, "extract_orthogroup", plasmid=True)
    if request.method != "POST":
        form = extract_form_class()
        return render(request, 'chlamdb/extract_orthogroup.html', my_locals(locals()))

    form = extract_form_class(request.POST)
    if not form.is_valid():
        form = extract_form_class()
        # add error message in web page
        return render(request, 'chlamdb/extract_orthogroup.html', my_locals(locals()))

    include_taxids, include_plasmids = form.get_include_choices()
    exclude_taxids, exclude_plasmids = form.get_exclude_choices()
    n_missing = form.get_n_missing()
    single_copy = "checkbox_single_copy" in request.POST

    sum_include_lengths = len(include_taxids)
    if include_plasmids is not None:
        sum_include_lengths += len(include_plasmids)

    if n_missing >= sum_include_lengths:
        wrong_n_missing = True
        return render(request, 'chlamdb/extract_orthogroup.html', my_locals(locals()))

    og_counts_in = db.get_og_count(include_taxids, plasmids=include_plasmids)
    if not single_copy:
        og_counts_in["presence"] = og_counts_in[og_counts_in > 0].count(axis=1)
        og_counts_in["selection"] = og_counts_in.presence >= (
            sum_include_lengths - n_missing)
    else:
        og_counts_in["presence"] = og_counts_in[og_counts_in == 1].count(
            axis=1)
        og_counts_in["absence"] = og_counts_in[og_counts_in == 0].count(axis=1)
        og_counts_in["selection"] = ((og_counts_in.presence >= (sum_include_lengths - n_missing))
                                     & (og_counts_in.absence + og_counts_in.presence == sum_include_lengths))

    sum_exclude_lengths = len(exclude_taxids)
    if exclude_plasmids is not None:
        sum_exclude_lengths += len(exclude_plasmids)
    if sum_exclude_lengths > 0:
        mat_exclude = db.get_og_count(
            exclude_taxids, plasmids=exclude_plasmids)
        mat_exclude["presence"] = mat_exclude[mat_exclude > 0].count(axis=1)
        mat_exclude["exclude"] = mat_exclude.presence > 0
        neg_index = mat_exclude[mat_exclude.exclude].index
    else:
        neg_index = pd.Index([])

    pos_index = og_counts_in[og_counts_in.selection].index
    selection = pos_index.difference(neg_index).tolist()
    if len(selection) == 0:
        no_match = True
        return render(request, 'chlamdb/extract_orthogroup.html', my_locals(locals()))

    count_all_genomes = db.get_og_count(selection, search_on="orthogroup")
    number_orth = count_all_genomes.index
    number_orth = len(number_orth)

    if not single_copy:
        orthogroup2count_all = count_all_genomes[count_all_genomes > 0].count(
            axis=1)
    else:
        orthogroup2count_all = count_all_genomes[count_all_genomes == 1].count(
            axis=1)
    max_n = orthogroup2count_all.max()
    match_groups_data = []

    all_taxids = include_taxids
    if include_plasmids is not None:
        all_taxids += include_plasmids
    annotations = db.get_genes_from_og(orthogroups=selection, taxon_ids=all_taxids,
                                       terms=["gene", "product", "locus_tag"])
    if annotations.empty:
        no_match = True
        return render(request, 'chlamdb/extract_orthogroup.html', my_locals(locals()))

    opt_header, optional_annotations = get_optional_annotations(
        db, seqids=annotations.index.tolist())
    details_header, details_data = get_table_details(db, annotations)
    annotations = annotations.join(optional_annotations)
    grouped = annotations.groupby("orthogroup")
    genes = grouped["gene"].apply(list)
    products = grouped["product"].apply(list)

    if "COG" in opt_header:
        cogs = grouped["cog"].apply(list)

    if "KO" in opt_header:
        kos = grouped["ko"].apply(list)

    table_headers = ["Orthogroup", "Genes", "Products"]
    table_headers.extend(opt_header)
    table_headers.extend(
        [f"Present in {sum_include_lengths}", f"Freq complete database (/{max_n})"])

    for row, count in orthogroup2count_all.items():
        cnt_in = og_counts_in.presence.loc[row]
        g = genes.loc[row]
        gene_data = format_lst_to_html(
            "-" if pd.isna(entry) else entry for entry in g)
        prod_data = format_lst_to_html(products.loc[row])
        column_header = format_orthogroup(row)
        optional = []
        if "KO" in opt_header and row in kos:
            optional.append(format_lst_to_html(
                kos.loc[row], add_count=True, format_func=format_ko_url))
        if "COG" in opt_header:
            optional.append(format_lst_to_html(
                cogs.loc[row], add_count=True, format_func=format_cog_url))
        entry = [column_header, gene_data, prod_data, *optional, cnt_in, count]
        match_groups_data.append(entry)

    ref_genomes = db.get_genomes_description(
    ).loc[include_taxids].reset_index()

    envoi_extract = True
    return render(request, 'chlamdb/extract_orthogroup.html', my_locals(locals()))


def venn_orthogroup(request):
    biodb_path = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb_path)
    page_title = page2title["venn_orthogroup"]

    venn_form_class = make_venn_from(db, limit=6)
    if request.method != "POST":
        form_venn = venn_form_class()
        return render(request, 'chlamdb/venn_orthogroup.html', my_locals(locals()))

    form_venn = venn_form_class(request.POST)
    if not form_venn.is_valid():
        return render(request, 'chlamdb/venn_orthogroup.html', my_locals(locals()))

    targets = form_venn.get_taxids()
    genomes = db.get_genomes_description()
    og_count = db.get_og_count(targets)
    fmt_data = []
    for taxon in og_count:
        ogs = og_count[taxon]
        ogs_str = ",".join(f"{to_s(format_orthogroup(og))}"
                           for og, cnt in ogs.items() if cnt > 0)
        genome = genomes.loc[int(taxon)].description
        fmt_data.append(f"{{name: {to_s(genome)}, data: [{ogs_str}]}}")
    series = "[" + ",".join(fmt_data) + "]"

    og_list = og_count.index.tolist()
    annotations = db.get_genes_from_og(
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
    envoi_venn = True
    return render(request, 'chlamdb/venn_orthogroup.html', my_locals(locals()))


def format_pfam(pfam_id, base=None, to_url=False):
    if base is None:
        fmt_entry = f"PF{pfam_id:04d}"
    else:
        fmt_entry = base
    if to_url:
        return f"<a href=/fam_pfam/{fmt_entry}>{fmt_entry}</a>"
    return fmt_entry


def index_comp(request, type):
    page_title = page2title[f"index_comp_{type}"]
    if type == 'Pfam':
        return render(request, 'chlamdb/index_pfam.html', my_locals(locals()))
    if type == 'COG':
        return render(request, 'chlamdb/index_cog.html', my_locals(locals()))
    if type == 'ko':
        return render(request, 'chlamdb/index_ko.html', my_locals(locals()))
    if type == 'orthology':
        return render(request, 'chlamdb/index_orthology.html', my_locals(locals()))


def entry_list_ko(request,):
    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    page_title = page2title["entry_list_ko"]

    # retrieve taxid list
    genomes_data = db.get_genomes_infos()
    taxids = [str(i) for i in genomes_data.index.to_list()]

    # retrieve entry list
    ko_all = db.get_ko_hits(taxids,
                            search_on="taxid",
                            indexing="taxid")
    # retrieve annotations
    ko_desc = db.get_ko_desc(ko_all.index.to_list())
    ko_mod = db.get_ko_modules(ko_all.index.to_list())
    ko_path = db.get_ko_pathways(ko_all.index.to_list())

    # count frequency and n genomes
    combined_df = pd.DataFrame(ko_all.sum(axis=1).rename('count'))
    ko_freq = ko_all[ko_all > 0].count(axis=1).to_dict()

    combined_df["accession"] = [format_ko(ko) for ko in combined_df.index]
    combined_df["modules"] = [format_ko_modules(ko_mod, ko) if ko in ko_mod else '-' for ko in combined_df.index]
    combined_df["description"] = [ko_desc[ko] for ko in combined_df.index]
    combined_df["pathways"] = [format_ko_path(ko_path, ko) if ko in ko_path else '-' for ko in combined_df.index]
    combined_df["freq"] = [ko_freq[ko] for ko in combined_df.index]

    combined_df = combined_df.sort_values(["count", "freq"], ascending=False)

    # combine into df
    return render(request, 'chlamdb/entry_list_ko.html', my_locals(locals()))


def entry_list_cog(request,):
    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    page_title = page2title["entry_list_cog"]
    # retrieve taxid list
    genomes_data = db.get_genomes_infos()
    taxids = [str(i) for i in genomes_data.index.to_list()]

    # retrieve entry list
    cog_all = db.get_cog_hits(taxids,
                              search_on="taxid",
                              indexing="taxid")
    # retrieve annotations
    cogs_summaries = db.get_cog_summaries(cog_all.index.tolist(), as_df=True, only_cog_desc=True)

    # count frequency and n genomes
    cog_count = cog_all.sum(axis=1)
    cog_freq = cog_all[cog_all > 0].count(axis=1)
    cogs_summaries["accession"] = [format_cog(cog) for cog in cogs_summaries.index]

    # combine into df
    combined_df = cogs_summaries.merge(cog_count.rename('count'),
                                       left_index=True,
                                       right_index=True).merge(cog_freq.rename('freq'),
                                                               left_index=True,
                                                               right_index=True).sort_values(["count"], ascending=False)

    return render(request, 'chlamdb/entry_list_cog.html', my_locals(locals()))


def entry_list_pfam(request,):
    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    page_title = page2title["entry_list_pfam"]
    # retrieve taxid list
    genomes_data = db.get_genomes_infos()
    taxids = [str(i) for i in genomes_data.index.to_list()]

    # retrieve entry list
    pfam_all = db.get_pfam_hits(taxids,
                                search_on="taxid",
                                indexing="taxid")
    # retrieve annotations
    pfam_annot = db.get_pfam_def(pfam_all.index.to_list())

    # count frequency and n genomes
    pfam_count = pfam_all.sum(axis=1)
    pfam_freq = pfam_all[pfam_all > 0].count(axis=1)
    pfam_annot["accession"] = [format_pfam(pfam) for pfam in pfam_annot.index]

    # combine into df
    combined_df = pfam_annot.merge(pfam_count.rename('count'),
                                   left_index=True,
                                   right_index=True)\
                            .merge(pfam_freq.rename('freq'),
                                   left_index=True,
                                   right_index=True)\
                            .sort_values(["count"],
                                         ascending=False)

    return render(request, 'chlamdb/entry_list_pfam.html', my_locals(locals()))


def extract_pfam(request, classification="taxon_id"):
    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    page_title = page2title["extract_pfam"]

    extract_form_class = make_extract_form(db, "extract_pfam", plasmid=True, label="Pfam domains")
    if request.method != "POST":
        form = extract_form_class()
        return render(request, 'chlamdb/extract_Pfam.html', my_locals({"form": form, "page_title": page_title}))

    form = extract_form_class(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/extract_Pfam.html', my_locals(locals()))

    include, include_plasmids = form.get_include_choices()
    exclude, exclude_plasmids = form.get_exclude_choices()
    n_missing = form.get_n_missing()

    sum_include_length = len(include)
    if include_plasmids is not None:
        sum_include_length += len(include_plasmids)

    sum_exclude_length = len(exclude)
    if exclude_plasmids is not None:
        sum_exclude_length += len(exclude_plasmids)

    if n_missing >= sum_include_length:
        ctx = {"wrong_n_missing": True, "form": form, "page_title": page_title}
        return render(request, 'chlamdb/extract_Pfam.html', my_locals(ctx))

    pfam_include = db.get_pfam_hits(include, plasmids=include_plasmids,
                                    search_on="taxid", indexing="taxid")
    if sum_exclude_length > 0:
        pfam_exclude = db.get_pfam_hits(exclude, plasmids=exclude_plasmids,
                                        search_on="taxid", indexing="taxid")
        pfam_exclude["sum_pos"] = pfam_exclude[pfam_exclude > 0].count(axis=1)
        pfam_exclude["exclude"] = pfam_exclude.sum_pos > 0
        neg_index = pfam_exclude[pfam_exclude.exclude].index
    else:
        neg_index = pd.Index([])

    pfam_include["sum_pos"] = pfam_include[pfam_include > 0].count(axis=1)
    pfam_include["selection"] = pfam_include.sum_pos >= len(include) - n_missing
    pos_index = pfam_include[pfam_include.selection].index
    selection = pos_index.difference(neg_index).tolist()
    pfam_defs = db.get_pfam_def(selection)

    if len(selection) == 0:
        ctx = {"no_match": True, "form": form, "page_title": page_title}
        return render(request, 'chlamdb/extract_Pfam.html', my_locals(ctx))

    all_database = db.get_pfam_hits(pfam_include.index.tolist(), search_on="pfam", indexing="taxid")
    sums = all_database.sum(axis=1)
    sum_group = len(selection)

    match_groups_data = []
    for no, pfam in enumerate(selection):
        count = sums.loc[pfam]
        pfam_def = pfam_defs["def"].loc[pfam]
        data = [no + 1, format_pfam(pfam), pfam_def,
                pfam_include.sum_pos.loc[pfam], sums.loc[pfam]]
        match_groups_data.append(data)

    ctx = {"envoi_extract": True,
           "sum_group": sum_group,
           "n_genomes": sum_include_length,
           "max_n": sums.max(),
           "match_groups_data": match_groups_data,
           "form": form,
           "sum_include_length": sum_include_length,
           "sum_exclude_length": sum_exclude_length,
           "n_missing": n_missing,
           "page_title": page_title}
    return render(request, 'chlamdb/extract_Pfam.html', my_locals(ctx))


def format_ko(ko_id, as_url=False, base=None):
    if base is None:
        base = f"K{int(ko_id):05d}"
    if not as_url:
        return base
    return f"<a href=\"/fam_ko/{base}\">{base}</a>"


def format_ko_url(ko_id):
    return format_ko(ko_id, as_url=True)


def format_ko_path(hsh_pathways, ko, as_list=False, with_taxid=None):
    pathways = hsh_pathways.get(ko, [])
    if len(pathways) == 0:
        if as_list:
            return []
        return "-"
    if with_taxid is None:
        fmt_lst = (f"<a href=\"/KEGG_mapp_ko/map{i:05d}\">{d}</a>" for i, d in pathways)
    else:
        fmt_lst = (f"<a href=\"/KEGG_mapp_ko/map{i:05d}/{with_taxid}\">{d}</a>" for i, d in pathways)

    if as_list:
        return list(fmt_lst)
    return "<br>".join(fmt_lst)


def format_ko_module(module_id, module_desc=None):
    if module_desc is None:
        return f"<a href=\"/KEGG_module_map/M{module_id:05d}\">M{module_id:05d}</a>"
    else:
        return f"<a href=\"/KEGG_module_map/M{module_id:05d}\">{module_desc}</a>"


def format_ko_modules(hsh_modules, ko):
    modules = hsh_modules.get(ko, [])
    if len(modules) == 0:
        return "-"
    return "<br>".join([format_ko_module(i, d) for i, d in modules])


def extract_ko(request):
    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    page_title = page2title["extract_ko"]

    extract_form_class = make_extract_form(db, "extract_ko", plasmid=True, label="Kegg Orthologs")

    if request.method != "POST":
        form = extract_form_class()
        return render(request, 'chlamdb/extract_ko.html', my_locals(locals()))

    form = extract_form_class(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/extract_ko.html', my_locals(locals()))

    include, include_plasmids = form.get_include_choices()
    exclude, exclude_plasmids = form.get_exclude_choices()
    n_missing = form.get_n_missing()

    sum_include_length = len(include)
    if include_plasmids is not None:
        sum_include_length += len(include_plasmids)

    sum_exclude_length = len(exclude)
    if exclude_plasmids is not None:
        sum_exclude_length += len(exclude_plasmids)

    if n_missing >= sum_include_length:
        hsh_var = {"wrong_n_missing": True, "form": form}
        return render(request, 'chlamdb/extract_ko.html', my_locals(hsh_var))

    mat_include = db.get_ko_hits(include, plasmids=include_plasmids)
    if len(exclude) > 0:
        mat_exclude = db.get_ko_hits(exclude, plasmids=exclude_plasmids)
        mat_exclude["sum_pos"] = mat_exclude[mat_exclude > 0].count(axis=1)
        mat_exclude["exclude"] = mat_exclude.sum_pos > 0
        neg_index = mat_exclude[mat_exclude.exclude].index
    else:
        neg_index = pd.Index([])

    mat_include["sum_pos"] = mat_include[mat_include > 0].count(axis=1)
    mat_include["selection"] = mat_include.sum_pos >= len(include) - n_missing
    pos_index = mat_include[mat_include.selection].index
    selection = pos_index.difference(neg_index).tolist()
    if len(selection) == 0:
        hsh_var = {"no_match": True, "form": form}
        return render(request, 'chlamdb/extract_ko.html', my_locals(hsh_var))

    all_database = db.get_ko_hits(selection, search_on="ko", indexing="taxid")
    ko_total_count = all_database.sum(axis=1)
    ko_desc = db.get_ko_desc(selection)
    ko_mod = db.get_ko_modules(selection)
    ko_path = db.get_ko_pathways(selection)
    match_groups_data = []
    for ko in selection:
        kof = format_ko(ko)
        kod = ko_desc.get(ko, "-")
        kop = format_ko_path(ko_path, ko)
        kom = format_ko_modules(ko_mod, ko)
        kot = ko_total_count.loc[ko]
        data = [kof, kod, kop, kom, mat_include.sum_pos.loc[ko], kot]
        match_groups_data.append(data)

    max_n = ko_total_count.max()
    envoi_extract = True
    return render(request, 'chlamdb/extract_ko.html', my_locals(locals()))


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


def extract_cog(request):
    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    extract_form_class = make_extract_form(db, "extract_cog", plasmid=True, label="COG")
    page_title = page2title["extract_cog"]

    if request.method != "POST":
        form = extract_form_class()
        return render(request, 'chlamdb/extract_cogs.html', my_locals(locals()))

    form = extract_form_class(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/extract_cogs.html', my_locals(locals()))

    include, include_plasmids = form.get_include_choices()
    exclude, exclude_plasmids = form.get_exclude_choices()
    n_missing = form.get_n_missing()

    sum_include_length = len(include)
    if include_plasmids is not None:
        sum_include_length += len(include_plasmids)

    sum_exclude_length = len(exclude)
    if exclude_plasmids is not None:
        sum_exclude_length += len(exclude_plasmids)

    if n_missing >= sum_include_length:
        wrong_n_missing = True
        return render(request, 'chlamdb/extract_cogs.html', my_locals(locals()))

    cog_include = db.get_cog_hits(include, plasmids=include_plasmids,
                                  search_on="taxid", indexing="taxid")
    if sum_exclude_length > 0:
        cog_exclude = db.get_cog_hits(exclude, plasmids=exclude_plasmids,
                                      search_on="taxid", indexing="taxid")
        cog_exclude["sum_pos"] = cog_exclude[cog_exclude > 0].count(axis=1)
        cog_exclude["exclude"] = cog_exclude.sum_pos > 0
        neg_index = cog_exclude[cog_exclude.exclude].index
    else:
        neg_index = pd.Index([])

    cog_include["sum_pos"] = cog_include[cog_include > 0].count(axis=1)
    cog_include["selection"] = cog_include.sum_pos >= len(include) - n_missing
    pos_index = cog_include[cog_include.selection].index
    selection = pos_index.difference(neg_index).tolist()
    if len(selection) == 0:
        no_match = True
        return render(request, 'chlamdb/extract_cogs.html', my_locals(locals()))

    all_database = db.get_cog_hits(cog_include.index.tolist(), search_on="cog", indexing="taxid")
    sums = all_database.sum(axis=1)

    cat_count = {}
    cogs_summaries = db.get_cog_summaries(sums.index.tolist())
    cogs_funct = db.get_cog_code_description()
    cog_data = []
    for cog_id in selection:
        count = sums.loc[cog_id]

        # some cogs do not have a description, skip those
        if cog_id not in cogs_summaries:
            continue

        data = [format_cog(cog_id)]
        func_acc = []
        for func, func_descr, cog_descr in cogs_summaries[cog_id]:
            func_acc.append((func, func_descr))
            inc, not_incl = cat_count.get(func, (0, 0))
            cat_count[func] = (inc + cog_include.sum_pos.loc[cog_id], not_incl)
        funcs = "<br>".join(f"{func} ({func_desc})" for func, func_desc in func_acc)
        data = (format_cog(cog_id), funcs, cog_descr,
                cog_include.sum_pos.loc[cog_id], str(count))
        cog_data.append(data)

    # get the categories for all cogs
    for cog_id, details_lst in cogs_summaries.items():
        for func, func_descr, cog_descr in details_lst:
            inc, not_incl = cat_count.get(func, (0, 0))
            cat_count[func] = (inc, not_incl + sums.loc[cog_id])

    max_n = sums.max(axis=0)
    sum_group = len(selection)

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
    envoi_extract = True
    return render(request, 'chlamdb/extract_cogs.html', my_locals(locals()))


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


def format_cog(cog_id, as_url=False, base=None):
    if base is None:
        base = f"COG{int(cog_id):04d}"
    if as_url is False:
        return base
    return f"<a href=\"/fam_cog/{base}\">{base}</a>"


def format_cog_url(cog_id):
    return format_cog(cog_id, as_url=True)


def escape_quotes(unsafe):
    return unsafe.replace("\"", "\\\"")


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
    cog_hits = db.get_cog_hits(targets, indexing="taxid", search_on="taxid")
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


def genomes(request):
    biodb_path = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb_path)
    page_title = page2title["genomes"]
    genomes_data = db.get_genomes_infos()
    genomes_descr = db.get_genomes_description()

    genomes_data = genomes_data.join(genomes_descr)

    genomes_data.gc = genomes_data.gc.apply(round)
    genomes_data.coding_density = genomes_data.coding_density.apply(lambda x: round(100 * x))
    genomes_data.length = genomes_data.length.apply(lambda x: round(x / pow(10, 6), 2))

    filenames_tax_id = db.get_filenames_to_taxon_id()
    filenames_tax_id_db = pd.DataFrame.from_dict(list(filenames_tax_id.items()))
    filenames_tax_id_db.columns = ['filename', 'taxon_id']
    filenames_tax_id_db.index = list(filenames_tax_id_db['taxon_id'])
    filenames_list = list(filenames_tax_id_db["filename"])

    path_faa = [settings.BLAST_DB_PATH + "/faa/" + filename + ".faa" for filename in filenames_list]
    path_fna = [settings.BLAST_DB_PATH + "/fna/" + filename + ".fna" for filename in filenames_list]
    path_ffn = [settings.BLAST_DB_PATH + "/ffn/" + filename + ".ffn" for filename in filenames_list]
    path_gbk = [settings.BLAST_DB_PATH + "/gbk/" + filename + ".gbk" for filename in filenames_list]

    filenames_tax_id_db['path_to_faa'] = path_faa
    filenames_tax_id_db['path_to_fna'] = path_fna
    filenames_tax_id_db['path_to_ffn'] = path_ffn
    filenames_tax_id_db['path_to_gbk'] = path_gbk
    filenames_tax_id_db = filenames_tax_id_db[["path_to_faa", "path_to_fna", "path_to_ffn", "path_to_gbk"]]
    genomes_data = genomes_data.join(filenames_tax_id_db, on="taxon_id")
    data_table_header = [
        "Name",
        "%GC",
        "N proteins",
        "N contigs",
        "Size (Mbp)",
        "Percent coding",
        "N plasmid contigs",
        "faa seq",
        "fna seq",
        "ffn seq",
        "gbk file"
    ]
    data_table = genomes_data[[
        "id",
        "description",
        "gc",
        "n_prot",
        "n_contigs",
        "length",
        "coding_density",
        "has_plasmid",
        "path_to_faa",
        "path_to_fna",
        "path_to_ffn",
        "path_to_gbk"]].values.tolist()
    return render(request, 'chlamdb/genomes.html', my_locals(locals()))


def extract_contigs(request, genome):
    page_title = page2title["extract_contigs"]

    taxid = int(genome)
    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)

    descr = db.get_genomes_description().description.to_dict()
    prot_infos = db.get_proteins_info([taxid], search_on="taxid",
                                      as_df=True, to_return=["locus_tag", "product", "gene"])
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


def format_lst(lst):
    hsh_values = {}
    for item in lst:
        val = hsh_values.get(item, 0)
        hsh_values[item] = val + 1
    return hsh_values


def tab_homologs(db, infos, hsh_organism, ref_seqid=None, og=None):
    if ref_seqid is None:
        orthogroup_title = f"Homologs in group_{og}"
    else:
        locus_tag = infos.loc[ref_seqid].locus_tag
        orthogroup_title = f"Homologs of {locus_tag}"

    headers = ["", "Locus tag", "Source", "Gene", "Product"]
    identities = None
    if ref_seqid is not None:
        identities = db.get_og_identity(og=og, ref_seqid=ref_seqid)
        headers.insert(2, "Identity")

    homologues = []
    index = 0
    orga_set = set()
    for seqid, data in infos.iterrows():
        organism = hsh_organism[seqid]
        locus_fmt = format_locus(data.locus_tag, to_url=True)
        entry = [index + 1, locus_fmt, organism, format_gene(data.gene), data["product"]]
        if ref_seqid is not None:
            if seqid == ref_seqid:
                continue
            else:
                orga_set.add(organism)
                ident = round(identities.loc[seqid].identity, 1)
            if ident == 0:
                ident = "-"
            entry.insert(2, ident)
        homologues.append(entry)
        index += 1

    n_genomes = len(orga_set) if ref_seqid is not None else len(set(hsh_organism.values()))
    return {"orthogroup": orthogroup_title,
            "n_genomes": "1 genome" if n_genomes == 1 else f"{n_genomes} genomes",
            "headers": headers,
            "homologues": homologues}


def make_div(figure_or_data, include_plotlyjs=False, show_link=False, div_id=None):
    from plotly import offline
    div = offline.plot(
        figure_or_data,
        include_plotlyjs=include_plotlyjs,
        show_link=show_link,
        output_type="div",
    )
    if ".then(function ()" in div:
        div = """{div.partition(".then(function ()")[0]}</script>"""
    if div_id:
        import re

        try:
            existing_id = re.findall(r'id="(.*?)"|$', div)[0]
            div = div.replace(existing_id, div_id)
        except IndexError:
            pass
    return div


def tab_lengths(n_homologues, annotations):
    import plotly.figure_factory as ff

    length_distrib = n_homologues > 1
    if not length_distrib:
        return {"length_distrib": False}

    lengths = annotations["length"]
    max_protein_length = lengths.max()
    std_protein_length = f"{lengths.std():.1f}"
    min_protein_length = lengths.min()
    mean_protein_length = f"{lengths.mean():.1f}"
    median_protein_length = f"{lengths.median():.1f}"
    if len(lengths.unique()) > 1:
        fig1 = ff.create_distplot([lengths.tolist()], ["Sequence length"], bin_size=20)
        fig1.update_xaxes(range=[0, max_protein_length])
        fig1.layout.margin.update({"l": 80, "r": 20, "b": 40, "t": 20, "pad": 10, })
        html_plot_prot_length = make_div(fig1, div_id="distplot")
    else:
        return {"length_distrib": True, "single_length": True, "prot_length": lengths.iloc[0]}

    return {"length_distrib": True,
            "max_protein_length": max_protein_length,
            "std_protein_length": std_protein_length,
            "min_protein_length": min_protein_length,
            "mean_protein_length": mean_protein_length,
            "median_protein_length": median_protein_length,
            "html_plot_prot_length": html_plot_prot_length}


class SimpleTextColumn(Column):
    def __init__(self, header=None):
        super().__init__(header)

    def get_face(self, index):
        return TextFace(index, fsize=7)


class PfamColumn(Column):
    def __init__(self, header, pfam_col, pfam_cmap):
        super().__init__(header)
        self.pfam_col = pfam_col
        self.pfam_cmap = pfam_cmap

    def get_face(self, index):
        prot_length, pfam_infos = self.pfam_col[index]
        dummy_seq = "-" * prot_length
        pfam_entries = []
        for pfam, start, end in pfam_infos:
            fmt_entry = f"arial|6|white|{format_pfam(pfam)}"
            entry = [start, end, "[]", None, 8, "black", self.pfam_cmap[pfam], fmt_entry]
            pfam_entries.append(entry)
        return SeqMotifFace(dummy_seq, motifs=pfam_entries, seq_format="line")


def og_tab_get_swissprot_homologs(db, annotations):
    homologs = db.get_swissprot_homologs(annotations.index.tolist(), indexing="accession")
    summary = homologs.groupby("definition").count()
    swissprot = []
    for definition, data in summary.iterrows():
        swissprot.append([definition, data.accession])
    return {"reviewed": swissprot, "n_swissprot_hits": len(homologs)}


def prepare_default_tree(og_phylogeny):
    tree = Tree(og_phylogeny)
    R = tree.get_midpoint_outgroup()
    if R is not None:
        root = "(midpoint rooted)"
        tree.set_outgroup(R)
    tree.ladderize()

    return tree, root


def tab_og_phylogeny(db, og_id, compare_to=None):
    og_phylogeny = db.get_og_phylogeny(og_id)
    pfam_col = None
    ident_col = None
    if optional2status.get("pfam", False):
        annots = db.get_genes_from_og(orthogroups=[og_id], terms=["locus_tag", "length"])
        pfams = db.get_pfam_hits_info(annots.index.tolist())
        unique_pfams = pfams.pfam.unique()
        color_palette = (mpl_col.to_hex(col) for col in sns.color_palette(None, len(unique_pfams)))
        pfam_cmap = dict(zip(unique_pfams, color_palette))
        tmp_hsh_infos = collections.defaultdict(list)
        hsh_pfam_infos = {}
        for index, infos in pfams.iterrows():
            tmp_hsh_infos[infos.seqid].append([infos.pfam, infos.start, infos.end])
        for seqid, data in annots.iterrows():
            pfam_entries = tmp_hsh_infos.get(seqid, [])
            hsh_pfam_infos[data.locus_tag] = [data.length, pfam_entries]
        pfam_col = PfamColumn("Pfam domains", hsh_pfam_infos, pfam_cmap)

    if compare_to is not None:
        identity_matrix = db.get_og_identity(og=og_id, ref_seqid=compare_to)
        seqids = identity_matrix.index.tolist()
        seqids.append(compare_to)
        seqid_to_locus = db.get_proteins_info(seqids, to_return=["locus_tag"], as_df=True)
        all_infos = identity_matrix.join(seqid_to_locus).set_index("locus_tag").round(1)
        all_infos.loc[seqid_to_locus.loc[compare_to].locus_tag] = 100.0

        ident_col = SimpleColorColumn.fromSeries(all_infos.identity, color_gradient=True,
                                                 header="Identity", default_val="-", is_str_index=True)

    tree, root = prepare_default_tree(og_phylogeny)
    locuses = [branch.name for branch in tree.iter_leaves()]
    locus_to_genome = db.get_locus_to_genomes(locuses)

    og_filename = f"OG{og_id: 07}_mafft.faa"

    e_tree = EteTree(tree)
    e_tree.add_column(SimpleTextColumn("Locus tag"))
    if ident_col is not None:
        e_tree.add_column(ident_col)
    if pfam_col is not None:
        e_tree.add_column(pfam_col)
    e_tree.rename_leaves(locus_to_genome, leaf_name_type=str)

    asset_path = f"/temp/og_phylogeny{og_id}.svg"
    path = settings.BASE_DIR + '/assets/' + asset_path
    e_tree.render(path, dpi=1200)

    algn_file = f"/alignments/{og_filename}"
    return {"og_phylogeny": asset_path, "root": root, "og_alignment": algn_file}


def tab_og_conservation_tree(db, group, compare_to=None):
    ref_phylogeny = db.get_reference_phylogeny()
    leaf_to_name = db.get_genomes_description().description.to_dict()

    # Note: as the orthogroup may either be in a plasmid or in the chromosome
    # of the bacteria, we need to index by taxon to group them (index on taxon)
    count = db.get_og_count([group], search_on="orthogroup")

    tree = Tree(ref_phylogeny)
    R = tree.get_midpoint_outgroup()
    if R is not None:
        tree.set_outgroup(R)
    tree.ladderize()
    e_tree = EteTree(tree)

    e_tree.add_column(SimpleColorColumn.fromSeries(count.loc[group], header="Number of homologs"))
    if compare_to is not None:
        identity_matrix = db.get_og_identity(og=group, ref_seqid=compare_to)
        seqids = identity_matrix.index.tolist()

        # to get the taxid of the reference seqid, so as to exclude it from
        # the phylogenetic tree
        seqids.append(compare_to)
        seqid_to_taxon = db.get_taxid_from_seqid(seqids)
        identity_matrix["taxid"] = identity_matrix.index.map(seqid_to_taxon)
        max_identity = identity_matrix.groupby("taxid").max().round(1)
        max_identity.loc[seqid_to_taxon[compare_to]] = 100.0
        col = SimpleColorColumn.fromSeries(max_identity.identity, color_gradient=True,
                                           header="Identity", default_val="-")
        e_tree.add_column(col)

    e_tree.rename_leaves(leaf_to_name)

    dpi = 1200
    asset_path = f"/temp/og_conservation{group}.svg"
    path = settings.BASE_DIR + '/assets/' + asset_path
    e_tree.render(path, dpi=dpi)
    return {"asset_path": asset_path}


def og_tab_get_kegg_annot(db, seqids, from_taxid=None):
    ko_hits = db.get_ko_hits(seqids, search_on="seqid", indexing="seqid")
    if ko_hits.empty:
        return {}

    n_occurences = ko_hits["ko"].value_counts()
    ko_ids = n_occurences.index.tolist()
    ko_descr = db.get_ko_desc(ko_ids)
    ko_pathways = db.get_ko_pathways(ko_ids)
    ko_modules = db.get_ko_modules(ko_ids)

    ko_entries = []
    for ko_id, count in n_occurences.items():
        entry = [format_ko(ko_id, as_url=True), count]
        descr = ko_descr.get(ko_id, "-")
        entry.append(descr)
        entry.append(format_ko_path(ko_pathways, ko_id, with_taxid=from_taxid))
        entry.append(format_ko_modules(ko_modules, ko_id))
        ko_entries.append(entry)

    ko_header = ["KO", "Occurences", "Description", "Pathways", "Modules"]
    return {
        "ko_entries": ko_entries,
        "ko_header": ko_header
    }


def og_tab_get_amr_annot(db, seqids):
    amr_hits = db.get_amr_hits_from_seqids(seqids)
    if amr_hits.empty:
        return {}

    col_titles = {"closest_seq": "Closest Sequence"}

    amr_hits["closest_seq"] = amr_hits["closest_seq"].map(format_refseqid_to_ncbi)
    amr_hits["gene"] = amr_hits[["gene", "hmm_id"]].apply(format_gene_to_ncbi_hmm, axis=1)
    return {
        "amr_entries": amr_hits.values,
        "amr_header": [col_titles.get(col, col.capitalize())
                       for col in amr_hits.columns]
    }


def og_tab_get_cog_annot(db, seqids):
    cog_hits = db.get_cog_hits(seqids, indexing="seqid", search_on="seqid")

    if cog_hits.empty:
        return {}

    n_entries = cog_hits["cog"].value_counts()
    cog_summ = db.get_cog_summaries(n_entries.index.tolist())
    cog_entries = []
    for cog_id, count in n_entries.items():
        if cog_id not in cog_summ:
            # should add a warning on the web page
            continue

        entry = [format_cog(cog_id, as_url=True), count]
        funcs = []
        func_descrs = []
        cog_descrs = []
        for func, func_descr, cog_descr in cog_summ[cog_id]:
            funcs.append(func)
            func_descrs.append(func_descr)
            cog_descrs.append(cog_descr)
        entry.append(cog_descrs.pop())
        entry.append("<br>".join(funcs))
        entry.append("<br>".join(func_descrs))
        cog_entries.append(entry)

    if len(cog_entries) == 0:
        return {}

    cog_header = ["COG", "Occurences", "Description", "Category", "Category description"]
    return {
        "cog_header": cog_header,
        "cog_entries": cog_entries
    }


def og_tab_get_pfams(db, annotations):
    seqids = annotations.index.tolist()

    pfam_hits = db.get_pfam_hits(seqids, search_on="seqid", indexing="seqid")
    pfam_hits_count = pfam_hits.groupby(["pfam"]).count()
    pfam_defs = db.get_pfam_def(pfam_hits_count.index.tolist())
    all_infos = pfam_hits_count.join(pfam_defs)

    all_infos_tab = []
    for pfam, data in all_infos.iterrows():
        all_infos_tab.append([format_pfam(pfam, to_url=True), data.seqid, data["def"]])
    return {
        "pfam_def": all_infos_tab
    }


def tab_og_best_hits(db, orthogroup, locus=None):
    try:
        refseq_newick = db.get_refseq_phylogeny(orthogroup)
    except Exception:
        # no phylogeny for that orthogroup
        return {"has_refseq_phylo": False}
    ete_tree = Tree(refseq_newick)
    loci = list(leaf.name.split(".")[0] for leaf in ete_tree.iter_leaves())
    match_infos = db.get_refseq_matches_info(loci, search_on="accession")
    zdb_taxids = db.get_taxid_from_accession(loci)
    orgas = db.get_genomes_description().description.to_dict()
    acc_to_orga = match_infos.set_index("accession")["organism"]

    R = ete_tree.get_midpoint_outgroup()
    if R is not None:
        ete_tree.set_outgroup(R)
    ete_tree.ladderize()

    for leaf in ete_tree.iter_leaves():
        shortened = leaf.name.split(".")[0]
        if shortened in acc_to_orga.index:
            orga_name = acc_to_orga.loc[shortened]
            leaf.add_face(TextFace(f"{leaf.name} | {orga_name}"), 0, "branch-right")
            continue

        color = "red"
        if locus is not None and shortened == locus:
            color = "green"
        taxid = zdb_taxids.loc[shortened].taxid
        orga_name = orgas[taxid]
        leaf.add_face(TextFace(f"{leaf.name} | {orga_name}", fgcolor=color),
                      0, "branch-right")

    asset_path = f"/temp/og_best_hit_phylogeny_{orthogroup}.svg"
    path = settings.BASE_DIR + '/assets/' + asset_path
    ts = TreeStyle()
    ts.show_leaf_name = False
    ete_tree.render(path, tree_style=ts, dpi=1200)
    return {"best_hits_phylogeny": asset_path, "has_refseq_phylo": True}


def orthogroup(request, og):
    tokens = og.split("_")
    try:
        og_id = int(tokens[1])
    except Exception:
        menu = True
        invalid_id = True
        return render(request, "chlamdb/og.html", my_locals(locals()))
    else:
        valid_id = True

    biodb = settings.BIODB_DB_PATH
    db = DB.load_db(biodb, settings.BIODB_CONF)

    og_counts = db.get_og_count([og_id], search_on="orthogroup")

    if len(og_counts.index) == 0:
        valid_id = False
        return render(request, "chlamdb/og.html", my_locals(locals()))

    annotations = db.get_genes_from_og(orthogroups=[og_id],
                                       terms=["locus_tag", "gene", "product", "length"])

    hsh_organisms = db.get_organism(annotations.index.tolist())
    hsh_genes = format_lst(annotations["gene"].tolist())
    hsh_products = format_lst(annotations["product"].tolist())
    n_homologues = og_counts.loc[og_id].sum()

    gene_annotations = []
    for index, values in enumerate(hsh_genes.items()):
        gene, cnt = values
        if pd.isna(gene):
            gene = "-"
        gene_annotations.append([index + 1, gene, cnt])

    product_annotations = []
    for index, values in enumerate(hsh_products.items()):
        product, cnt = values
        if pd.isna(product):
            product = "-"
        product_annotations.append([index + 1, product, cnt])

    swissprot, cog_ctx, kegg_ctx, pfam_ctx, amr_ctx = {}, {}, {}, {}, {}
    best_hit_phylo = {}
    if optional2status.get("COG", False):
        cog_ctx = og_tab_get_cog_annot(db, annotations.index.tolist())

    if optional2status.get("KEGG", False):
        kegg_ctx = og_tab_get_kegg_annot(db, annotations.index.tolist())

    if optional2status.get("pfam", False):
        pfam_ctx = og_tab_get_pfams(db, annotations)

    if optional2status.get("BLAST_swissprot", False):
        swissprot = og_tab_get_swissprot_homologs(db, annotations)

    try:
        og_phylogeny_ctx = tab_og_phylogeny(db, og_id)
    except NoPhylogenyException:
        og_phylogeny_ctx = {}

    if optional2status.get("BBH_phylogenies", False):
        best_hit_phylo = tab_og_best_hits(db, og_id)

    if optional2status.get("AMR", False):
        amr_ctx = og_tab_get_amr_annot(db, annotations.index.tolist())

    og_conserv_ctx = tab_og_conservation_tree(db, og_id)
    length_tab_ctx = tab_lengths(n_homologues, annotations)
    homolog_tab_ctx = tab_homologs(db, annotations, hsh_organisms, og=og_id)
    context = {
        "valid_id": valid_id,
        "n_homologues": n_homologues,
        "og": og,
        "menu": True,
        "gene_annotations": gene_annotations,
        "product_annotations": product_annotations,
        **homolog_tab_ctx,
        **length_tab_ctx,
        **og_conserv_ctx, **best_hit_phylo,
        **cog_ctx, **kegg_ctx, **pfam_ctx, **og_phylogeny_ctx, **swissprot,
        **amr_ctx
    }
    return render(request, "chlamdb/og.html", my_locals(context))

def tab_general(db, seqid):
    hsh_infos = db.get_proteins_info([seqid], to_return=["locus_tag", "gene", "product"],
            inc_non_CDS=True, inc_pseudo=True)
    hsh_organism = db.get_organism([seqid])
    gene_loc = db.get_gene_loc([seqid], as_hash=False)

    gene_pos = []
    nucl_length = 0
    for index, row in gene_loc.iterrows():
        gene_pos.append((row.start, row.end, row.strand))
        nucl_length += row.end-row.start+1

    locus_tag, gene, product = hsh_infos[seqid]
    if pd.isna(gene):
        gene = "-"
    organism = hsh_organism[seqid]
    return {
        "locus_tag": locus_tag,
        "organism": organism,
        "gene_pos" : gene_pos,
        "gene": gene,
        "nucl_length": nucl_length,
        "prot": product
    }

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

        color = colors.linearlyInterpolatedColor(colors.gray,
                                                 colors.firebrick, self.min_val, self.max_val, val)
        text_face = TextFace(str(int(val)).center(12 - len(str(int(val)))))
        text_face.inner_background.color = to_color_code(color)
        self.set_default_params(text_face)
        return text_face


def get_sequence(db, seqid, flanking=0):
    loc = db.get_gene_loc([seqid], as_hash=False)
    bioentry, accession, length, seq = db.get_bioentry_list(seqid, search_on="seqid")

    if len(loc) == 2:
        # Need to handle the special case where a gene is overlapping both ends
        # of a circular contig.
        # This code assumes that the gene overlaps both ends of a circular contig
        # and won't work otherwise
        loc0 = loc.loc[0]
        loc1 = loc.loc[1]
        if loc0.strand != loc1.strand:
            raise Exception("Unsupported case of fragment gene on different strands")

        _, strand, start, stop = (int(i) for i in loc0.tolist())
        if start==1:
            fet1 = SeqFeature(FeatureLocation(start-1, stop+flanking, strand=strand))
            fet0 = SeqFeature(FeatureLocation(int(loc1.start-flanking-1), int(loc1.end), strand=strand))
        else:
            fet0 = SeqFeature(FeatureLocation(start-1-flanking, stop, strand=strand))
            fet1 = SeqFeature(FeatureLocation(int(loc1.start)-1, int(loc1.end)+flanking, strand=strand))
        extracted0 = fet0.extract(seq)
        extracted1 = fet1.extract(seq)
        extracted = extracted0+extracted1
        red_start = flanking
        red_stop = len(extracted0)+(fet1.location.end-fet1.location.start-flanking)
    elif len(loc) == 1:
        _, strand, start, stop = (int(i) for i in loc.loc[0].tolist())
        start -= 1
        if start < 50:
            start_w_flank = 0
            red_start = start
        else:
            start_w_flank = start - flanking
            red_start = 50

        if stop + flanking > len(seq):
            stop_w_flank = len(seq) - 1
        else:
            stop_w_flank = stop + flanking
        red_stop = red_start + stop - start
        fet = SeqFeature(FeatureLocation(start_w_flank, stop_w_flank, strand=strand))
        extracted = fet.extract(seq)
    else:
        raise Exception("Unsupported case of fragmented gene")
    return extracted[0:red_start] + "<font color='red'>" + \
        extracted[red_start:red_stop] + "</font>" + extracted[red_stop:]

def tab_get_pfam_annot(db, seqid):
    pfam_hits = db.get_pfam_hits_info(seqid)
    feature_viewer_fet = []
    pfam_grouped = pfam_hits.groupby(["pfam"])
    pfam_starts = pfam_grouped["start"].apply(list)
    pfam_ends = pfam_grouped["end"].apply(list)
    pfam_defs_df = db.get_pfam_def(pfam_hits.pfam.tolist())

    pfam_defs = []
    for pfam, starts in pfam_starts.items():
        ends = pfam_ends.loc[pfam]
        name = format_pfam(pfam)
        data = "[" + ",".join(f"{{x: {start}, y: {end}}}" for start, end in zip(starts, ends)) + "]"
        feature = (
            f"{{data: {data}, "
            f" name: \"{name}\", "
            "  color: \"#0F8292\","
            "  type : \"rect\","
            "}"
        )
        pfam_def = pfam_defs_df["def"].loc[pfam]
        pfam_defs.append((format_pfam(pfam, to_url=True), pfam_def))
        feature_viewer_fet.append(feature)
    return {"pfam_domains": "[" + ",".join(feature_viewer_fet) + "]",
            "pfam_def": pfam_defs}


def genomic_region_df_to_js(df, start, end, name=None):
    features = []
    for curr_seqid, data in df.iterrows():
        feature_name = ""
        if "gene" in data and not pd.isna(data.gene):
            feature_name = data.gene
        feature_type = data.type
        if data.is_pseudo:
            feature_type = "pseudo"

        prod = to_s(data["product"])
        features.append((
            f"{{start: {data.start_pos}, gene: {to_s(feature_name)}, end: {data.end_pos},"
            f"strand: {data.strand}, type: {to_s(feature_type)}, product: {prod},"
            f"locus_tag: {to_s(data.locus_tag)}}}"
        ))
    features_str = "[" + ",".join(features) + "]"
    genome_name = ""
    if name is not None:
        genome_name = f"name: {to_s(name)}, "
    return f"{{{genome_name} start: {start}, end: {end}, features: {features_str}}}"


def locusx_genomic_region(db, seqid, window):
    hsh_loc = db.get_gene_loc([seqid])
    strand, start, end = hsh_loc[seqid]
    window_start, window_stop = start - window, start + window

    hsh_organism = db.get_organism([seqid], id_type="seqid")
    bioentry, _, contig_size, _ = db.get_bioentry_list(seqid, search_on="seqid")
    qualifiers = db.get_bioentry_qualifiers(bioentry)
    is_circular = "circular" in qualifiers["value"].values
    df_seqids = db.get_features_location(bioentry, search_on="bioentry_id")

    if 2*window >= contig_size:
        window_start = 0
        window_stop = contig_size
    elif window_start<0 and not is_circular:
        window_start = 0
        df_seqids = df_seqids[df_seqids.start_pos < window_stop]
    elif window_stop>contig_size and not is_circular:
        window_stop = contig_size
        df_seqids = df_seqids[df_seqids.end_pos>window_start]
    elif window_start<0:
        # circular contig
        diff = contig_size+window_start
        mask_circled = (df_seqids.end_pos >= diff)
        mask_same = (df_seqids.start_pos <= window_stop)

        df_seqids.loc[mask_same, "start_pos"] -= window_start
        df_seqids.loc[mask_same, "end_pos"] -= window_start
        df_seqids.loc[mask_circled, "start_pos"] -= diff
        df_seqids.loc[mask_circled, "end_pos"] -= diff
        df_seqids = pd.concat([df_seqids.loc[mask_same],
            df_seqids.loc[mask_circled]])
        window_stop -= window_start
        window_start = 0
    elif window_stop > contig_size:
        # circular contig
        diff = window_stop-contig_size

        mask_same = (df_seqids.end_pos >= window_start)
        mask_circled = (df_seqids.start_pos <= diff)

        df_seqids.loc[mask_same, "start_pos"] -= diff
        df_seqids.loc[mask_same, "end_pos"] -= diff
        df_seqids.loc[mask_circled, "start_pos"] += (contig_size-diff)
        df_seqids.loc[mask_circled, "end_pos"] += (contig_size-diff)
        df_seqids = pd.concat([df_seqids.loc[mask_same],
            df_seqids.loc[mask_circled]])
        window_start -= diff
        window_stop = contig_size
    else:
        df_seqids = df_seqids[(df_seqids.end_pos>window_start)]
        df_seqids = df_seqids[(df_seqids.start_pos < window_stop)]

    if len(df_seqids)!=len(df_seqids["seqfeature_id"].unique()):
        # This case may happen when a gene overlaps the break of a circular contig.
        # The location of this gene will be coded as join(...,...) in the gbk file
        # and stored as two separate genes with the same seqid in BioSQL.
        # If we want to display the whole contig as a continuous sequence, it is necessary
        # to detect this and manually merge this overlapping gene.
        grouped = df_seqids[["seqfeature_id", "strand", "end_pos",
            "start_pos"]].groupby("seqfeature_id")
        start = grouped["start_pos"].min()
        end = grouped["end_pos"].max()
        strands = df_seqids[["seqfeature_id", "strand"]].drop_duplicates("seqfeature_id")
        df_seqids = start.to_frame().join(end).join(strands.set_index("seqfeature_id"))
    else:
        df_seqids = df_seqids.set_index("seqfeature_id")

    # Some parts are redundant with get_features_location
    # those two function should be merged at some point
    infos = db.get_proteins_info(df_seqids.index.tolist(),
                                 to_return=["gene", "locus_tag", "product"], as_df=True,
                                 inc_non_CDS=True, inc_pseudo=True)
    cds_type = db.get_CDS_type(df_seqids.index.tolist())
    all_infos = infos.join(cds_type)
    all_infos = all_infos.join(df_seqids)
    return all_infos, window_start, window_stop

def tab_get_refseq_homologs(db, seqid):
    refseq_hits = db.get_refseq_hits([seqid]).set_index("match_id")
    refseq_hits_infos = db.get_refseq_matches_info(refseq_hits.index.tolist())
    all_infos = refseq_hits.join(refseq_hits_infos)

    header = ["Refseq accession", "Evalue", "Score", "ID(%)", "# gaps", "Len", "Description", "Organism"]
    entries = []
    for match_id, data in all_infos.iterrows():
        to_ncbi = format_refseqid_to_ncbi(data.accession)
        entries.append((to_ncbi, data.evalue, data.bitscore,
                        data.pident, data.gaps, data.length, data.description, data.organism))
    return {"n_refseq_homologs": len(refseq_hits),
            "refseq_headers": header,
            "blast_data": entries}


def locusx(request, locus=None, menu=True):
    biodb = settings.BIODB_DB_PATH
    db = DB.load_db(biodb, settings.BIODB_CONF)

    if locus is None:
        return render(request, 'chlamdb/locus.html', my_locals({"valid_id": False}))
    try:
        seqid, feature_type, is_pseudogene = db.get_seqid(locus_tag=locus,
                                                          feature_type=True)
    except Exception:
        return render(request, 'chlamdb/locus.html', my_locals({"valid": False}))
    else:
        valid_id = True

    page_title = f'Locus tag: {locus}'
    sequence = get_sequence(db, seqid, flanking=50)
    window_size = 8000
    all_infos, wd_start, wd_end = locusx_genomic_region(db, seqid, window=window_size)
    region_js = genomic_region_df_to_js(all_infos, wd_start, wd_end)
    genomic_region_ctx = {"genomic_region": region_js,
            "window_size": window_size*2}
    general_tab = tab_general(db, seqid)

    if feature_type!="CDS" or is_pseudogene:
        if is_pseudogene:
            feature_type = "Pseudogene"
        context = {
            "valid_id": valid_id,
            "menu": True,
            "seq": sequence,
            "feature_type": feature_type,
            "page_title": page_title,
            **general_tab,
            **genomic_region_ctx
        }
        return render(request, 'chlamdb/locus.html', my_locals(context))

    translation = db.get_translation(seqid)

    # a bit of an hack
    general_tab["length"] = len(translation)
    og_inf   = db.get_og_count([seqid], search_on="seqid")
    og_id    = int(og_inf.loc[seqid].orthogroup) # need to convert from numpy64 to int
    og_annot = db.get_genes_from_og(orthogroups=[og_id],
                                    terms=["locus_tag", "gene", "product", "length"])
    all_og_c = db.get_og_count([og_id], search_on="orthogroup")
    all_org = db.get_organism(og_annot.index.tolist())

    n_homologues = all_og_c.loc[og_id].sum() - 1
    og_size = n_homologues + 1
    og_num_genomes = len(set(all_org.values()))

    if n_homologues>1:
        og_conserv_ctx  = tab_og_conservation_tree(db, og_id, compare_to=seqid)
        homolog_tab_ctx = tab_homologs(db, og_annot, all_org, seqid, og_id)
        try:
            og_phylogeny_ctx = tab_og_phylogeny(db, og_id, compare_to=seqid)
        except NoPhylogenyException:
            og_phylogeny_ctx = {}
    else:
        og_conserv_ctx = {}
        homolog_tab_ctx = {"n_genomes": "1 genome"}
        og_phylogeny_ctx = {}

    kegg_ctx, cog_ctx, pfam_ctx, amr_ctx = {}, {}, {}, {}
    diamond_matches_ctx = {}
    swissprot_ctx = {}
    best_hit_phylo = {}
    if optional2status.get("KEGG", False):
        taxids = db.get_organism([seqid], as_taxid=True)
        kegg_ctx = og_tab_get_kegg_annot(db, [seqid], from_taxid=taxids[seqid])

    if optional2status.get("COG", False):
        cog_ctx = og_tab_get_cog_annot(db, [seqid])

    if optional2status.get("pfam", False):
        pfam_ctx = tab_get_pfam_annot(db, [seqid])

    if optional2status.get("BLAST_swissprot", False):
        swissprot_ctx = locus_tab_swissprot_hits(db, seqid)

    if optional2status.get("BLAST_database", False):
        diamond_matches_ctx = tab_get_refseq_homologs(db, seqid)

    if optional2status.get("BBH_phylogenies", False):
        best_hit_phylo = tab_og_best_hits(db, og_id, locus=locus)

    if optional2status.get("AMR", False):
        amr_ctx = og_tab_get_amr_annot(db, [seqid])

    context = {
        "valid_id": valid_id,
        "menu": True,
        "n_homologues": n_homologues,
        "og_id": format_orthogroup(og_id, to_url=True),
        "og_size": og_size,
        "og_num_genomes": og_num_genomes,
        "translation": translation,
        "seq": sequence,
        "sequence_type": feature_type,
        "locus": locus,
        "feature_type": "CDS",
        "page_title": page_title,
        **cog_ctx,
        **kegg_ctx,
        **homolog_tab_ctx,
        **general_tab,
        **og_conserv_ctx,
        **og_phylogeny_ctx,
        **pfam_ctx,
        **genomic_region_ctx,
        **diamond_matches_ctx,
        **swissprot_ctx,
        **best_hit_phylo,
        **amr_ctx
    }
    return render(request, 'chlamdb/locus.html', my_locals(context))


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
        {"label": f"{i.name}: {i.description} ({i.entry_type})",
         "value": f"{i.name}: {i.description}"}
        for i in results
    ]
    return JsonResponse(data, safe=False)

# NOTE: should refactor this code to avoid duplicated code


def search_bar(request):
    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    option2status = db.get_config_table()
    index = sb.ChlamdbIndex.use_index(settings.SEARCH_INDEX)
    user_query = request.GET.get("search2")

    results = list(index.search(user_query, limit=None))

    if len(results) == 0:
        ctx = {"search_failed": True, "search_term": user_query}
        return render(request, "chlamdb/search.html", my_locals(ctx))

    has_ko = option2status.get("KEGG", False)
    has_cog = option2status.get("COG", False)
    has_pfam = option2status.get("pfam", False)

    genes, cog, ko, pfam, pat, mod = [], [], [], [], [], []
    for result in results:
        if result.entry_type == sb.EntryTypes.Gene:
            locus_tag = format_locus(result.locus_tag, to_url=True)
            gene = str_if_none(result.name)
            product = str_if_none(result.description)
            genes.append([locus_tag, gene, product, result.organism])
        elif result.entry_type == sb.EntryTypes.COG and has_cog:
            cog.append([format_cog(None, base=result.name, as_url=True), result.description])
        elif result.entry_type == sb.EntryTypes.KO and has_ko:
            ko.append([format_ko(None, base=result.name, as_url=True), result.description])
        elif result.entry_type == sb.EntryTypes.PFAM and has_pfam:
            pfam.append([format_pfam(None, base=result.name, to_url=True), result.description])
        elif result.entry_type == sb.EntryTypes.Module and has_ko:
            mod.append([format_module(None, base=result.name, to_url=True), result.description])
        elif result.entry_type == sb.EntryTypes.Pathway and has_ko:
            pat.append([format_pathway(None, base=result.name, to_url=True), result.description])

    gene_active, cogs_active, ko_active, pfam_active, pat_active, mod_active = "active", "", "", "", "", ""
    if len(genes) == 0:
        gene_active = ""
        if len(cog) > 0:
            cogs_active = "active"
        elif len(ko) > 0:
            ko_active = "active"
        elif len(pfam) > 0:
            pfam_active = "active"
        elif len(pat) > 0:
            pat_active = "active"
        elif len(mod) > 0:
            mod_active = "active"

    genes_headers = ["Accession", "Gene", "Product", "Organism"]
    cog_headers = ["COG", "Description"]
    ko_headers = ["KO", "Description"]
    pfam_headers = ["PFAM domain", "Description"]
    pat_headers = ["KEGG Pathway", "Description"]
    mod_headers = ["KEGG Module", "Description"]
    ctx = {"search_term": user_query,
           "gene_active": gene_active,
           "cogs_active": cogs_active,
           "ko_active": ko_active,
           "pfam_active": pfam_active,
           "pat_active": pat_active,
           "mod_active": mod_active,
           "genes_headers": genes_headers,
           "genes": genes,
           "cog_headers": cog_headers,
           "pfam_headers": pfam_headers,
           "ko_headers": ko_headers,
           "pat_headers": pat_headers,
           "mod_headers": mod_headers,
           "modules": mod,
           "pathways": pat,
           "cogs": cog,
           "pfam": pfam,
           "ko": ko,
           "pat": pat,
           "mod": mod}
    return render(request, "chlamdb/search.html", my_locals(ctx))


def get_all_prot_infos(db, seqids, orthogroups):
    hsh_gene_locs = db.get_gene_loc(seqids)
    hsh_prot_infos = db.get_proteins_info(seqids)
    hsh_organisms = db.get_organism(seqids)
    group_count = set()
    all_locus_data = []

    for index, seqid in enumerate(seqids):
        # NOTE: all seqids are attributed an orthogroup, the case where
        # seqid is not in orthogroups should therefore not arise.
        og = orthogroups.loc[seqid].orthogroup
        fmt_orthogroup = format_orthogroup(og, to_url=True)
        group_count.add(fmt_orthogroup)
        strand, start, end = hsh_gene_locs[seqid]
        organism = hsh_organisms[seqid]
        locus, prot_id, gene, product = hsh_prot_infos[seqid]
        if gene is None:
            gene = ""
        data = (index, fmt_orthogroup, locus, prot_id, start, end, strand, gene, product, organism)
        all_locus_data.append(data)
    return all_locus_data, group_count


def format_orthogroup(og, to_url=False, from_str=False):
    base_str = og
    if not from_str:
        base_str = f"group_{og}"
    if to_url:
        return f"<a href=\"/orthogroup/{base_str}\">{base_str}</a>"
    return base_str


def format_locus(locus, to_url=True):
    if to_url:
        return f"<a href=\"/locusx/{locus}\">{locus}</a>"
    return locus


class FamCogColorFunc:
    def __init__(self, og, red_color):
        self.og = og
        self.red_color = red_color

    def get_color(self, taxid):
        if (self.og, taxid) in self.red_color:
            return "#FA5858"
        else:
            return EteTree.GREEN


def tab_gen_profile_tree(db, main_series, header, intersect):
    """
     Generate the tree from the profiles tab in the pfam/ko/cog pages:
     -ref_tree: the phylogenetic tree
     - main_series: the cog/ko/pfam count per taxid
     - header: the header of the main_series in the tree
     -intersect: a dataframe containing the seqid, taxid and orthogroups of the pfam/cog/ko hits
    """

    # the (group, taxid) in this dataframe are those that should be colored in red
    # in the profile (correspondance between a cog entry and an orthogroup)
    unique_og = intersect.orthogroup.unique().tolist()
    red_color = set(tuple(entry) for entry in intersect.to_numpy())
    df_og_count = db.get_og_count(list(unique_og), search_on="orthogroup").T
    ref_tree = db.get_reference_phylogeny()
    ref_names = db.get_genomes_description().description.to_dict()

    tree = Tree(ref_tree)
    R = tree.get_midpoint_outgroup()
    if R is not None:
        tree.set_outgroup(R)
    tree.ladderize()
    e_tree = EteTree(tree)
    e_tree.rename_leaves(ref_names)

    face_params = {"color": EteTree.RED}
    e_tree.add_column(SimpleColorColumn.fromSeries(main_series,
                                                   header=header, face_params=face_params))

    for og in df_og_count:
        og_serie = df_og_count[og]
        color_chooser = FamCogColorFunc(og, red_color)
        col_column = SimpleColorColumn(og_serie.to_dict(), header=format_orthogroup(og),
                                       col_func=color_chooser.get_color)
        e_tree.add_column(col_column)
    return e_tree


# TODO : add error handling
def fam_cog(request, cog_id):
    page_title = page2title["fam_cog"]

    biodb_path = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb_path)
    cog_id = int(cog_id[3:])

    if request.method != "GET":
        return render(request, 'chlamdb/fam.html', my_locals(locals()))

    df_seqid_to_cog = db.get_cog_hits([cog_id], indexing="seqid", search_on="cog", keep_taxid=True)
    if len(df_seqid_to_cog) == 0:
        return render(request, 'chlamdb/fam.html', {"msg": f"No entry for {format_cog(cog_id)}"})

    seqids = df_seqid_to_cog.index.tolist()

    orthogroups = db.get_og_count(seqids, search_on="seqid", keep_taxid=True)
    cog_info = db.get_cog_summaries([cog_id], only_cog_desc=True, as_df=True)
    all_locus_data, group_count = get_all_prot_infos(db, seqids, orthogroups)
    ref_tree = db.get_reference_phylogeny()
    cog_func = db.get_cog_code_description()

    df_cog_count = df_seqid_to_cog.groupby(["taxid"]).count()
    fam = format_cog(cog_id)
    e_tree = tab_gen_profile_tree(db, df_cog_count.cog, format_cog(cog_id), orthogroups)
    asset_path = f"/temp/fam_tree_{cog_id}.svg"
    path = settings.BASE_DIR + "/assets/" + asset_path
    e_tree.render(path, dpi=500)

    func, cog_description = cog_info.loc[cog_id]
    info_func = "<br>".join((cog_func[code] for code in func))
    type = "cog"

    info = [info_func, cog_description]
    menu = True
    envoi = True
    return render(request, 'chlamdb/fam.html', my_locals(locals()))


def format_module(mod_id, base=None, to_url=False):
    if base is None:
        formated = f"M{mod_id:05d}"
    else:
        formated = base

    if to_url:
        return f"<a href=/KEGG_module_map/{formated}>{formated}</a>"
    return formated


def fam_ko(request, ko_str):
    page_title = page2title["fam_ko"]

    ko_id = int(ko_str[len("K"):])
    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)

    df_ko_hits = db.get_ko_hits([ko_id], search_on="ko", indexing="seqid", keep_taxid=True)
    seqids = df_ko_hits.index.tolist()
    seqid_to_og = db.get_og_count(seqids, search_on="seqid", keep_taxid=True)
    red_color = set(tuple(entry) for entry in seqid_to_og.to_numpy())

    pathways = db.get_ko_pathways([ko_id])
    modules = db.get_ko_modules([ko_id])
    modules_id = [mod_id for key, values in modules.items() for mod_id, desc in values]
    modules_data = db.get_modules_info(modules_id)
    ko_desc = db.get_ko_desc([ko_id])[ko_id]
    all_locus_data, group_count = get_all_prot_infos(db, seqids, seqid_to_og)

    pathway_data = format_ko_path(pathways, ko_id, as_list=True)
    module_data = [(format_ko_module(mod_id), cat, mod_desc)
                   for mod_id, mod_desc, mod_def, path, cat in modules_data]
    df_ko_count = df_ko_hits.groupby(["taxid"]).count()

    fam = format_ko(ko_id)
    e_tree = tab_gen_profile_tree(db, df_ko_count.ko, fam, seqid_to_og)
    asset_path = f"/temp/fam_tree_{fam}.svg"
    path = settings.BASE_DIR + f"/assets/{asset_path}"

    e_tree.render(path, dpi=500)
    type = "ko"
    menu = True
    envoi = True
    return render(request, 'chlamdb/fam.html', my_locals(locals()))


def fam_pfam(request, pfam):
    page_title = page2title["fam_pfam"]

    context = {
        "type": "pfam",
        "menu": True,
        "envoi": True,
        "fam": pfam
    }
    if len(pfam) < 2 or not pfam.startswith("PF"):
        # add error message
        return render(request, 'chlamdb/fam.html', my_locals(context))
    try:
        pfam_id = int(pfam[len("PF"):])
    except Exception:
        return render(request, 'chlamdb/fam.html', my_locals(context))

    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    pfam_hits = db.get_pfam_hits([pfam_id], search_on="pfam", indexing="seqid", keep_taxid=True)
    seqids = pfam_hits.seqid.unique().tolist()
    orthogroups = db.get_og_count(seqids, search_on="seqid", keep_taxid=True)
    all_locus_data, group_count = get_all_prot_infos(db, seqids, orthogroups)
    infos = db.get_pfam_def([pfam_id])

    e_tree = tab_gen_profile_tree(db, pfam_hits.groupby(["taxid"])["pfam"].count(),
                                  pfam, orthogroups)

    asset_path = f"/temp/fam_tree_{pfam_id}.svg"
    path = settings.BASE_DIR + "/assets/" + asset_path
    e_tree.render(path, dpi=500)

    context["all_locus_data"] = all_locus_data
    context["group_count"] = group_count
    context["info"] = [infos.loc[pfam_id]["def"]]
    context["asset_path"] = asset_path
    context["page_title"] = page_title
    return render(request, 'chlamdb/fam.html', my_locals(context))


def COG_phylo_heatmap(request, frequency):
    biodb = settings.BIODB_DB_PATH

    page_title = page2title["COG_phylo_heatmap"]

    if request.method != "GET":
        return render(request, 'chlamdb/COG_phylo_heatmap.html', my_locals(locals()))
    freq = frequency != "False"

    db = DB.load_db_from_name(biodb)
    tree = db.get_reference_phylogeny()
    descr = db.get_genomes_description()
    all_taxids = descr.index.tolist()

    all_cog_hits = db.get_cog_hits(all_taxids, indexing="taxid", search_on="taxid")
    all_cog_funcs = db.get_cog_summaries(all_cog_hits.index.unique().tolist(),
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

    funcs_descr = db.get_cog_code_description()
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
        col = SimpleColorColumn.fromSeries(func_count,
                                           header=detailed_func + "(" + func + ")", color_gradient=True)
        e_tree.add_column(col)

    freq = frequency
    path = settings.BASE_DIR + f"/assets/temp/COG_tree_{freq}.svg"
    asset_path = f"/temp/COG_tree_{freq}.svg"
    e_tree.render(path, dpi=600)
    envoi = True
    return render(request, 'chlamdb/COG_phylo_heatmap.html', my_locals(locals()))


class KOModuleChooser:
    def __init__(self, hsh):
        self.hsh = hsh

    def get_color(self, index):
        return self.hsh[index]


def KEGG_module_map(request, module_name):
    page_title = page2title["KEGG_module_map"]

    if request.method != "GET":
        return render(request, 'chlamdb/KEGG_module_map.html', my_locals(locals()))
    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)

    try:
        module_id = int(module_name[len("M"):])
    except Exception:
        # add error message: module not formated correctly
        valid_id = False
        return render(request, 'chlamdb/KEGG_module_map.html', my_locals(locals()))

    module_infos = db.get_modules_info([module_id])
    if len(module_infos) != 1:
        # add error message
        return render(request, 'chlamdb/KEGG_module_map.html', my_locals(locals()))
    else:
        mod_id, module_descr, module_def, cat, sub_cat = module_infos[0]

    parser = ModuleParser(module_def)
    expr_tree = parser.parse()
    ko_ids = db.get_module_kos(module_id)
    mat = db.get_ko_hits(ko_ids, search_on="ko", indexing="taxid").T
    map_data = [(format_ko(ko_id), ko_desc) for ko_id, ko_desc in db.get_ko_desc(ko_ids).items()]
    if mat.empty:
        # should add an error message: no gene was associated for any
        # of the KO of the current module
        envoi = True
        menu = True
        valid_id = True
        return render(request, 'chlamdb/KEGG_module_map.html', my_locals(locals()))

    seqids = db.get_ko_hits(ko_ids, search_on="ko", indexing="seqid")
    associated_ogs = db.get_og_count(seqids.index.tolist(), search_on="seqid")
    og_taxid = db.get_og_count(associated_ogs.orthogroup.unique().tolist(), search_on="orthogroup").T
    ko_to_og_mapping = seqids.join(associated_ogs).groupby(["ko"])["orthogroup"].unique()
    leaf_to_name = db.get_genomes_description().description.to_dict()

    tree = Tree(db.get_reference_phylogeny())
    R = tree.get_midpoint_outgroup()
    tree.set_outgroup(R)
    tree.ladderize()
    e_tree = EteTree(tree)

    hsh_pres = collections.defaultdict(dict)
    for ko in ko_ids:
        if ko not in mat.columns:
            e_tree.add_column(SimpleColorColumn({}, header=format_ko(ko), default_val="-"))
            continue
        if ko not in ko_to_og_mapping.index:
            e_tree.add_column(SimpleColorColumn.fromSeries(mat[ko], header=format_ko(ko)))
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
    completeness = ModuleCompletenessColumn(hsh_n_missing, "", add_missing=False)
    e_tree.add_column(completeness)
    e_tree.rename_leaves(leaf_to_name)

    big = len(mat.columns) >= 40
    dpi = 800 if big else 1200
    path = settings.BASE_DIR + '/assets/temp/KEGG_tree_%s.svg' % module_name
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
    og_taxid = db.get_og_count(associated_ogs.orthogroup.unique().tolist(), search_on="orthogroup").T
    ko_to_og_mapping = seqids.join(associated_ogs).groupby(["ko"])["orthogroup"].unique()
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
            e_tree.add_column(SimpleColorColumn.fromSeries(mat[ko], header=format_ko(ko)))
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
                   "error_title": "Error"}
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
        ctx = {"error": True, "error_title": "No hits for this pathway"}
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
    path = settings.BASE_DIR + f"/assets/temp/{map_name}.svg"
    asset_path = f"/temp/{map_name}.svg"
    e_tree.render(path, dpi=800)
    ctx = {"pathway_num": kos.iloc[0].pathway,
           "pathway": kos.iloc[0].description,
           "header": header,
           "data": data,
           "asset_path": asset_path,
           "url": map_name + "+" + "+".join(all_kos),
           "envoi": True,
           "error": False,
           "page_titleA": page_title
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


def cog_venn_subset(request, category):
    # Note: add error handling

    biodb = settings.BIODB_DB_PATH
    db = DB.load_db(biodb, settings.BIODB_CONF)

    targets = [int(i) for i in request.GET.getlist('h')]
    if len(targets) > 5:
        targets = targets[0:6]

    cog_hits = db.get_cog_hits(targets, indexing="taxid", search_on="taxid")
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
    ko_modules = db.get_ko_modules(ko_counts["KO"].values.tolist(), as_pandas=True, compact=True)
    ko_modules_info = db.get_modules_info(ko_modules["module_id"].unique().tolist(), as_pandas=True)
    filtered_modules = ko_modules_info[ko_modules_info.subcat == category]
    selected_kos = filtered_modules.merge(ko_modules, left_on="module_id",
                                          right_on="module_id", how="inner")["ko_id"].unique()
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


def to_s(f):
    return "\"" + str(f) + "\""


def js_bioentries_to_description(hsh):
    taxon_map = 'var taxon2description = { '
    mid = ",".join(f"{to_s(bioentry)}: {to_s(description)}" for bioentry, description in hsh.items())
    return taxon_map + mid + "};"


def module_barchart(request):
    biodb = settings.BIODB_DB_PATH
    db = DB.load_db(biodb, settings.BIODB_CONF)
    page_title = page2title["module_barchart"]

    venn_form_class = make_venn_from(db)
    if request.method != "POST":
        form = venn_form_class()
        return render(request, 'chlamdb/module_barplot.html', my_locals(locals()))

    form = venn_form_class(request.POST)
    if not form.is_valid():
        form = venn_form_class()
        return render(request, 'chlamdb/module_barplot.html', my_locals(locals()))

    taxids = form.get_taxids()
    taxon2description = db.get_genomes_description().description.to_dict()

    ko_counts = db.get_ko_count(taxids, keep_seqids=True, as_multi=False)
    ko_ids = ko_counts.KO.unique()
    ko_module_ids = db.get_ko_modules(ko_ids.tolist(), as_pandas=True, compact=True)
    ko_modules_info = db.get_modules_info(ko_module_ids["module_id"].unique().tolist(), as_pandas=True)

    merged = ko_counts.merge(ko_module_ids, left_on="KO", right_on="ko_id", how="inner")
    merged = merged.merge(ko_modules_info, left_on="module_id", right_on="module_id", how="inner")
    cat_count = merged[["taxid", "subcat", "seqid"]].groupby(["taxid", "subcat"]).nunique()
    subcategories_list = cat_count.index.unique(level="subcat").to_list()
    subcategories = ",".join(f"{to_s(cat)}" for cat in subcategories_list)
    labels = f"[{subcategories}]"

    taxon_map = js_bioentries_to_description(taxon2description)

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

    taxids = "?" + "&".join((f"h={i}" for i in taxids))
    series = "[" + ",".join(series_data) + "]"
    envoi = True
    form = venn_form_class()
    return render(request, 'chlamdb/module_barplot.html', my_locals(locals()))


def orthogroup_list_cog_barchart(request, accessions=False):
    biodb = settings.BIODB
    if accessions == 'False' or accessions == 'F':
        accessions = False

    orthogroup_list = [i for i in request.GET.getlist('h')]

    reference = request.GET.getlist('ref')[0]

    series, \
        labels, \
        serie_all_counts, \
        serie_target_counts, \
        series_counts, \
        labels_counts, \
        category_description, \
        category_map, \
        n_missing_cog, \
        missing_cog_list = orthogroup2cog_series(orthogroup_list, reference_taxon=reference, accessions=accessions)

    no_cogs_url = "?g=" + ('&g=').join(missing_cog_list)
    orthogroups_url = '?h=' + ('&h=').join(orthogroup_list)

    return render(request, 'chlamdb/orthogroup_list_cog_barchart.html', my_locals(locals()))


def cog_barchart(request):
    # Not too happy with this code
    # should be reviewed and the old hashes replaced by pandas
    # operations
    page_title = page2title["cog_barchart"]

    biodb = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb)
    venn_form_class = make_venn_from(db)

    if request.method != 'POST':
        form = venn_form_class()
        return render(request, 'chlamdb/cog_barplot.html', my_locals(locals()))

    form = venn_form_class(request.POST)

    if not form.is_valid():
        # add error message
        return render(request, 'chlamdb/cog_barplot.html', my_locals(locals()))

    target_bioentries = form.get_taxids()

    hsh_counts = db.get_cog_counts_per_category(target_bioentries)
    taxon2description = db.get_genomes_description().description.to_dict()
    category_dico = db.get_cog_code_description()

    # create a dictionnary to convert cog category description and one letter code
    category_map = 'var description2category = {'
    map_lst = (f"\"{func_descr}\": \"{func}\"" for func, func_descr in category_dico.items())
    category_map = category_map + ",".join(map_lst) + '};'

    taxon_map = 'var taxon2description = {'
    taxon_map_lst = (f"\"{target}\": \"{taxon2description[target]}\"" for target in target_bioentries)
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

        all_series_templates.append(one_serie_template % (taxon, ','.join(one_category_list)))

    taxids = "?" + "&".join((f"h={i}" for i in target_bioentries))
    series = serie_template % ''.join(all_series_templates)
    labels = labels_template % ('"' + '","'.join(all_categories) + '"')
    envoi = True
    return render(request, 'chlamdb/cog_barplot.html', my_locals(locals()))


def pan_genome(request, type):
    biodb = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb)
    page_title = page2title[f"pan_genome_{type}"]
    venn_form_class = make_venn_from(db, plasmid=False)

    if request.method != "POST":
        form = venn_form_class()
        return render(request, 'chlamdb/pan_genome.html', my_locals(locals()))

    form = venn_form_class(request.POST)
    if not form.is_valid():
        # should add an error message
        form = venn_form_class()
        return render(request, 'chlamdb/pan_genome.html', my_locals(locals()))
    taxids = form.get_taxids()

    if type == "COG":
        df_hits = db.get_cog_hits(taxids, search_on="taxid", indexing="taxid")
        type_txt = "COG orthologs"
    elif type == "orthology":
        df_hits = db.get_og_count(taxids, search_on="taxid")
        type_txt = "orthologs"
    elif type == "ko":
        df_hits = db.get_ko_hits(taxids, search_on="taxid")
        type_txt = "KEGG orthologs"
    elif type == "Pfam":
        df_hits = db.get_pfam_hits(taxids, search_on="taxid")
        type_txt = "PFAM domains"
    else:
        form = venn_form_class()
        return render(request, 'chlamdb/pan_genome.html', my_locals(locals()))

    unique, counts = np.unique(np.count_nonzero(df_hits, axis=1), return_counts=True)
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

    envoi = True
    return render(request, 'chlamdb/pan_genome.html', my_locals(locals()))


blast_input_dir = {"blastp": "faa",
                   "tblastn": "ffn",
                   "blastn_fna": "ffn",
                   "blastn_ffn": "fna",
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
            hsp = hit.hsps[0]
            scores = hits[query]
            scores.append((hit.accession, 100.0 * hsp.identities / hsp.align_length))
            accessions.add(hit.accession)

    if blast_type in ["blastp", "blastx", "blastn_ffn"]:
        acc_to_taxid = db.get_taxid_from_accession(list(accessions), look_on="locus_tag")
    elif blast_type in ["tblastn", "blastn_fna"]:
        acc_to_taxid = db.get_taxid_from_accession(list(accessions), look_on="contig")
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
        col = SimpleColorColumn(hsh_values, header=query, color_gradient=True,
                                default_val="-", gradient_value_range=[min_val, max_val])
        e_tree.add_column(col)

    base_file_name = time.strftime("blast_%d_%m_%y_%H_%M.svg", time.gmtime())
    path = settings.BASE_DIR + f"/assets/temp/{base_file_name}"
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
        return render(request, 'chlamdb/blast.html', my_locals({"form": form, "page_title": page_title}))

    form = blast_form_class(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/blast.html', my_locals({"form": form}))

    input_sequence = form.cleaned_data['blast_input']
    number_blast_hits = form.cleaned_data['max_number_of_hits']
    blast_type = form.cleaned_data['blast']
    target = form.get_target()
    no_query_name = False

    if '>' in input_sequence:
        try:
            records = [i for i in SeqIO.parse(StringIO(input_sequence), 'fasta')]
            for record in records:
                if len(record.seq) == 0:
                    context = {"error_message": "Empty sequence in input",
                               "error_title": "Query format error", "envoi": True,
                               "form": form, "wrong_format": True}
                    return render(request, 'chlamdb/blast.html', my_locals(context))
        except Exception:
            context = {"error_message": "Error while parsing the fasta query",
                       "error_title": "Query format error",
                       "envoi": True, "form": form, "wrong_format": True}
            return render(request, 'chlamdb/blast.html', my_locals(context))
    else:
        no_query_name = True
        input_sequence = "".join(input_sequence.split()).upper()
        records = [SeqRecord(Seq(input_sequence))]

    dna = set("ATGCNRYKMSWBDHV")
    prot = set('ACDEFGHIKLMNPQRSTVWYXZJOU')
    sequence_set = set()
    for rec in records:
        sequence_set = sequence_set.union(set(rec.seq.upper()))
    check_seq_DNA = sequence_set - dna
    check_seq_prot = sequence_set - prot

    if check_seq_prot and blast_type in ["blastp", "tblastn"]:
        wrong_chars = ", ".join(check_seq_prot)
        if len(check_seq_prot) > 1:
            error_message = f"Unexpected characters in amino-acid query: {wrong_chars}"
        else:
            error_message = f"Unexpected character in amino-acid query: {wrong_chars}"
        context = {"error_message": error_message, "error_title": "Query format error",
                   "envoi": True, "form": form, "wrong_format": True}
        return render(request, 'chlamdb/blast.html', my_locals(context))
    elif check_seq_DNA and blast_type in ["blastn", "blastn_ffn", "blast_fna", "blastx"]:
        wrong_chars = ", ".join(check_seq_DNA)
        if len(check_seq_DNA) > 1:
            error_message = f"Unexpected characters in nucleotide query: {wrong_chars}"
        else:
            error_message = f"Unexpected character in nucleotide query: {wrong_chars}"
        context = {"error_message": error_message, "wrong_format": True,
                   "error_title": "Query format error", "envoi": True, "form": form}
        return render(request, 'chlamdb/blast.html', my_locals(context))

    query_file = NamedTemporaryFile(mode='w')
    SeqIO.write(records, query_file, "fasta")
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
        context = {"error_message": "Error in blast " + str(e), "wrong_format": True,
                   "error_title": "Blast error", "envoi": True, "form": form}
        return render(request, 'chlamdb/blast.html', my_locals(context))

    if blast_stdout.find("No hits found") != -1:
        blast_no_hits = blast_stdout
    elif len(blast_stderr) != 0:
        blast_err = blast_stderr
    else:
        if target == "all":
            asset_path = gen_blast_heatmap(db, blast_stdout, blast_type, no_query_name)
        rand_id = id_generator(6)
        blast_file_l = settings.BASE_DIR + '/assets/temp/%s.xml' % rand_id
        f = open(blast_file_l, 'w')
        f.write(blast_stdout)
        f.close()
        asset_blast_path = '/assets/temp/%s.xml' % rand_id
        js_out = True
        phylo_distrib = target == "all"
    envoi = True
    return render(request, 'chlamdb/blast.html', my_locals(locals()))


def format_gene(gene):
    if pd.isna(gene):
        return "-"
    else:
        return gene


def format_taxid_to_ncbi(organism, taxid):
    val = (
        f"""<a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}" target="_top">"""
        f"""{organism}</a>"""
    )
    return val


def format_refseqid_to_ncbi(seqid):
    return f"<a href=\"http://www.ncbi.nlm.nih.gov/protein/{seqid}\">{seqid}</a>"


def format_gene_to_ncbi_hmm(gene_and_hmmid):
    gene, hmm_id = gene_and_hmmid
    if hmm_id:
        # The hmm_id contains a version number, which is not in the ncbi URL.
        hmm_id = hmm_id.rsplit(".", 1)[0]
        return f"<a href=\"https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/{hmm_id}\">{gene}</a>"  # noqa
    return gene


def locus_tab_swissprot_hits(db, seqid):
    swissprot_homologs = db.get_swissprot_homologs([seqid])

    blast_data = []
    header = ["Swissprot accession", "Eval", "Score", "ID (%)", "N gaps",
              "Alignment length", "Annot score", "Gene", "Description", "Organism"]

    for num, data in swissprot_homologs.iterrows():
        line = [data.accession, data.evalue, data.bitscore, data.perc_id,
                data.gaps, data.match_len, data.pe, data.gene, data.definition,
                format_taxid_to_ncbi(data.organism, data.taxid)]
        blast_data.append(line)

    ctx = {
        "n_swissprot_hits": len(swissprot_homologs),
        "swissprot_blast_data": blast_data,
        "swissprot_headers": header
    }
    return ctx


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


def plot_region(request):
    max_region_size = 100000
    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    page_title = page2title["plot_region"]
    form_class = make_plot_form(db)

    if request.method != "POST":
        form = form_class()
        return render(request, 'chlamdb/plot_region.html', my_locals({"form": form, "page_title": page_title}))

    form = form_class(request.POST)
    if not form.is_valid():
        # add error message
        return render(request, 'chlamdb/plot_region.html', my_locals({"form": form, "page_title": page_title}))

    accession = form.get_accession()
    prot_info = db.get_proteins_info(ids=[accession], search_on="locus_tag", as_df=True)
    if prot_info.empty:
        context = {"form": form,
                   "errors": ["Accession not found"],
                   "error": True,
                   "page_title": page_title}
        return render(request, 'chlamdb/plot_region.html', my_locals(context))

    genomes = form.get_genomes()
    all_homologs = form.get_all_homologs()
    seqid = int(prot_info.index[0])
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

    try:
        region_size = form.get_region_size()
        if region_size > max_region_size or region_size < 5000:
            context = {
                "form": form,
                "error": True,
                "errors": [f"Region size should be between 5000 and {max_region_size} bp"],
                "page_title": page_title
            }
            return render(request, 'chlamdb/plot_region.html', my_locals(context))
    except ValueError:
        context = {"form": form,
                   "error": True,
                   "errors": ["Wrong format for region size"],
                   "page_title": page_title}
        return render(request, 'chlamdb/plot_region.html', my_locals(context))

    locus_tags = db.get_proteins_info(seqids, to_return=["locus_tag"], as_df=True)
    to_highlight = locus_tags["locus_tag"].tolist()
    all_regions = []
    connections = []
    prev_infos = None
    all_identities = []

    genomic_regions = []
    for seqid in seqids:
        region, start, end = locusx_genomic_region(db, int(seqid), region_size / 2)
        genomic_regions.append([seqid, region, start, end])

    # remove overlapping regions (e.g. if two matches are on the same
    # region, avoid displaying this region twice).
    # XXX : coalesce_regions does not do what it is intended to do, to fix
    filtered_regions = genomic_regions # coalesce_regions(genomic_regions, seqids)

    for genomic_region in filtered_regions:
        seqid, region, start, end = genomic_region
        if region["strand"].loc[seqid] * ref_strand == -1:
            mean_val = (end + start) / 2
            region["strand"] *= -1
            dist_vector_start = region["start_pos"]-mean_val
            len_vector = region["end_pos"]-region["start_pos"]
            region["end_pos"] = region["start_pos"]-2*dist_vector_start
            region["start_pos"] = region["end_pos"] - len_vector
        ogs = db.get_og_count(region.index.tolist(), search_on="seqid")
        # need to reset index to keep seqid in the next merge
        region = region.join(ogs).reset_index()
            
        if not prev_infos is None:
            # BM: horrible code (I wrote it, I should know).
            # would be nice to refactor it in a more efficient and clean way.

            # We need to drop na from orthogroup, as some genes, like rRNA or tRNA
            # are not assigned any orthogroup
            common_og = region.dropna(subset=["orthogroup"]).merge(
                    prev_infos, on="orthogroup")[["locus_tag_x",
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
            related = (f"{to_s(loc)}: [" + ",".join(values) + "]" for loc, values in hsh_agg.items())
            connections.append("{" + ",".join(related) + "}")

        taxid = organisms.loc[seqid].taxid
        genome_name = hsh_description[taxid]
        js_val = genomic_region_df_to_js(region, start, end, genome_name)
        all_regions.append(js_val)
        prev_infos = region[["orthogroup", "locus_tag", "seqid"]]

    ctx = {"form": form, "genomic_regions": "[" + "\n,".join(all_regions) + "]",
           "window_size": max_region_size, "to_highlight": to_highlight, "envoi": True,
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

        js_code = get_circos_data(int(reference_taxon), [int(i) for i in target_taxons], highlight_og=og_list)

        envoi = True
        envoi_region = True

        return render(request, 'chlamdb/circos_main.html', my_locals(locals()))

    if request.method == 'GET':
        server, db = manipulate_biosqldb.load_db(biodb)

        reference_taxon = int(request.GET.getlist('ref')[0])

        if request.GET.getlist('t')[0] == '':
            # if no target list given, get the 10 closest genomes
            try:
                sql_order = 'select taxon_2 from comparative_tables_core_orthogroups_identity_msa '\
                            'where taxon_1=%s order by identity desc;' % (reference_taxon)
                ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order)]
                target_taxons = ordered_taxons[0:10]
            except Exception:
                sql_order = 'select taxon_2 from comparative_tables_shared_orthogroups '\
                            'where taxon_1=%s order by n_shared_orthogroups DESC;' % (reference_taxon)

                ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order)]
                target_taxons = ordered_taxons[0:10]
        else:
            target_taxons = [int(i) for i in request.GET.getlist('t')]
        highlight = request.GET.getlist('h')

        task = run_circos_main.delay(reference_taxon, target_taxons, highlight)
        task_id = task.id
        envoi_circos = True
        envoi_region = True

        # return HttpResponse(json.dumps({'task_id': task.id}), content_type='application/json')

    return render(request, 'chlamdb/circos_main.html', my_locals(locals()))


def get_circos_data(reference_taxon, target_taxons, highlight_og=False):
    biodb_path = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb_path)

    if highlight_og:
        df_genes = db.get_genes_from_og(highlight_og, taxon_ids=[reference_taxon], terms=["locus_tag"])
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

    seq_og = db.get_og_count(df_feature_location.index.to_list(), search_on="seqid")
    count_all_genomes = db.get_og_count(seq_og["orthogroup"].to_list(), search_on="orthogroup")
    orthogroup2count_all = count_all_genomes[count_all_genomes > 0].count(axis=1)
    homologs_count = df_feature_location.loc[df_feature_location.term_name == 'CDS']\
        .join(seq_og)\
        .reset_index()\
        .set_index("orthogroup")\
        .merge(orthogroup2count_all.rename('value'),
               left_index=True,
               right_index=True)[["bioentry_id", "start_pos", "end_pos", "value"]]

    df_identity = db.get_identity_closest_homolog(reference_taxon, target_taxons)\
                    .set_index(["target_taxid"])

    c = circosjs.CircosJs()

    c.add_contigs_data(df_bioentry)

    # sort taxons by number of homologs (from mot similar to most dissmilar)
    target_taxon_n_homologs = df_identity.groupby(["target_taxid"])\
                                         .count()["seqfeature_id_1"]\
                                         .sort_values(ascending=False)
    locus2seqfeature_id = db.get_hsh_locus_to_seqfeature_id()
    seqfeature_id2locus_tag = {value: key for key, value in locus2seqfeature_id.items()}
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
    c.add_histogram_track(homologs_count, "n_genomes", radius_diff=0.1, sep=0.005, outer=True)
    c.add_gene_track(minus_strand, "minus", sep=0, radius_diff=0.04, highlight_list=locus_list)
    c.add_gene_track(plus_strand, "plus", sep=0, radius_diff=0.04, highlight_list=locus_list)

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
        c.add_heatmap_track(df_combined, f"target_{target_taxon}", color="comp", sep=sep, radius_diff=0.04)

    c.add_line_track(df_bioentry, "GC_content", radius_diff=0.12, fillcolor="green")

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
    # local_vars["form"] = form

    return render(request, 'chlamdb/circos.html', local_vars)


def alignment(request, input_fasta):

    handle = open(input_fasta, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        pass
    return render(request, 'chlamdb/alignment.html', my_locals(locals()))


def plot_heatmap(request, type):
    import plotly.graph_objects as go
    import scipy.cluster.hierarchy as shc
    from scipy.cluster import hierarchy

    page_title = page2title[f"plot_heatmap_{type}"]

    biodb = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb)
    form_class = make_venn_from(db)

    if request.method != "POST":
        form_venn = form_class()
        ctx = {"form_venn": form_venn, "type": type, "page_title": page_title}
        return render(request, 'chlamdb/plot_heatmap.html', my_locals(ctx))

    form_venn = form_class(request.POST)
    if not form_venn.is_valid():
        ctx = {"form_venn": form_venn, "type": type, "page_title": page_title}
        return render(request, 'chlamdb/plot_heatmap.html', my_locals(ctx))

    taxon_ids = form_venn.get_taxids()

    if len(taxon_ids) <= 1:
        error_message = "Please select at least two genomes"
        error_title = "Wrong input"
        ctx = {"error": True, "type": type, "type": type, "error_title": error_title,
               "form_venn": form_venn, "error_message": error_message, "page_title": page_title}
        return render(request, 'chlamdb/plot_heatmap.html', my_locals(ctx))

    if type == "COG":
        mat = db.get_cog_hits(taxon_ids, indexing="taxid", search_on="taxid")
        mat.index = [format_cog(i) for i in mat.index]
    elif type == "orthology":
        mat = db.get_og_count(taxon_ids)
        mat.index = [format_orthogroup(i) for i in mat.index]
    elif type == "ko":
        mat = db.get_ko_hits(taxon_ids)
        mat.index = [format_ko(i) for i in mat.index]
    elif type == "Pfam":
        mat = db.get_pfam_hits(taxon_ids)
        mat.index = [format_pfam(i) for i in mat.index]
    else:
        form_venn = form_class()
        return render(request, 'chlamdb/plot_heatmap.html', my_locals(locals()))

    target2description = db.get_genomes_description().description.to_dict()
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
    fig.update_xaxes(visible=False)

    html_plot = make_div(fig, div_id="heatmap")
    envoi_heatmap = True
    return render(request, 'chlamdb/plot_heatmap.html', my_locals(locals()))


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
        return render(request, 'chlamdb/kegg_genomes.html', my_locals(locals()))

    form = single_genome_form(request.POST)
    if not form.is_valid():
        form = single_genome_form()
        return render(request, 'chlamdb/kegg_genomes.html', my_locals(locals()))

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
        entry = (format_pathway(element, to_url=True, taxid=taxid), descr, count)
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
        return render(request, 'chlamdb/kegg_genomes_modules.html', my_locals(locals()))

    form = single_genome_form(request.POST)
    if not form.is_valid():
        form = single_genome_form()
        return render(request, 'chlamdb/kegg_genomes_modules.html', my_locals(locals()))

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
    modules_info = db.get_modules_info(grouped_count.index.tolist(), as_pandas=True)
    all_infos = modules_info.set_index("module_id").join(grouped_count)

    header = ["Module", "Category", "Sub-category", "Description"]
    entries = []

    taxons = []
    for taxid in taxids:
        taxons.append(genomes[taxid])

    for module_id, data in all_infos.iterrows():
        line = [format_module(module_id, to_url=True), data["cat"], data.subcat, data.descr]
        for taxid in taxids:
            line.append(data[taxid])
        entries.append(line)
    envoi_comp = True
    return render(request, 'chlamdb/kegg_genomes_modules.html', my_locals(locals()))


def kegg(request):
    page_title = page2title["kegg"]

    db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    module_overview_form = make_module_overview_form(db)
    form_cat = module_overview_form(request.POST)
    form_cat = module_overview_form()

    single_genome_form = make_single_genome_form(db)
    form_genome = single_genome_form(request.POST)
    form_genome = single_genome_form()

    module_overview_form = make_module_overview_form(db, True)
    form_subcat = module_overview_form(request.POST)
    form_subcat = module_overview_form()

    comp_metabo_form = make_metabo_from(db)
    form_module = comp_metabo_form(request.POST)
    form_module = comp_metabo_form()

    pathway_overview_form = make_pathway_overview_form(db)
    form_pathway = pathway_overview_form(request.POST)
    form_pathway = pathway_overview_form()

    return render(request, 'chlamdb/kegg.html', my_locals(locals()))


def kegg_module_subcat(request):
    page_title = page2title["kegg_module_subcat"]

    biodb = settings.BIODB_DB_PATH
    db = DB.load_db(biodb, settings.BIODB_CONF)

    module_overview_form = make_module_overview_form(db, True)
    if request.method != "POST":
        form = module_overview_form()
        return render(request, 'chlamdb/module_subcat.html', my_locals(locals()))

    form = module_overview_form(request.POST)
    if not form.is_valid():
        # TODO: add error message
        form = module_overview_form()
        return render(request, 'chlamdb/module_subcat.html', my_locals(locals()))

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

    unique_module_ids = ko_count_subcat.index.get_level_values("module_id").unique().tolist()
    if len(unique_module_ids) == 0:
        # add error message : no module found
        envoi = True
        return render(request, 'chlamdb/module_subcat.html', my_locals(locals()))

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
    path = settings.BASE_DIR + '/assets/temp/metabo_tree.svg'
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
        return render(request, 'chlamdb/module_overview.html', my_locals(locals()))

    form = module_overview_form(request.POST)
    if not form.is_valid():
        # TODO: add error message
        form = module_overview_form()
        return render(request, 'chlamdb/module_overview.html', my_locals(locals()))

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
        return render(request, 'chlamdb/module_overview.html', my_locals(locals()))

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
    path = settings.BASE_DIR + '/assets/temp/metabo_tree.svg'
    asset_path = '/temp/metabo_tree.svg'
    e_tree.render(path, dpi=500, w=800)

    envoi = True
    return render(request, 'chlamdb/module_overview.html', my_locals(locals()))


class TabularComparisonViewBase(View):

    template = 'chlamdb/tabular_comparison.html'
    hist_colour_index_shift = 0
    tab_name = "comp"

    def dispatch(self, request, *args, **kwargs):
        biodb_path = settings.BIODB_DB_PATH
        self.db = DB.load_db_from_name(biodb_path)
        self.page_title = page2title[self.view_name]
        self.comp_metabo_form = make_metabo_from(self.db)
        self.show_comparison_table = False
        self._hash_to_taxon_dict = None
        return super(TabularComparisonViewBase, self).dispatch(request, *args, **kwargs)

    def get(self, request, *args, **kwargs):
        self.form = self.comp_metabo_form()
        return render(request, self.template, self.context)

    def post(self, request, *args, **kwargs):
        self.form = self.comp_metabo_form(request.POST)
        if self.form.is_valid():
            self.show_comparison_table = True
            self.set_table_data()

        return render(request, self.template, self.context)

    @property
    def view_name(self):
        return f"{self.view_type}_comparison"

    @property
    def context(self):
        context = {
            "page_title": self.page_title,
            "form_title": self.form_title,
            "form_help": self.form_help,
            "form": self.form,
            "show_comparison_table": self.show_comparison_table,
            "tab_name": self.tab_name,
            "view_name": self.view_name,
            "view_type": self.view_type,
            }
        if self.show_comparison_table:
            context["table_headers"] = self.table_headers
            context["table_rows"] = self.table_rows
            context["table_title"] = self.table_title
            context["table_help"] = self.table_help
            context["first_coloured_row"] = self.first_coloured_row
            context["n_data_columns"] = len(self.base_info_headers) + \
                len(self.targets)
            context["hist_colour_index_shift"] = self.hist_colour_index_shift
        return my_locals(context)

    @property
    def hash_to_taxon_dict(self):
        if self._hash_to_taxon_dict is None:
            self._hash_to_taxon_dict = self.db.get_genomes_description().description.to_dict()
        return self._hash_to_taxon_dict

    def set_table_data(self):
        self.targets = self.form.get_choices()
        self.n_selected = len(self.targets)

        taxon_list = [self.hash_to_taxon_dict[target_id]
                      for target_id in self.targets]
        self.table_headers = self.base_info_headers + taxon_list

        self.table_rows = self.get_table_rows()
        self.n_rows = len(self.table_rows)

    def get_table_rows(self):
        """This method should return an iteratable object with each iteration
        yielding a row of the table. Apart from the data, rows can contain
        values used to color the cells, these are stored in the rows after
        the data and will be matched to the corresponding data value with a
        shift in index (row[i] colored by row[i + hist_colour_index_shift]).
        """
        raise NotImplementedError

    @property
    def table_title(self):
        return "Number of {} present at least once in 1 of the {} selected "\
               "genomes: <strong>{}</strong>".format(self.compared_obj_name,
                                                     self.n_selected,
                                                     self.n_rows)

    @property
    def form_title(self):
        return "Compare the distribution of shared {}.".format(
            self.compared_obj_name)

    @property
    def form_help(self):
        return "Compare the size of the {} shared by selected genomes (targets).".format(
            self.compared_obj_name)

    @property
    def first_coloured_row(self):
        return len(self.base_info_headers)


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
        form = comp_metabo_form()
        return render(request, 'chlamdb/module_comp.html', my_locals(locals()))

    taxids = form.get_choices()

    if len(taxids) == 0:
        form = comp_metabo_form()
        return render(request, 'chlamdb/module_comp.html', my_locals(locals()))

    module_hits = db.get_ko_count_cat(taxon_ids=taxids)
    tipo = type(taxids)
    grouped_count = module_hits.groupby(["taxon_id", "module_id"]).sum()
    grouped_count = grouped_count.unstack(level="taxon_id")
    grouped_count.columns = grouped_count["count"].columns
    grouped_count.fillna(0, inplace=True, downcast="infer")
    taxids = grouped_count.columns
    genomes = db.get_genomes_description().description.to_dict()
    modules_info = db.get_modules_info(grouped_count.index.tolist(), as_pandas=True)
    all_infos = modules_info.set_index("module_id").join(grouped_count)

    header = ["Module", "Category", "Sub-category", "Description"]
    entries = []

    taxons = []
    for taxid in taxids:
        taxons.append(genomes[taxid])

    for module_id, data in all_infos.iterrows():
        line = [format_module(module_id, to_url=True), data["cat"], data.subcat, data.descr]
        for taxid in taxids:
            line.append(data[taxid])
        entries.append(line)
    envoi_comp = True
    return render(request, 'chlamdb/module_comp.html', my_locals(locals()))


class PfamComparisonView(TabularComparisonViewBase):

    view_type = "pfam"
    base_info_headers = ["Domain ID", "Description", "nDomain"]

    table_help = """
    The ouput table contains the list of shared Pfam domains and the number of
    times each of them was identified in the selected genomes.
    <br>nDomain: total number of occurence of this domain in the
    complete database.
    <br>Click on Pfam accession to get detailed phylogenetic profile of the
    corresponding Pfam entry.
    """

    compared_obj_name = "domains"

    def get_table_rows(self):
        pfam_hits = self.db.get_pfam_hits(ids=self.targets)
        pfam_defs = self.db.get_pfam_def(pfam_hits.index.tolist(),
                                         add_ttl_count=True)

        table_rows = []
        for key, values in pfam_hits.iterrows():
            entry_infos = pfam_defs.loc[key]
            base_infos = [format_pfam(key, to_url=True), entry_infos["def"],
                          entry_infos.ttl_cnt]
            table_rows.append(base_infos + values.values.tolist())

        return table_rows


class CogComparisonView(TabularComparisonViewBase):

    view_type = "cog"
    base_info_headers = ["COG accession", "Description", "# complete DB", "# genomes"]

    table_help = """
    The ouput table contains the list of COG annotated in selected genomes and
    the number of times each of them was identified in each genome.
    <br>Click on COG accession to get detailed phylogenetic profile of the
    corresponding COG entry.
    """

    compared_obj_name = "COG"

    def get_table_rows(self):
        cog_hits = self.db.get_cog_hits(
            ids=self.targets, search_on="taxid", indexing="taxid")
        # retrieve entry list
        cog_all = self.db.get_cog_hits(
            ids=list(self.hash_to_taxon_dict.keys()),
            search_on="taxid",
            indexing="taxid")

        # count frequency and n genomes
        cog_count = cog_all.sum(axis=1)
        cog_freq = cog_all[cog_all > 0].count(axis=1)
        cog_freq = cog_freq.loc[cog_hits.index]
        cog_count = cog_count.loc[cog_hits.index]
        # retrieve annotations
        cogs_summaries = self.db.get_cog_summaries(
            cog_hits.index.tolist(), as_df=True, only_cog_desc=True)

        cog2description = cogs_summaries.description.to_dict()
        cog_hits["accession"] = [format_cog(cog, as_url=True)
                                 for cog in cog_hits.index]
        cog_hits["description"] = [cog2description[cog] if cog in cog2description else '-'
                                   for cog in cog_hits.index]

        # combine into df
        combined_df = cog_hits.merge(
            cog_count.rename('count'),
            left_index=True,
            right_index=True).merge(
                cog_freq.rename('freq'),
                left_index=True,
                right_index=True).sort_values(["count"], ascending=False)

        cols = combined_df.columns.to_list()
        ordered_cols = cols[self.n_selected:] + cols[:self.n_selected]
        return combined_df[ordered_cols].values

    @property
    def first_coloured_row(self):
        return 4


class OrthogroupComparisonView(TabularComparisonViewBase):

    view_type = "orthology"
    base_info_headers = ["Orthogroup", "Annotaion"]

    table_help = """
    The ouput table contains the number of homologs in the shared orthogroups
    of the selected genomes. Interesting for comparing the size of orthogroups
    within genomes.
    <br> Homolog counts can be reordrered by clicking on column headers.<br>
    <br>Click on Orthologous group to get all the homologs identified in the
    database and the phylogenetic profile.
    """

    compared_obj_name = "orthogroups"

    @property
    def view_name(self):
        return "orthogroup_comparison"

    def get_table_rows(self):
        og_count = self.db.get_og_count(self.targets)
        annotations = self.db.get_genes_from_og(
            orthogroups=og_count.index.tolist(),
            taxon_ids=self.targets, terms=["product"])

        products = annotations.groupby("orthogroup")["product"].apply(list)

        og_data = []
        for og, items in og_count.iterrows():
            row = [format_orthogroup(og, to_url=True)]
            if og in products.index:
                row.append(format_lst_to_html(products.loc[og]))
            else:
                row.append("-")
            row.extend(items)

            og_data.append(row)

        return og_data


class KoComparisonView(TabularComparisonViewBase):

    view_type = "ko"

    table_help = """
    The ouput table contains the number of homologs in the shared Kegg
    Orthologs of the selected genomes and the total number of homologs
    of each Kegg Orthologs identified in the whole collection. Interesting
    for comparing the size of Kegg Orthologs within genomes.
    <br> Homolog counts can be reordrered by clicking on column headers.<br>
    <br>Click on the Ko entry and list the Ko modules and pathways of which it
    is part.
    """

    base_info_headers = ["KO", "Annot", "tot"]
    compared_obj_name = "KO"

    def get_table_rows(self):
        hits = self.db.get_ko_count(self.targets).unstack(level=0, fill_value=0)
        hits.columns = [col for col in hits["count"].columns.values]

        ko2annot = self.db.get_ko_desc(hits.index.tolist())
        df_ttl = self.db.get_ko_count(hits.index.tolist(), search_on="ko_id")
        ko2total_count = df_ttl.groupby("KO").sum()["count"].to_dict()
        table_rows = []
        for key, values in hits.iterrows():
            row = [format_ko(key, as_url=True), ko2annot[key], ko2total_count[key]]
            row.extend(values.values)
            table_rows.append(row)
        return table_rows


class AmrClassComparisonView(TabularComparisonViewBase):

    view_type = "amr"
    compared_obj_name = "AMR"

    base_info_headers = ["Class"]
    group_by = "class"

    _table_help = """
    The ouput table contains the number of times a given AMR {} appears
    in the selected genomes, color coded according to the quality
    (coverage*identity) of the best hit for that genome.<br>
    <br> Counts can be reordrered by clicking on column headers.<br>
    """

    @property
    def table_help(self):
        return self._table_help.format(self.group_by)

    def get_row_data(self, groupid, data):
        return [groupid]

    def get_table_rows(self):
        hits = self.db.get_amr_hits_from_taxonids(self.targets)

        table_rows = []
        hits["quality"] = hits["coverage"] * hits["identity"] / 10000
        for groupid, data in hits.groupby(self.group_by):
            row = self.get_row_data(groupid, data)
            taxonids = data["bioentry.taxon_id"]
            values = [len(taxonids[taxonids == target_id])
                      for target_id in self.targets]
            colours = [data[taxonids == target_id]["quality"].max() if value
                       else 0 for value, target_id in zip(values, self.targets)]
            row.extend(values)
            row.extend(colours)
            table_rows.append(row)
        return table_rows

    @property
    def hist_colour_index_shift(self):
        return len(self.targets)

    @property
    def view_name(self):
        return f"{self.view_type}_{self.group_by}_comparison"

    @property
    def tab_name(self):
        return f"{self.group_by}_comp"


class AmrGeneComparisonView(AmrClassComparisonView):

    base_info_headers = ["Gene", "scope", "Class", "Subclass", "Annotation"]
    group_by = "gene"

    _scope_hint = """
    <br> Note that genes are split into "core" and "plus" scopes, where
    "core" proteins are expected to have an effect on resistance while
    "plus" proteins are included with a less stringent criteria.<br>
    """

    @property
    def table_help(self):
        return self._table_help.format(self.group_by) + self._scope_hint

    def get_row_data(self, groupid, data):
        return [format_gene_to_ncbi_hmm((groupid, data.iloc[0].hmm_id)),
                data.iloc[0]["scope"],
                safe_replace(data.iloc[0]["class"], "/", " / "),
                safe_replace(data.iloc[0]["subclass"], "/", " / "),
                data.iloc[0]["seq_name"]]


class AmrSubclassComparisonView(AmrClassComparisonView):

    base_info_headers = ["Subclass", "Class"]
    group_by = "subclass"

    def get_row_data(self, groupid, data):
        return [groupid, safe_replace(data.iloc[0]["class"], "/", " / ")]


def faq(request):
    a = 2
    return render(request, 'chlamdb/FAQ.html', my_locals(locals()))


def phylogeny(request):
    page_title = page2title["phylogeny_intro"]

    biodb_path = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb_path)

    genomes_data = db.get_genomes_infos()
    genomes_descr = db.get_genomes_description()

    asset_path = "/temp/species_tree.svg"
    path = settings.BASE_DIR + '/assets/temp/species_tree.svg'

    core = db.get_n_orthogroups(only_core=True)

    genomes_data = genomes_data.join(genomes_descr)

    genomes_data.gc = genomes_data.gc.apply(round)
    genomes_data.coding_density = genomes_data.coding_density.apply(lambda x: round(100 * x))
    genomes_data.length = genomes_data.length.apply(lambda x: round(x / pow(10, 6), 2))

    data_table_header = ["Name", "%GC", "N proteins", "N contigs", "Size (Mbp)", "Percent coding"]
    data_table = genomes_data[["description", "gc", "n_prot", "n_contigs", "length", "coding_density"]].values.tolist()

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
        ["completeness", "Completeness", "#d7191c", False],
        ["contamination", "Contamination", "black", False]]

    for serie_name, header, col, is_relative in tree_params:
        data = genomes_data[serie_name]
        e_tree.add_column(SimpleColorColumn.fromSeries(data, header=None, use_col=False))
        stack = StackedBarColumn(data.to_dict(), colours=[col, "white"],
                                 relative=is_relative, header=header,
                                 header_params=header_params, face_params=stacked_face_params)
        e_tree.add_column(stack)

    e_tree.rename_leaves(genomes_descr.description.to_dict())
    e_tree.render(path, dpi=500)
    return render(request, 'chlamdb/phylogeny_intro.html', my_locals(locals()))


def genomes_intro(request):
    page_title = page2title["genomes_intro"]

    biodb_path = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb_path)

    genomes_data = db.get_genomes_infos()
    genomes_descr = db.get_genomes_description()

    asset_path = "/temp/species_tree.svg"
    path = settings.BASE_DIR + '/assets/temp/species_tree.svg'

    genomes_data = genomes_data.join(genomes_descr)

    genomes_data.gc = genomes_data.gc.apply(round)
    genomes_data.coding_density = genomes_data.coding_density.apply(lambda x: round(100 * x))
    genomes_data.length = genomes_data.length.apply(lambda x: round(x / pow(10, 6), 2))

    data_table_header = ["Name", "%GC", "N proteins", "N contigs", "Size (Mbp)", "Percent coding"]
    data_table = genomes_data[["description", "gc", "n_prot", "n_contigs", "length", "coding_density"]].values.tolist()

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
        ["completeness", "Completeness", "#d7191c", False],
        ["contamination", "Contamination", "black", False]]

    for serie_name, header, col, is_relative in tree_params:
        data = genomes_data[serie_name]
        e_tree.add_column(SimpleColorColumn.fromSeries(data, header=None, use_col=False))
        stack = StackedBarColumn(data.to_dict(), colours=[col, "white"],
                                 relative=is_relative, header=header,
                                 header_params=header_params, face_params=stacked_face_params)
        e_tree.add_column(stack)

    e_tree.rename_leaves(genomes_descr.description.to_dict())
    e_tree.render(path, dpi=500)
    return render(request, 'chlamdb/genomes_intro.html', my_locals(locals()))
