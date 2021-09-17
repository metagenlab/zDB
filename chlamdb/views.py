# -*- coding: utf-8 -*-

# todo circos gc file curently written in home directory, move it to other place
# todo save temp files in temp folder




from django.shortcuts import render
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


import seaborn as sns
import matplotlib.colors as mpl_col

import collections

import numpy 
import re
import os

from django.http import HttpResponseRedirect
from django.http import HttpResponse
from django.http import JsonResponse

from chlamdb.forms import make_contact_form
from chlamdb.forms import make_plot_form
from chlamdb.forms import SearchForm
from chlamdb.forms import BiodatabaseForm
from chlamdb.forms import make_circos_form
from chlamdb.forms import make_circos2genomes_form
from chlamdb.forms import make_mummer_form
from chlamdb.forms import make_blast_form
from chlamdb.forms import make_crossplot_form
from chlamdb.forms import ConnexionForm
from chlamdb.forms import DBForm
from chlamdb.forms import make_motif_form
from chlamdb.forms import PCRForm
from chlamdb.forms import make_extract_form
from chlamdb.forms import make_circos_orthology_form
from chlamdb.forms import make_interpro_from
from chlamdb.forms import make_metabo_from
from chlamdb.forms import make_module_overview_form
from chlamdb.forms import make_extract_region_form
from chlamdb.forms import make_venn_from
from chlamdb.forms import make_priam_form
from chlamdb.forms import AnnotForm
from chlamdb.forms import hmm_sets_form
from chlamdb.forms import hmm_sets_form_circos
from chlamdb.forms import make_blastnr_form
from chlamdb.forms import make_genome_selection_form
from chlamdb.forms import make_comment_from
from chlamdb.forms import LocusInt
from chlamdb.forms import get_LocusAnnotForm
from chlamdb.forms import BlastProfileForm
from chlamdb.forms import make_pairwiseid_form
from chlamdb.forms import make_locus2network_form
from chlamdb.forms import heatmap_form
from chlamdb.forms import blast_sets_form
from chlamdb.forms import make_pairwiseCDS_length_form
from django.contrib.auth import logout
from django.conf import settings
from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django.forms.utils import flatatt


from django.core.cache import cache
from tempfile import NamedTemporaryFile
from Bio import SeqIO
import simplejson
import string
import random
import json
from django.shortcuts import render
from celery.result import AsyncResult
from django.http import HttpResponse
from django.http import StreamingHttpResponse
from chlamdb.forms import GenerateRandomUserForm

from metagenlab_libs import db_utils
from metagenlab_libs.ete_phylo import EteTree, SimpleColorColumn, ModuleCompletenessColumn
from metagenlab_libs.ete_phylo import KOAndCompleteness
from metagenlab_libs.ete_phylo import Column
from metagenlab_libs.KO_module import ModuleParser

from metagenlab_libs.chlamdb import search_bar as sb

from ete3 import Tree

from reportlab.lib import colors

from ete3 import TextFace, StackedBarFace, SeqMotifFace

import pandas as pd

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram


# could also be extended to cache the results of frequent queries
# (e.g. taxid -> organism name) to avoid db queries
with db_utils.DB.load_db_from_name(settings.BIODB_DB_PATH) as db:
    hsh_config = db.get_config_table(ret_mandatory=True)
    optional2status = { name: value for name, (mandatory, value) in hsh_config.items() if not mandatory}
    missing_mandatory = [name for name, (mandatory, value) in hsh_config.items()
            if mandatory and not value]
    
def my_locals(local_dico):
    local_dico["optional2status"] = optional2status
    local_dico["missing_mandatory"] = missing_mandatory
    return local_dico


def get_task_info(request):
    
    task_id = request.GET.get('task_id', None)
    
    if task_id is not None:
        task = AsyncResult(task_id)
        data = {
            'state': task.state,
            'result': task.result,
        }
        
        return HttpResponse(json.dumps(data), content_type='application/json')
    else:
        return HttpResponse('No job id given.')


if settings.DEBUG == True:
    debug_mode = 1
else:
    debug_mode = 0





def id_generator(size=6, chars=string.ascii_uppercase + string.ascii_lowercase + string.digits):
   return ''.join(random.choice(chars) for _ in range(size))


def extract_alphanumeric(input_string):
    from string import ascii_letters, digits
    import string
    return "".join([ch for ch in input_string if ch in (ascii_letters + digits + '_-.')])



def help(request):
    return render(request, 'chlamdb/help.html', my_locals(locals()))


def about(request):
    import re
    import bibtexparser
    path = settings.BASE_DIR + '/assets/bibliography/references.bib'
    with open(path) as bibtex_file:
        bib_database = bibtexparser.load(bibtex_file)

    entry_list = []

    for entry in bib_database.entries:

        """
        Sophie S Abby, Jean Cury, Julien Guglielmini, Bertrand Néron, Marie Touchon, and Eduardo PC
        Rocha. Identification of protein secretion systems in bacterial genomes. Scientific reports, 6, 2016.
        53, 57, 61, 164

        Peter JA Cock, Tiago  Antao, Jeffrey T Chang, Brad
        Biopython: freely    available    python    tools    for computational molecular biology and bioinformatics.Bioinformatics, 25(11):1422–1423, 2009.

        {'publisher': 'Oxford University Press', 'year': '2009', 'pages': '1422--1423', 'number': '11', 'volume': '25',
         'journal': 'Bioinformatics',
         'author': 'Cock, Peter JA and Antao, Tiago and Chang, Jeffrey T and Chapman, Brad A and Cox, Cymon J and Dalke, Andrew and Friedberg, Iddo and Hamelryck, Thomas and Kauff, Frank and Wilczynski, Bartek and others',
         'title': 'Biopython: freely available Python tools for computational molecular biology and bioinformatics',
         'ENTRYTYPE': 'article', 'ID': 'cock2009biopython'}

        """
        string = ("<b>%s</b></br> %s, %s, %s(%s):%s, %s" % (re.sub('[{}]','', entry["title"]),
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

def create_user(username, mail, password, first_name, last_name, staff=True):
    from django.contrib.auth.models import User
    user = User.objects.create_user(username, mail, password)
    user.first_name, user.last_name = first_name, last_name
    user.is_staff = staff
    user.save()


class StackedBarColumn(Column):
    def __init__(self, values, header, colours=None, relative=False,
            face_params=None, header_params=None):
        super().__init__(header, face_params, header_params)
        self.values=values
        self.colours = colours
        self.relative = relative
        if relative:
            self.max = max(values.values())
            self.min = min(values.values())

    def get_face(self, index):
        val = self.values[int(index)]

        if self.relative and self.max!=self.min:
            val = 100*float(val-self.min)/(self.max-self.min)
        elif self.relative:
            val = 100

        face = StackedBarFace([val, 100-val], width=50, height=9, colors=self.colours)
        self.set_default_params(face)
        face.inner_border.color = "black"
        face.inner_border.width = 0
        return face


def home(request):
    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb_path)

    genomes_data = db.get_genomes_infos()
    genomes_descr = db.get_genomes_description()

    asset_path = "/temp/species_tree.svg"
    path = settings.BASE_DIR + '/assets/temp/species_tree.svg'

    genomes_data = genomes_data.join(genomes_descr)

    genomes_data.gc = genomes_data.gc.apply(round)
    genomes_data.coding_density = genomes_data.coding_density.apply(lambda x: round(100*x))
    genomes_data.length = genomes_data.length.apply(lambda x: round(x/pow(10,6), 2))

    data_table_header = ["Name", "%GC", "N proteins", "N contigs", "Size (Mbp)", "Percent coding"]
    data_table = genomes_data[["description", "gc", "n_prot", "n_contigs", "length", "coding_density"]].values.tolist()

    # plot phylo only of not already in assets
    tree = db.get_reference_phylogeny()
    t1 = Tree(tree)
    R = t1.get_midpoint_outgroup()
    if not R is None:
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
        stack = StackedBarColumn(data.to_dict(), colours = [col, "white"],
                relative=is_relative, header=header,
                header_params=header_params, face_params=stacked_face_params)
        e_tree.add_column(stack)

    e_tree.rename_leaves(genomes_descr.description.to_dict())
    e_tree.render(path, dpi=500)
    return render(request, 'chlamdb/home.html', my_locals(locals()))


def format_lst_to_html(lst_elem, add_count=True, format_func=None):
    if format_func==None:
        format_func = lambda x:x

    dict_elem = {}
    for elem in lst_elem:
        if pd.isna(elem):
            elem = "-"
        cnt = dict_elem.get(elem, 0)
        dict_elem[elem] = cnt+1

    elems = []
    for k, v in dict_elem.items():
        if k != "-":
            token = format_func(k)
        else:
            token = k
        if add_count and k!="-":
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

    if len(annotations)==2:
        return header, annotations[0].join(annotations[1], how="outer")
    elif len(annotations)==1:
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
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

    extract_form_class = make_extract_form(db, plasmid=True)
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
    if not include_plasmids is None:
        sum_include_lengths += len(include_plasmids)

    if n_missing>=sum_include_lengths:
        wrong_n_missing = True
        return render(request, 'chlamdb/extract_orthogroup.html', my_locals(locals()))

    og_counts_in = db.get_og_count(include_taxids, plasmids=include_plasmids)
    if not single_copy:
        og_counts_in["presence"] = og_counts_in[og_counts_in > 0].count(axis=1)
    else:
        og_counts_in["presence"] = og_counts_in[og_counts_in == 1].count(axis=1)

    og_counts_in["selection"] = og_counts_in.presence >= (sum_include_lengths-n_missing)

    sum_exclude_lengths = len(exclude_taxids)
    if not exclude_plasmids is None:
        sum_exclude_lengths += len(exclude_plasmids)
    if sum_exclude_lengths>0:
        mat_exclude = db.get_og_count(exclude_taxids, plasmids=exclude_plasmids)
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

    if not single_copy:
        orthogroup2count_all = count_all_genomes[count_all_genomes > 0].count(axis=1)
    else:
        orthogroup2count_all = count_all_genomes[count_all_genomes == 1].count(axis=1)
    max_n = orthogroup2count_all.max()
    match_groups_data = []

    all_taxids = include_taxids
    if not include_plasmids is None:
        all_taxids += include_plasmids
    annotations = db.get_genes_from_og(orthogroups=selection, taxon_ids=all_taxids,
        terms=["gene", "product", "locus_tag"])
    if annotations.empty:
        no_match = True
        return render(request, 'chlamdb/extract_orthogroup.html', my_locals(locals()))
    

    opt_header, optional_annotations = get_optional_annotations(db, seqids=annotations.index.tolist())
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
    table_headers.extend([f"Present in {sum_include_lengths}", f"Freq complete database (/{max_n})"])

    for row, count in orthogroup2count_all.iteritems():
        cnt_in = og_counts_in.presence.loc[row]
        g = genes.loc[row]
        gene_data = format_lst_to_html("-" if pd.isna(entry) else entry for entry in g)
        prod_data = format_lst_to_html(products.loc[row])
        column_header = format_orthogroup(row)
        optional = []
        if "KO" in opt_header and row in kos:
            optional.append(format_lst_to_html(kos.loc[row], add_count=True, format_func=format_ko_url))
        if "COG" in opt_header and row in cogs:
            optional.append(format_lst_to_html(cogs.loc[row], add_count=True, format_func=format_cog_url))
        entry = [column_header, gene_data, prod_data, *optional, cnt_in, count]
        match_groups_data.append(entry)

    ref_genomes = db.get_genomes_description().loc[include_taxids].reset_index()

    envoi_extract = True
    return render(request, 'chlamdb/extract_orthogroup.html', my_locals(locals()))


def venn_orthogroup(request):
    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb_path)

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
        ogs_str = ",".join(f"{to_s(format_orthogroup(og))}" for og, cnt in ogs.iteritems() if cnt > 0)
        genome = genomes.loc[int(taxon)].description
        fmt_data.append(f"{{name: {to_s(genome)}, data: [{ogs_str}] }}")
    series = "[" + ",".join(fmt_data) + "]"

    og_list = og_count.index.tolist()
    annotations = db.get_genes_from_og(orthogroups=og_list, taxon_ids=genomes.index.tolist())
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
        og_item = f"h[{to_s(format_orthogroup(og))}] = {og_info};"
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


def extract_pfam(request, classification="taxon_id"):
    db = db_utils.DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    extract_form_class = make_extract_form(db, plasmid=True, label="Pfam domains")

    if request.method != "POST":
        form = extract_form_class()
        return render(request, 'chlamdb/extract_Pfam.html', my_locals(locals()))

    form = extract_form_class(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/extract_Pfam.html', my_locals(locals()))

    include, include_plasmids = form.get_include_choices()
    exclude, exclude_plasmids = form.get_exclude_choices()
    n_missing = form.get_n_missing()

    sum_include_length = len(include)
    if not include_plasmids is None:
        sum_include_length += len(include_plasmids)

    sum_exclude_length = len(exclude)
    if not exclude_plasmids is None:
        sum_exclude_length += len(exclude_plasmids)

    if n_missing>=sum_include_length:
        ctx = {"wrong_n_missing" : True}
        return render(request, 'chlamdb/extract_Pfam.html', my_locals(ctx))

    pfam_include = db.get_pfam_hits(include, plasmids=include_plasmids, 
            search_on="taxid", indexing="taxid")
    if sum_exclude_length > 0:
        pfam_exclude = db.get_pfam_hits(exclude, plasmids=exclude_plasmids,
                search_on="taxid", indexing="taxid")
        pfam_exclude["sum_pos"] = pfam_exclude[cog_exclude > 0].count(axis=1)
        pfam_exclude["exclude"] = pfam_exclude.sum_pos > 0
        neg_index = pfam_exclude[pfam_exclude.exclude].index
    else:
        neg_index = pd.Index([])

    pfam_include["sum_pos"] = pfam_include[pfam_include > 0].count(axis=1)
    pfam_include["selection"] = pfam_include.sum_pos >= len(include)-n_missing
    pos_index = pfam_include[pfam_include.selection].index
    selection = pos_index.difference(neg_index).tolist()
    pfam_defs = db.get_pfam_def(selection)

    if len(selection) == 0:
        ctx = {no_match: True}
        return render(request, 'chlamdb/extract_ko.html', my_locals(ctx))

    all_database = db.get_pfam_hits(pfam_include.index.tolist(), search_on="pfam", indexing="taxid")
    sums = all_database.sum(axis=1)
    sum_group = len(selection)

    match_groups_data = []
    for no, pfam in enumerate(selection):
        count = sums.loc[pfam]
        pfam_def = pfam_defs["def"].loc[pfam]
        data = [no+1, format_pfam(pfam), pfam_def,
                pfam_include.sum_pos.loc[pfam], sums.loc[pfam]]
        match_groups_data.append(data)

    ctx = {"envoi_extract": True,
            "sum_group": sum_group,
            "n_genomes": sum_include_length,
            "max_n": sums.max(),
            "match_groups_data": match_groups_data,
            "form": form}
    return render(request, 'chlamdb/extract_Pfam.html', my_locals(ctx))


def format_ko(ko_id, as_url=False, base=None):
    if base is None:
        base = f"K{int(ko_id):05d}"
    if not as_url:
        return base
    return f"<a href=\"/fam_ko/{base}\">{base}</a>"

def format_ko_url(ko_id):
    return format_ko(ko_id, as_url=True)


def format_ko_path(hsh_pathways, ko, as_list=False):
    pathways = hsh_pathways.get(ko, [])
    if len(pathways) == 0:
        if as_list:
            return []
        return "-"
    fmt_lst = (f"<a href=\"/KEGG_mapp_ko/map{i:05d}\">{d}</a>" for i, d in pathways)
    
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
    db = db_utils.DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    extract_form_class = make_extract_form(db, plasmid=True, label="Kegg Orthologs")

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
    if not include_plasmids is None:
        sum_include_length += len(include_plasmids)

    sum_exclude_length = len(exclude)
    if not exclude_plasmids is None:
        sum_exclude_length += len(exclude_plasmids)

    if n_missing>=sum_include_length:
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
    mat_include["selection"] = mat_include.sum_pos >= len(include)-n_missing
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

    # should contain the loci of the KO in the reference genome
    max_n = ko_total_count.max()
    locus_list = [] # [i[0] for i in server.adaptor.execute_and_fetchall(locus_list_sql,)]
    target_circos_taxons = include + exclude

    # url to get the barchart of selected KO
    taxons_in_url = "?i="+("&i=").join(map(str, include)) + '&m=%s' % str(n_missing)
    taxon_out_url = "&o="+("&o=").join(map(str, exclude))
    envoi_extract = True
    mm = 'module'
    pp = 'pathway'

    locus2annot = {}
    locus_tag2cog_catego = {}
    locus_tag2cog_name = {}
    locus_tag2ko = {}
    pathway2category = {}
    module2category = {}
    ko2ko_pathways = {}
    ko2ko_modules = {}
    locus2interpro = {}
    return render(request, 'chlamdb/extract_ko.html', my_locals(locals()))



def venn_pfam(request):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

    venn_form_class = make_venn_from(db, limit=6)
    if request.method != "POST":
        form_venn = venn_form_class()
        return render(request, 'chlamdb/venn_Pfam.html', my_locals({"form_venn": form_venn}))

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
        str_fmt = ",".join(f"\"{format_pfam(pfam)}\"" for pfam, _ in non_zero.iteritems())
        series_tab.append(f"{{name: \"{genomes_desc[target]}\", data: [{str_fmt}]}}")
    series = "[" + ",".join(series_tab) + "]"

    descriptions = []
    for pfam, pfam_info in data.iterrows():
        pfam_def = pfam_info["def"]
        descriptions.append(f"h[\"{format_pfam(pfam)}\"] = \"{pfam_def}\"")

    ctx = {"envoi_venn": True,
            "series": series,
            "pfam2description": ";".join(descriptions)}
    return render(request, 'chlamdb/venn_Pfam.html', my_locals(ctx))


def extract_cog(request):
    db = db_utils.DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    extract_form_class = make_extract_form(db, plasmid=True, label="COG")

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
    if not include_plasmids is None:
        sum_include_length += len(include_plasmids)

    sum_exclude_length = len(exclude)
    if not exclude_plasmids is None:
        sum_exclude_length += len(exclude_plasmids)

    if n_missing>=sum_include_length:
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
    cog_include["selection"] = cog_include.sum_pos >= len(include)-n_missing
    pos_index = cog_include[cog_include.selection].index
    selection = pos_index.difference(neg_index).tolist()
    if len(selection) == 0:
        no_match = True
        return render(request, 'chlamdb/extract_ko.html', my_locals(locals()))

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
            cat_count[func] = (inc+cog_include.sum_pos.loc[cog_id], not_incl)
        funcs = "<br>".join(f"{func} ({func_desc})" for func, func_desc in func_acc)
        data = (format_cog(cog_id), funcs, cog_descr,
                cog_include.sum_pos.loc[cog_id], str(count))
        cog_data.append(data)

    # get the categories for all cogs
    for cog_id, details_lst in cogs_summaries.items():
        for func, func_descr, cog_descr in details_lst:
            inc, not_incl = cat_count.get(func, (0, 0))
            cat_count[func] = (inc, not_incl+sums.loc[cog_id])

    max_n = sums.max(axis=0)
    sum_group = len(selection)

    # Code to generate the barchart diagrams
    cat_map_str = ",".join([f"\"{func}\": \"{descr}\"" for func, descr in cogs_funct.items()])
    category_map = f"var category_description = {{ {cat_map_str} }};"
    ttl_sel = sum([c1 for func, (c1, c2) in cat_count.items()])
    ttl_all = sum([c2 for func, (c1, c2) in cat_count.items()])

    cat_count_comp = ",".join([f"\"{func}\" : [\"{c2}\", \"{c1}\"]" for func, (c1, c2) in cat_count.items()])
    category_count_complete = f"var category_count_complete = {{ {cat_count_comp} }};"

    serie_selection_val = [str(round(float(c1)/ttl_sel, 2)) for func, (c1, c2) in cat_count.items()]
    serie_all_val = [str(round(float(c2)/ttl_all, 2)) for func, (c1, c2) in cat_count.items()]
    serie_selection = f"{{ labels: \"selection\", values: [{','.join(serie_selection_val)}] }}"
    serie_all = f"{{ labels: \"complete genomes\", values: [{','.join(serie_all_val)}] }}"
    series_str = ",".join([serie_all, serie_selection])
    series = f"[ {series_str} ]"
    labels_str = ",".join([f"\"{funct}\"" for funct in cat_count.keys()])
    labels = f"[ {labels_str} ]"

    locus2annot = {}
    locus_tag2cog_catego = {}
    locus_tag2cog_name = {}
    locus_tag2ko = {}
    pathway2category = {}
    module2category = {}
    ko2ko_pathways = {}
    ko2ko_modules = {}
    locus2interpro = {} #  get_locus_annotations(biodb, locus_list)
    envoi_extract = True
    return render(request, 'chlamdb/extract_cogs.html', my_locals(locals()))


def venn_ko(request):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

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

    taxids    = form_venn.get_taxids()
    genomes   = db.get_genomes_description().description.to_dict()
    ko_counts = db.get_ko_count(taxids)

    fmt_data = []
    ko_list = ko_counts.index.get_level_values("KO").unique().to_list()
    for taxid in taxids:
        kos = ko_counts.loc[taxid].index.values
        kos_str = ",".join(f"{to_s(format_ko(ko))}" for ko in kos)
        genome = genomes[taxid]
        fmt_data.append(f"{{name: {to_s(genome)}, data: [{kos_str}] }}")
    series = "[" + ",".join(fmt_data) + "]"

    ko_descriptions = db.get_ko_desc(ko_list)
    ko2description = []
    for ko, ko_desc in ko_descriptions.items():
        forbidden = "\""
        ko_item = f"h[{to_s(format_ko(ko))}] = [{forbidden}{ko_desc}{forbidden}];"
        ko2description.append(ko_item)
    ko2description = "\n".join(ko2description)
    envoi_venn = True
    return render(request, 'chlamdb/venn_ko.html', my_locals(locals()))


def format_cog(cog_id, as_url=False, base=None):
    if base is None:
        base = f"COG{int(cog_id):04d}"
    if as_url==False:
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
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

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
        str_fmt = ",".join(f"\"{format_cog(cog)}\"" for cog, count in non_zero_cogs.iteritems())
        series_tab.append( f"{{name: \"{genome_desc[target]}\", data: [{str_fmt}]}}" )
    series = "[" + ",".join(series_tab) + "]"

    cog2description_l = []
    cog_codes = db.get_cog_code_description()
    for cog, data in data.iterrows():
        name = escape_quotes(data.description)
        func = data.function
        functions = ",".join(f"\"{abbr}\"" for abbr in func)
        cog2description_l.append(f"h[\"{format_cog(cog)}\"] = [[{functions}], \"{name}\"]")

    cog_func_dict = (f"\"{func}\": \"{descr}\"" for func, descr in cog_codes.items())
    cog_func_dict = "{"+",".join(cog_func_dict)+"}"
    cog2description = ";".join(cog2description_l)
    envoi_venn = True
    return render(request, 'chlamdb/venn_cogs.html', my_locals(locals()))



def format_lst(lst):
    hsh_values = {}
    for item in lst:
        val = hsh_values.get(item, 0)
        hsh_values[item] = val+1
    return hsh_values


def tab_homologs(db, infos, hsh_organism, ref_seqid=None, og=None):
    n_genomes = len(set(hsh_organism.values()))
    if n_genomes==1:
        n_genomes = "1 genome"
    else:
        n_genomes = f"{n_genomes} genomes"

    headers = ["", "Locus tag", "Source", "Gene", "Product"]
    identities = None
    if ref_seqid!=None:
        identities = db.get_og_identity(og, ref_seqid)
        headers.insert(2 , "Identity")

    homologues = []
    index = 0
    for seqid, data in infos.iterrows():
        organism = hsh_organism[seqid]
        locus_fmt = format_locus(data.locus_tag, to_url=True)
        entry = [index+1, locus_fmt, organism, format_gene(data.gene), data["product"]]
        if ref_seqid!=None:
            if seqid==ref_seqid:
                ident = "N/A"
            else:
                ident = round(identities.loc[seqid].identity, 1)
            if ident==0:
                ident = "-"
            entry.insert(2, ident)

        homologues.append(entry)
        index += 1

    return {"orthogroup": orthogroup,
            "n_genomes": n_genomes,
            "headers": headers,
            "homologues": homologues }


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
            print(existing_id, div_id)
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
        return { "length_distrib": True, "single_length": True, "prot_length": lengths.iloc[0]}

    return {"length_distrib": True,
            "max_protein_length": max_protein_length,
            "std_protein_length": std_protein_length,
            "min_protein_length": min_protein_length,
            "mean_protein_length": mean_protein_length,
            "median_protein_length": median_protein_length,
            "html_plot_prot_length": html_plot_prot_length }


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
        dummy_seq = "-"*prot_length
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


def tab_og_phylogeny(db, og_id):
    og_phylogeny = db.get_og_phylogeny(og_id)
    pfam_col = None
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

    tree = Tree(og_phylogeny)
    locuses = [branch.name for branch in tree.iter_leaves()]
    locus_to_genome = db.get_locus_to_genomes(locuses)
    R = tree.get_midpoint_outgroup()
    root = "(unrooted)"
    if not R is None:
        root = "(midpoint rooted)"
        tree.set_outgroup(R)
    tree.ladderize()
    e_tree = EteTree(tree)

    e_tree.add_column(SimpleTextColumn("Locus tag"))
    if not pfam_col is None:
        e_tree.add_column(pfam_col)
    e_tree.rename_leaves(locus_to_genome, leaf_name_type=str)

    asset_path = f"/temp/og_phylogeny{og_id}.svg"
    path = settings.BASE_DIR + '/assets/' + asset_path
    e_tree.render(path, dpi=1200)
    return {"og_phylogeny": asset_path, "root": root}


def tab_og_conservation_tree(db, group, compare_to=None):
    from metagenlab_libs.ete_phylo import EteTree, SimpleColorColumn, ModuleCompletenessColumn
    from ete3 import Tree
    
    ref_phylogeny = db.get_reference_phylogeny()
    leaf_to_name = db.get_genomes_description().description.to_dict()
    
    # Note: as the orthogroup may either be in a plasmid or in the chromosome
    # of the bacteria, we need to index by taxon to group them (index on taxon)
    count = db.get_og_count([group], search_on="orthogroup")

    tree = Tree(ref_phylogeny)
    R = tree.get_midpoint_outgroup()
    if not R is None:
        tree.set_outgroup(R)
    tree.ladderize()
    e_tree = EteTree(tree)

    e_tree.add_column(SimpleColorColumn.fromSeries(count.loc[group], header=format_orthogroup(group)))
    if not compare_to is None:
        identity_matrix = db.get_og_identity(group, compare_to)
        seqids = identity_matrix.index.tolist()

        # to get the taxid of the reference seqid, so as to exclude it from 
        # the phylogenetic tree
        seqids.append(compare_to)
        seqid_to_taxon  = db.get_taxid_from_seqid(seqids)
        identity_matrix["taxid"] = identity_matrix.index.map(seqid_to_taxon)
        max_identity = identity_matrix.groupby("taxid").max()
        col = LocusHeatmapColumn.fromSeries(max_identity["identity"],
                header="Identity", cls=LocusHeatmapColumn)
        col.ref_taxon = seqid_to_taxon[compare_to]
        e_tree.add_column(col)

    e_tree.rename_leaves(leaf_to_name)

    dpi = 1200
    asset_path = f"/temp/og_conservation{group}.svg"
    path = settings.BASE_DIR + '/assets/' + asset_path
    e_tree.render(path, dpi=dpi)
    return {"asset_path": asset_path}


def og_tab_get_kegg_annot(db, seqids):
    ko_hits = db.get_ko_hits(seqids, search_on="seqid", indexing="seqid")
    if ko_hits.empty:
        return {}

    n_occurences = ko_hits["ko"].value_counts()
    ko_ids = n_occurences.index.tolist()
    ko_descr = db.get_ko_desc(ko_ids)
    ko_pathways = db.get_ko_pathways(ko_ids)
    ko_modules = db.get_ko_modules(ko_ids)
        
    ko_entries = []
    for ko_id, count in n_occurences.iteritems():
        entry = [format_ko(ko_id, as_url=True), count]
        descr = ko_descr.get(ko_id, "-")
        entry.append(descr)
        entry.append(format_ko_path(ko_pathways, ko_id))
        entry.append(format_ko_modules(ko_modules, ko_id))
        ko_entries.append(entry)

    ko_header = ["KO", "Occurences", "Description", "Pathways", "Modules"]
    return {
        "ko_entries": ko_entries,
        "ko_header": ko_header
    }


def og_tab_get_cog_annot(db, seqids):
    cog_hits = db.get_cog_hits(seqids, indexing="seqid", search_on="seqid")

    if cog_hits.empty:
        return {}

    n_entries = cog_hits["cog"].value_counts()
    cog_summ = db.get_cog_summaries(n_entries.index.tolist())
    cog_entries = []
    for cog_id, count in n_entries.iteritems():
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


def orthogroup(request, og):
    tokens = og.split("_")
    try:
        og_id = int(tokens[1])
    except:
        menu = True
        invalid_id = True
        return render(request, "chlamdb/og.html", my_locals(locals()))
    else:
        valid_id = True

    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

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
        gene_annotations.append([index+1, gene, cnt])

    product_annotations = []
    for index, values in enumerate(hsh_products.items()):
        product, cnt = values
        if pd.isna(product):
            product = "-"
        product_annotations.append([index+1, product, cnt])

    swissprot, cog_ctx, kegg_ctx, pfam_ctx = {}, {}, {}, {}
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
    except:
        og_phylogeny_ctx = {}

    og_conserv_ctx = tab_og_conservation_tree(db, og_id)
    length_tab_ctx = tab_lengths(n_homologues, annotations)
    homolog_tab_ctx = tab_homologs(db, annotations, hsh_organisms)
    context = {
        "valid_id": valid_id,
        "n_homologues": n_homologues,
        "og": og,
        "menu": True,
        "gene_annotations": gene_annotations,
        "product_annotations": product_annotations,
        **homolog_tab_ctx,
        **length_tab_ctx,
        **og_conserv_ctx,
        **cog_ctx, **kegg_ctx, **pfam_ctx, **og_phylogeny_ctx, **swissprot
    }
    return render(request, "chlamdb/og.html", my_locals(context))


def tab_general(seqid, hsh_organism, gene_loc, annot):
    organism = hsh_organism[seqid]
    strand, beg, end = gene_loc[seqid]
    gene = annot.loc[seqid].gene
    if pd.isna(gene):
        gene = "-"
    length = annot.loc[seqid].length+1
    product = annot.loc[seqid]["product"]
    locus_tag = annot.loc[seqid].locus_tag
    return {
        "locus_tag": locus_tag,
        "organism": organism,
        "strand": strand,
        "gene": gene,
        "start": beg,
        "end": end,
        "nucl_length": end-beg,
        "length": length,
        "prot": product
    }


# to be moved somewhere else at some point
def to_color_code(c):
    red = int(256*c.red)
    green = int(256*c.green)
    blue = int(256*c.blue)
    return f"#{red:x}{green:x}{blue:x}"


class LocusHeatmapColumn(SimpleColorColumn):
    def __init__(self, values, ref_taxon=None, header=None):
        super().__init__(values, header)
        self.ref_taxon = ref_taxon
        if len(values) > 0:
            self.min_val = min(v for k, v in values.items())
            self.max_val = max(v for k, v in values.items())


    def get_face(self, index):
        index = int(index)
        if index==self.ref_taxon:
            text_face = TextFace("-".center(11))
            text_face.inner_background.color = EteTree.GREEN
            self.set_default_params(text_face)
            return text_face

        val = self.values.get(index, None)
        if val is None:
            return TextFace("-")

        color = colors.linearlyInterpolatedColor(colors.gray,
                colors.firebrick, self.min_val, self.max_val, val)
        text_face = TextFace(str(int(val)).center(12-len(str(int(val)))))
        text_face.inner_background.color = to_color_code(color)
        self.set_default_params(text_face)
        return text_face


def get_sequence(db, seqid, flanking=0):
    loc      = db.get_gene_loc([seqid])
    bioentry, accession, length, seq = db.get_bioentry_list(seqid, search_on="seqid")
    strand, start, stop = loc[seqid]
    start -= 1

    if start < 50:
        start_w_flank = 0
        red_start = start
    else:
        start_w_flank = start-flanking
        red_start = 50

    if stop+flanking > len(seq):
        stop_w_flank = len(seq)-1
    else:
        stop_w_flank = stop+flanking
    red_stop = red_start + stop-start
    fet = SeqFeature(FeatureLocation(start_w_flank, stop_w_flank, strand=strand))
    extracted = fet.extract(seq)
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
    for pfam, starts in pfam_starts.iteritems():
        ends = pfam_ends.loc[pfam]
        name = format_pfam(pfam)
        data = "["+",".join(f"{{x:{start}, y:{end}}}" for start, end in zip(starts, ends))+"]"
        feature = (
            f"{{ data: {data}, "
            f" name: \"{name}\","
            "  color: \"#0F8292\","
            "  type : \"rect\","
            "}"
        )
        pfam_def = pfam_defs_df["def"].loc[pfam]
        pfam_defs.append((format_pfam(pfam, to_url=True), pfam_def))
        feature_viewer_fet.append(feature)
    return {"pfam_domains": "[" + ",".join(feature_viewer_fet) + "]",
            "pfam_def": pfam_defs}


def locusx_genomic_region(db, seqid, window):
    filename = f"genomic_region_{window}_{seqid}.svg"
    hsh_loc = db.get_gene_loc([seqid])
    strand, start, end = hsh_loc[seqid]
    window_start, window_stop = start-window, end+window

    if window_start<0:
        window_start=0
    df_seqids = db.get_seqid_in_neighborhood(seqid, window_start, window_stop)
    bioentry, _, _, _ = db.get_bioentry_list(seqid, search_on="seqid")
    contig_size = db.get_contig_size(bioentry)

    hsh_organism = db.get_organism([seqid], id_type="seqid")
    infos = db.get_proteins_info(df_seqids.index.tolist(),
            to_return=["gene", "locus_tag"], as_df=True, inc_non_CDS=True, inc_pseudo=True)
    cds_type = db.get_CDS_type(df_seqids.index.tolist())
    all_infos = infos.join(cds_type).join(df_seqids)

    features = []
    if window_stop > contig_size:
        window_stop = contig_size
    genome_range = [window_start, window_stop]

    for curr_seqid, data in all_infos.iterrows():
        feature_name = ""
        if "gene" in data and not pd.isna(data.gene):
            feature_name = data.gene
        feature_type = data.type
        if data.is_pseudo:
            feature_type = "pseudo"
        features.append((
            f"{{start: {data.start}, gene: {to_s(feature_name)}, end: {data.end},"
            f"strand: {data.strand}, type: {to_s(feature_type)},"
            f"locus_tag: {to_s(data.locus_tag)}}}"
        ))
    return {"genome_range": genome_range, "window_size": 2*window,
            "features": "[" + ",".join(features)+"]"}


def locusx_RNA(db, seqid, is_pseudogene):
    hsh_infos = db.get_proteins_info([seqid], to_return=["locus_tag", "gene", "product"],
            inc_non_CDS=True, inc_pseudo=True)
    hsh_loc = db.get_gene_loc([seqid])
    organism = db.get_organism([seqid])

    locus_tag, gene, product = hsh_infos[seqid]
    if gene is None:
        gene = "-"
    strand, start, stop = hsh_loc[seqid]
    return {
        "locus_tag": locus_tag,
        "gene": gene,
        "prot": product,
        "start": start,
        "end": stop,
        "strand": strand,
        "nucl_length": stop-start,
        "organism": organism[seqid]
    }


def tab_get_refseq_homologs(db, seqid):
    refseq_hits = db.get_refseq_hits([seqid]).set_index("match_id")
    refseq_hits_infos = db.get_refseq_matches_info(refseq_hits.index.tolist())
    all_infos = refseq_hits.join(refseq_hits_infos)

    header = ["Refseq accession", "Evalue", "Score", "ID(%)", "# gaps", "Len", "Description", "Organism"]
    entries = []
    taxids = db.get_refseq_taxonomy(all_infos.taxid.tolist())

    for match_id, data in all_infos.iterrows():
        if data.taxid not in taxids.index:
            taxonomic_name = "-"
        else:
            taxonomic_name = taxids.loc[data.taxid].taxonomic_name

        to_ncbi = f"<a href=\"http://www.ncbi.nlm.nih.gov/protein/{data.accession}\">{data.accession}</a>"
        to_taxo = f"<a href=\"http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={data.taxid}\">{taxonomic_name}</a>"
        entries.append((to_ncbi, data.evalue, data.bitscore,
            data.pident, data.gaps, data.length, data.description, to_taxo))
    return { "n_refseq_homologs": len(refseq_hits),
            "refseq_headers": header,
            "blast_data": entries }


def locusx(request, locus=None, menu=True):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)
    
    if locus is None:
        return render(request, 'chlamdb/locus.html', my_locals({"valid_id": False}))
    try:
        seqid, feature_type, is_pseudogene = db.get_seqid(locus_tag=locus,
                feature_type=True)
    except:
        return render(request, 'chlamdb/locus.html', my_locals({"valid": False}))
    else:
        valid_id = True

    sequence = get_sequence(db, seqid, flanking=50)
    genomic_region_ctx = locusx_genomic_region(db, seqid, window=8000)
    if feature_type!="CDS" or is_pseudogene:
        ctx_RNA = locusx_RNA(db, seqid, is_pseudogene)
        if is_pseudogene:
            feature_type = "Pseudogene"
        context = {
            "valid_id": valid_id,
            "menu": True,
            "seq": sequence,
            "feature_type": feature_type,
            **ctx_RNA,
            **genomic_region_ctx
        }
        return render(request, 'chlamdb/locus.html', my_locals(context))

    gene_loc = db.get_gene_loc([seqid])
    og_inf   = db.get_og_count([seqid], search_on="seqid")
    og_id    = int(og_inf.loc[seqid].orthogroup) # need to convert from numpy64 to int
    og_annot = db.get_genes_from_og(orthogroups=[og_id],
            terms=["locus_tag", "gene", "product", "length"])
    all_og_c = db.get_og_count([og_id], search_on="orthogroup")
    all_org  = db.get_organism(og_annot.index.tolist(), as_hash=True)
    n_homologues = all_og_c.loc[og_id].sum()
    translation = db.get_translation(seqid)
    homolog_tab_ctx = tab_homologs(db, og_annot, all_org, seqid, og_id)
    general_tab     = tab_general(seqid, all_org, gene_loc, og_annot)
    og_conserv_ctx  = tab_og_conservation_tree(db, og_id, compare_to=seqid)

    kegg_ctx, cog_ctx, pfam_ctx = {}, {}, {}
    diamond_matches_ctx = {}
    n_swissprot_hits = 0
    if optional2status.get("KEGG", False):
        kegg_ctx = og_tab_get_kegg_annot(db, [seqid])

    if optional2status.get("COG", False):
        cog_ctx = og_tab_get_cog_annot(db, [seqid])

    if optional2status.get("pfam", False):
        pfam_ctx = tab_get_pfam_annot(db, [seqid])

    if optional2status.get("BLAST_swissprot", False):
        n_swissprot_hits = db.get_n_swissprot_homologs(seqid)

    if optional2status.get("BLAST_database", False):
        diamond_matches_ctx = tab_get_refseq_homologs(db, seqid)

    try:
        og_phylogeny_ctx = tab_og_phylogeny(db, og_id)
    except:
        og_phylogeny_ctx = {}

    context = {
        "valid_id": valid_id,
        "menu": True,
        "n_homologues": n_homologues,
        "og_id": format_orthogroup(og_id, to_url=True),
        "translation": translation,
        "seq": sequence,
        "sequence_type": feature_type,
        "n_swissprot_hits": n_swissprot_hits,
        "locus": locus,
        "feature_type" : "CDS",
        **cog_ctx,
        **kegg_ctx,
        **homolog_tab_ctx,
        **general_tab,
        **og_conserv_ctx,
        **og_phylogeny_ctx,
        **pfam_ctx,
        **genomic_region_ctx,
        **diamond_matches_ctx
    }
    return render(request, 'chlamdb/locus.html', my_locals(context))


def str_if_none(s):
    if s is None:
        return "-"
    return s


def search_bar_helper(fam_name, fam_url, results):
    val = getattr(results, fam_name)
    if val is None:
        val = "-"
    else:
        val = f"<a href=\"/{fam_url}/{val}\">{val}</a>"
    return val


def search_bar(request):
    db = db_utils.DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    option2status = db.get_config_table()
    index = sb.ChlamdbIndex.use_index(settings.SEARCH_INDEX)
    user_query = request.GET.get("accession")

    results = list(index.search(user_query, limit=None))

    if len(results) == 0:
        ctx = {"search_failed": True, "search_term": user_query}
        return render(request, "chlamdb/search.html", my_locals(ctx))

    search_results = []
    has_ko = option2status.get("KEGG", False)
    has_cog = option2status.get("COG", False)
    has_pfam = option2status.get("pfam", False)

    genes, cog, ko, pfam = [], [], [], []
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

    gene_active, cogs_active, ko_active, pfam_active = "active", "", "", ""
    if len(genes)==0:
        gene_active = ""
        if len(cog)>0:
            cogs_active = "active"
        elif len(ko)>0:
            ko_active = "active"
        elif len(pfam)>0:
            pfam_active = "active"

    genes_headers = ["Accession", "Gene", "Product", "Organism"]
    cog_headers = ["COG", "Description"]
    ko_headers = ["KO", "Description"]
    pfam_headers = ["PFAM domain", "Description"]
    insert_index = 4
    ctx = {"search_term": user_query,
            "gene_active": gene_active,
            "cogs_active": cogs_active,
            "ko_active": ko_active, 
            "pfam_active": pfam_active,
            "genes_headers": genes_headers,
            "genes": genes,
            "cog_headers": cog_headers,
            "pfam_headers": pfam_headers, 
            "ko_headers": ko_headers,
            "cogs": cog,
            "pfam": pfam,
            "ko": ko}
    return render(request, "chlamdb/search.html", my_locals(ctx))


def hydropathy(request, locus):
    biodb = settings.BIODB
    print ('hydropathy -- %s --%s' % (biodb, locus))

    from chlamdb.plots import hydrophobicity_plots

    fig = hydrophobicity_plots.locus2hydrophobicity_plot(biodb, locus)

    path = settings.BASE_DIR + '/assets/temp/hydro.png'
    fig.savefig(path,dpi=500)
    asset_path = '/temp/hydro.png'

    return render(request, 'chlamdb/hydropathy.html', my_locals(locals()))


def get_all_prot_infos(db, seqids, orthogroups):
    hsh_gene_locs = db.get_gene_loc(seqids, as_hash=True)
    hsh_prot_infos = db.get_proteins_info(seqids)
    hsh_organisms = db.get_organism(seqids, as_hash=True)
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
        if gene==None:
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


def format_locus(locus, to_url=False):
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
    if not R is None:
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
    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb_path)
    cog_id = int(cog_id[3:])

    if request.method != "GET":
        return render(request, 'chlamdb/fam.html', my_locals(locals()))

    df_seqid_to_cog = db.get_cog_hits([cog_id], indexing="seqid", search_on="cog", keep_taxid=True)
    if len(df_seqid_to_cog)==0:
        return render(request, 'chlamdb/fam.html', {"msg": f"No entry for {format_cog(cog_id)}"})

    seqids      = df_seqid_to_cog.index.tolist()
    
    orthogroups = db.get_og_count(seqids, search_on="seqid", keep_taxid=True)
    cog_info    = db.get_cog_summaries([cog_id], only_cog_desc=True, as_df=True)
    all_locus_data, group_count = get_all_prot_infos(db, seqids, orthogroups)
    ref_tree  = db.get_reference_phylogeny()
    cog_func = db.get_cog_code_description()

    df_cog_count = df_seqid_to_cog.groupby(["taxid"]).count()
    e_tree = tab_gen_profile_tree(db, df_cog_count.cog, format_cog(cog_id), orthogroups)
    asset_path = f"/temp/fam_tree_{cog_id}.svg"
    path = settings.BASE_DIR+"/assets/"+asset_path
    e_tree.render(path, dpi=500)

    func, cog_description = cog_info.loc[cog_id]
    info_func = "<br>".join((cog_func[code] for code in func))
    type = "cog"

    info = [info_func, cog_description]
    menu = True
    envoi = True
    return render(request, 'chlamdb/fam.html', my_locals(locals()))


def format_pathway(pat_id):
    return f"map{pat_id:05d}"


def format_module(mod_id, to_url=False):
    formated = f"M{mod_id:05d}"
    if to_url:
        return f"<a href=/KEGG_module_map/{formated}>{formated}</a>"
    return formated


def fam_ko(request, ko_str):
    ko_id = int(ko_str[len("K"):])
    db = db_utils.DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)

    df_ko_hits  = db.get_ko_hits([ko_id], search_on="ko", indexing="seqid", keep_taxid=True)
    seqids      = df_ko_hits.index.tolist()
    seqid_to_og = db.get_og_count(seqids, search_on="seqid", keep_taxid=True)
    red_color   = set(tuple(entry) for entry in seqid_to_og.to_numpy())

    pathways = db.get_ko_pathways([ko_id])
    modules  = db.get_ko_modules([ko_id])
    modules_id = [mod_id for key, values in modules.items() for mod_id, desc in values]
    modules_data = db.get_modules_info(modules_id)
    ko_desc      = db.get_ko_desc([ko_id])[ko_id]
    all_locus_data, group_count = get_all_prot_infos(db, seqids, seqid_to_og)

    pathway_data = format_ko_path(pathways, ko_id, as_list=True)
    module_data = [(format_ko_module(mod_id), cat, mod_desc)
            for mod_id, mod_desc, mod_def, path, cat in modules_data]
    df_ko_count  = df_ko_hits.groupby(["taxid"]).count()

    fam  = format_ko(ko_id)
    e_tree = tab_gen_profile_tree(db, df_ko_count.ko, fam, seqid_to_og)
    asset_path = f"/temp/fam_tree_{fam}.svg"
    path = settings.BASE_DIR + f"/assets/{asset_path}"

    e_tree.render(path, dpi=500)
    type = "ko"
    menu = True
    envoi = True
    return render(request, 'chlamdb/fam.html', my_locals(locals()))


def fam_pfam(request, pfam):
    context = {
        "type": "pfam",
        "menu" : True,
        "envoi" : True,
        "fam": pfam
    }
    if len(pfam)<2 or not pfam.startswith("PF"):
        # add error message
        return render(request, 'chlamdb/fam.html', my_locals(context))
    try:
        pfam_id = int(pfam[len("PF"):])
    except:
        return render(request, 'chlamdb/fam.html', my_locals(context))

    db = db_utils.DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    pfam_hits = db.get_pfam_hits([pfam_id], search_on="pfam", indexing="seqid", keep_taxid=True)
    seqids = pfam_hits.seqid.unique().tolist()
    orthogroups = db.get_og_count(seqids, search_on="seqid", keep_taxid=True)
    all_locus_data, group_count = get_all_prot_infos(db, seqids, orthogroups)
    infos = db.get_pfam_def([pfam_id])

    e_tree = tab_gen_profile_tree(db, pfam_hits.groupby(["taxid"])["pfam"].count(),
        pfam, orthogroups)

    asset_path = f"/temp/fam_tree_{pfam_id}.svg"
    path = settings.BASE_DIR+"/assets/"+asset_path
    e_tree.render(path, dpi=500)

    context["all_locus_data"] = all_locus_data
    context["group_count"] = group_count
    context["info"] = [infos.loc[pfam_id]["def"]] 
    context["asset_path"] = asset_path
    return render(request, 'chlamdb/fam.html', my_locals(context))


def COG_phylo_heatmap(request, frequency):
    biodb = settings.BIODB_DB_PATH

    if request.method != "GET":
        return render(request, 'chlamdb/COG_phylo_heatmap.html', my_locals(locals()))
    freq = frequency!="False"

    db = db_utils.DB.load_db_from_name(biodb)
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
        if len(func)==1:
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
        col = SimpleColorColumn(func_count, header=detailed_func + "("+func+")")
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
    if request.method != "GET":
        return render(request, 'chlamdb/KEGG_module_map.html', my_locals(locals()))
    db = db_utils.DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)

    try:
        module_id = int(module_name[len("M"):])
    except:
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
        if not ko in mat.columns:
            e_tree.add_column(SimpleColorColumn({}, header=format_ko(ko), default_val="-"))
            continue
        if not ko in ko_to_og_mapping.index:
            e_tree.add_column(SimpleColorColumn.fromSeries(mat[ko], header=format_ko(ko)))
            continue
            
        associated_ogs = ko_to_og_mapping.loc[ko]
        n_homologs = og_taxid[associated_ogs].sum(axis=1)
        hsh_col, hsh_val = {}, {}
        for taxid, count in mat[ko].iteritems():
            hsh_curr = hsh_pres[taxid]
            if count>0:
                hsh_curr[ko] = 1
                hsh_val[taxid] = count
                hsh_col[taxid] = EteTree.RED
            elif taxid in n_homologs.index:
                cnt = n_homologs.loc[taxid]
                hsh_val[taxid] = cnt
                if cnt>0:
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

    e_tree.add_column(SimpleColorColumn({}, default_val = " "))
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
        if not ko in mat.columns:
            continue
        if not ko in ko_to_og_mapping.index:
            e_tree.add_column(SimpleColorColumn.fromSeries(mat[ko], header=format_ko(ko)))
            continue
            
        associated_ogs = ko_to_og_mapping.loc[ko]
        n_homologs = og_taxid[associated_ogs].sum(axis=1)
        hsh_col, hsh_val = {}, {}
        for taxid, count in mat[ko].iteritems():
            if count>0:
                hsh_val[taxid] = count
                hsh_col[taxid] = EteTree.RED
            elif taxid in n_homologs.index:
                cnt = n_homologs.loc[taxid]
                hsh_val[taxid] = cnt
                if cnt>0:
                    hsh_col[taxid] = EteTree.GREEN
            else:
                hsh_val[taxid] = 0
        col_chooser = KOModuleChooser(hsh_col)
        e_tree.add_column(SimpleColorColumn(hsh_val, header=format_ko(ko),
            col_func=col_chooser.get_color))
    e_tree.rename_leaves(leaf_to_name)
    return e_tree


def KEGG_mapp_ko(request, map_name, taxon_id=None):
    db = db_utils.DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    pathway = int(map_name[len("map"):])
    kos = db.get_ko_pathways([pathway], search_on="pathway", as_df=True)

    if not taxon_id is None:
        taxid = int(taxon_id)

    ko_list = kos["ko"].unique().tolist()
    ko_hits = db.get_ko_hits(ko_list, search_on="ko", indexing="taxid")
    ko_desc = db.get_ko_desc(ko_list)
    ko_ttl_count = ko_hits.sum(axis=1)
    hsh_organisms = db.get_genomes_description().description.to_dict()

    header = ["KO", "Description", "All occurrences"]
    data = []
    if not taxon_id is None:
        header.insert(2, "#in this genome")
    for ko_id, descr in ko_desc.items():
        if ko_id in ko_ttl_count.index:
            ttl = ko_ttl_count.loc[ko_id]
        else:
            ttl = 0
        if ko_id in ko_hits.index and not taxon_id is None:
            in_this_genome = ko_hits[taxid].loc[ko_id]
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

    all_kos = []
    for col in e_tree.columns:
        # hack: we should not have access to the inner details
        # of a class, but I am in an hurry
        if not taxon_id is None and taxid in col.values:
            all_kos.append(col.header)
    ctx = {
        "pathway": kos.iloc[0].description,
        "header": header,
        "data": data,
        "asset_path": asset_path,
        "url": map_name+"+"+"+".join(all_kos)
    }

    if not taxon_id is None:
        ctx["organism"] = hsh_organisms[taxid]
    return render(request, 'chlamdb/KEGG_map_ko.html', my_locals(ctx))


def get_cog(request, taxon_id, category):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb)

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
            print("Cog id unknown: ", cog_id)
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
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

    targets = [int(i) for i in request.GET.getlist('h')]
    if len(targets)>5:
        targets = targets[0:6]

    cog_hits    = db.get_cog_hits(targets, indexing="taxid", search_on="taxid")
    genome_desc = db.get_genomes_description().description.to_dict()
    cog_description = db.get_cog_summaries(cog_hits.index.tolist(), only_cog_desc=True, as_df=True)
    selected_cogs   = cog_description[cog_description.function.str.contains(category)]
    cog_codes       = db.get_cog_code_description()

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
        data = ",".join(f"\"{format_cog(cog)}\"" for cog, count in non_zero_cogs.iteritems())
        series_tab.append( f"{{name: \"{genome_desc[target]}\", data: [{data}]}}" )
    series = "[" + ",".join(series_tab) + "]"

    cog_func_dict = (f"\"{func}\": \"{descr}\"" for func, descr in cog_codes.items())
    cog_func_dict = "{"+",".join(cog_func_dict)+"}"
    display_form = False
    envoi_venn = True
    return render(request, 'chlamdb/venn_cogs.html', my_locals(locals()))


def ko_venn_subset(request, category):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

    category = category.replace("+", " ")
    try:
        targets = [int(i) for i in request.GET.getlist('h')]
    except:
        # add an error message
        return render(request, 'chlamdb/venn_ko.html', my_locals(locals()))

    if len(targets)>5:
        targets = targets[0:6]

    genomes   = db.get_genomes_description(targets).description.to_dict()
    ko_counts = db.get_ko_count_cat(taxon_ids=targets, subcategory_name=category, index=False)
    ko_count  = ko_counts.groupby("taxon_id")["KO"].apply(list)

    # shameful copy/cape from venn_ko
    fmt_data = []
    ko_set = set()
    for taxid in targets:
        if not taxid in ko_count.index:
            continue
        kos = ko_count.loc[taxid]
        kos_str = ",".join(f"{to_s(format_ko(ko))}" for ko in kos)
        genome = genomes[taxid]
        fmt_data.append(f"{{name: {to_s(genome)}, data: [{kos_str}] }}")
        ko_set = ko_set.union({ko for ko in kos })
    series = "[" + ",".join(fmt_data) + "]"

    ko_list = list(ko_set)
    ko_descriptions = db.get_ko_desc(ko_list)
    ko2description = []
    for ko, ko_desc in ko_descriptions.items():
        forbidden = "\""
        ko_item = f"h[{to_s(format_ko(ko))}] = {forbidden}{ko_desc}{forbidden};"
        ko2description.append(ko_item)

    display_form = False
    envoi_venn = True
    return render(request, 'chlamdb/venn_ko.html', my_locals(locals()))


def module_cat_info(request, taxid, category):
    # Not really efficient code: it would be better
    # to get the cat_id associated with category to select
    # in the category list to avoid multiple string comparison

    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb_path, settings.BIODB_CONF)

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
    mid = ",".join(f"{to_s(bioentry)}:{to_s(description)}" for bioentry, description in hsh.items())
    return taxon_map + mid + "};"


def module_barchart(request):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

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
        string = f"{{ label: {to_s(taxid)}, values : [" + ",".join(str_entry_data) + "]}"
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

    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb)
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
    map_lst = (f"\"{func_descr}\" : \"{func}\"" for func, func_descr in category_dico.items())
    category_map = category_map + ",".join(map_lst) + '};'

    taxon_map = 'var taxon2description = {'
    taxon_map_lst = (f"\"{target}\":\"{taxon2description[target]}\"" for target in target_bioentries)
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
    labels = labels_template % ('"'+'","'.join(all_categories) + '"')
    envoi = True
    return render(request, 'chlamdb/cog_barplot.html', my_locals(locals()))


def blastnr_cat_info(request, accession, rank, taxon):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    target_accessions = [i for i in request.GET.getlist('h')]
    counttype = request.GET.getlist('t')[0]
    top_n = request.GET.getlist('n')[0]

    biodb_id_sql = 'select biodatabase_id from biodatabase where name="%s"' % biodb
    biodb_id = server.adaptor.execute_and_fetchall(biodb_id_sql,)[0][0]

    if counttype == 'Majority':
        sql = 'select B.locus_tag, A.%s ,A.n from (select seqfeature_id,%s, count(*) as n from blastnr_blastnr A ' \
              ' inner join blastnr_blastnr_taxonomy B on A.subject_taxid=B.taxon_id where hit_number<=%s and query_bioentry_id=%s ' \
              ' group by seqfeature_id, %s order by seqfeature_id,n DESC) A ' \
              ' inner join custom_tables_locus2seqfeature_id B on A.seqfeature_id=B.seqfeature_id' % (rank,
                                                                                                      rank,
                                                                                                      top_n,
                                                                                                      accession,
                                                                                                      rank)

        data = server.adaptor.execute_and_fetchall(sql,)
        category2count = {}
        all_query_locus_list = []
        majority_locus_list = []

        for i in data:
            # keep only the majoritary taxon
            if i[0] not in all_query_locus_list:
                # keep only data for taxon of interest
                if i[1] == taxon:
                    majority_locus_list.append(i[0])
                all_query_locus_list.append(i[0])
        locus_list = majority_locus_list

    elif counttype == 'BBH':
        sql = ' select locus_tag from (select t2.*,t1.locus_tag from custom_tables_locus2seqfeature_id t1 ' \
              ' inner join blastnr_blastnr t2 on t1.seqfeature_id=t2.seqfeature_id' \
              ' where hit_number=1) A inner join blastnr_blastnr_taxonomy B on A.subject_taxid=B.taxon_id ' \
              ' where %s="%s"  and query_bioentry_id=%s;' % (rank, taxon, accession)

        locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]
    else:
        raise 'invalide type'

    # get number of hits for each kingdom
    sql_superkingdom = 'select locus_tag,superkingdom, count(*) as n from' \
    ' (select t2.*,t1.locus_tag from custom_tables_locus2seqfeature_id t1 ' \
    ' inner join blastnr_blastnr t2 on t1.seqfeature_id=t2.seqfeature_id' \
    ' where locus_tag in ("%s")) A' \
    ' inner join blastnr_blastnr_taxonomy B on A.subject_taxid=B.taxon_id  group by locus_tag,superkingdom;' % (  '","'.join(locus_list))

    superkingdom_data = server.adaptor.execute_and_fetchall(sql_superkingdom,)
    locus2superkingdom_counts = {}
    for row in superkingdom_data:
        if row[0] not in locus2superkingdom_counts:
            locus2superkingdom_counts[row[0]] = {}
            locus2superkingdom_counts[row[0]][row[1]] = row[2]
        else:
            locus2superkingdom_counts[row[0]][row[1]] = row[2]
    # calculate sum of hits
    for locus_tag in locus2superkingdom_counts:
        total = sum([int(i) for i in locus2superkingdom_counts[locus_tag].values()])
        lst = []
        for superkingdom in locus2superkingdom_counts[locus_tag]:
            lst.append(["%s_percent" % superkingdom, round(float(locus2superkingdom_counts[locus_tag][superkingdom])/total, 2)])
        for kd in lst:
            locus2superkingdom_counts[locus_tag][kd[0]] = kd[1]
        locus2superkingdom_counts[locus_tag]["TOTAL"] = total


    locus2annot, \
    locus_tag2cog_catego, \
    locus_tag2cog_name, \
    locus_tag2ko, \
    pathway2category, \
    module2category, \
    ko2ko_pathways, \
    ko2ko_modules, \
    locus2interpro = get_locus_annotations( locus_list)

    accession2taxon = manipulate_biosqldb.accession2taxon_id(server, biodb)
    sql = 'select bioentry_id, taxon_id from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id' \
          ' where t1.name="%s"' % biodb
    accession2taxon = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    circos_url = '?ref=%s&' % locus2annot[0][-1]
    target_taxons = [str(accession2taxon[i]) for i in target_accessions]
    reference_taxon = str(accession2taxon[accession])
    target_taxons.pop(target_taxons.index(reference_taxon))
    circos_url += "t="+('&t=').join((target_taxons)) + '&h=' + ('&h=').join(locus_list)

    locus_filter = '"' + '","'.join(locus_list) + '"'
    sql = 'select locus_tag,subject_accession,subject_kingdom,subject_scientific_name,subject_taxid,evalue,percent_identity ' \
          ' from custom_tables. locus2seqfeature_id_%s t1 ' \
          ' inner join blastnr_blastnr t2 on t1.seqfeature_id=t2.seqfeature_id ' \
          ' where t1.locus_tag in (%s) and hit_number=1;' % (  locus_filter)

    locus_tag2blastnr_BBH = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

    '''
    series, \
    labels, \
    serie_all_counts, \
    serie_target_counts, \
    series_counts, \
    labels_counts, \
    category_description, \
    category_map, \
    n_missing_cog, \
    missing_cog_list = locus_tag2cog_series( locus_list, reference_taxon=None)
    '''


    return render(request, 'chlamdb/blastnr_info.html', my_locals(locals()))


# TODO: implement with plotly to avoid having to save a picture at 
# every request + weird errors regarding Qt
def pan_genome(request, type):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb)

    venn_form_class = make_venn_from(db, plasmid=False)

    if request.method != "POST":
        form = venn_form_class()
        return render(request, 'chlamdb/pan_genome.html', my_locals(locals()))

    form = venn_form_class(request.POST)
    if not form.is_valid():
        # should add an error message
        form = venn_form_class()
        return render(request, 'chlamdb/pan_genome.html', my_locals(locals()))

    import seaborn as sn
    import matplotlib.pyplot as plt
    import pandas as pd
    import time

    taxids = form.get_taxids()

    if type == "COG":
        df_hits = db.get_cog_hits(taxids, search_on="taxid", indexing="taxid")
        type_txt = "COG"
    elif type == "orthology":
        df_hits = db.get_og_count(taxids, search_on="taxid")
        type_txt = "orthologs"
    elif type == "ko":
        df_hits = db.get_ko_hits(taxids, search_on="taxid")
        type_txt = "KO"
    elif type=="Pfam":
        df_hits = db.get_pfam_hits(taxids, search_on="taxid")
        type_txt = "PFAM"
    else:
        # should add an error message
        form = venn_form_class()
        return render(request, 'chlamdb/pan_genome.html', my_locals(locals()))

    target2description = db.get_genomes_description().description.to_dict()
    df_hits.columns = [target2description[i] for i in df_hits.columns.values]
    path2 = settings.BASE_DIR + '/assets/temp/pangenome_barplot.svg'
    asset_path2 = '/temp/pangenome_barplot.svg'
    n_entries_per_genome = df_hits[df_hits>0].count(axis=0).sort_values()

    fig, ax = plt.subplots()
    barplot = ax.bar(list(n_entries_per_genome.index), list(n_entries_per_genome.values))
    ax.set_ylabel(f"{type} count")
    ax.set_xticklabels(list(n_entries_per_genome.index), rotation=45, horizontalalignment="right")
    fig.tight_layout()
    fig.savefig(path2)

    hsh_shared_count = {}
    hsh_total_count = {}
    for entry, pres in df_hits[n_entries_per_genome.index[0]].items():
        if pres==0:
            continue
        hsh_shared_count[entry] = 1
        hsh_total_count[entry] = 1

    total_entries = [n_entries_per_genome.values[0]]
    shared_entries = []
    for i in range(1, len(n_entries_per_genome)):
        cur_total = total_entries[-1]
        cur_shared = 0
        for entry, pres in df_hits[n_entries_per_genome.index[i]].items():
            if pres==0:
                continue
            v = hsh_shared_count.get(entry, 0)
            if v==i:
                cur_shared += 1
                hsh_shared_count[entry] = i+1
            if entry not in hsh_total_count:
                cur_total += 1
                hsh_total_count[entry] = 1
        total_entries.append(cur_total)
        shared_entries.append(cur_shared)

    path = settings.BASE_DIR + '/assets/temp/pangenome.svg'
    asset_path = '/temp/pangenome.svg'
    fig, ax = plt.subplots()

    ax.plot(total_entries)
    ax2 = ax.twinx()
    ax2.plot([i for i in range(1, len(n_entries_per_genome))], shared_entries, color="red")

    ax.set_xticks([i for i in range(0, len(n_entries_per_genome))])
    ax.set_xticklabels(n_entries_per_genome.index, rotation=45, horizontalalignment="right")
    ax2.set_ylabel(f"Number of shared {type_txt}")
    ax.set_ylabel(f"Number of {type_txt} in pangenome")

    fig.tight_layout()
    fig.savefig(path)

    # Serie with the number of genomes having a given cog
    missing_entries = df_hits[df_hits == 0].count(axis=1)
    core = len(missing_entries[missing_entries == 0])

    # cogs present in all but one genome
    core_minus1 = len(missing_entries[missing_entries == 1])
    total = len(df_hits.index)

    genome_count_list = []
    for i, count in enumerate(n_entries_per_genome):
        genome_count_list.append([i+1, count])
    envoi = True
    return render(request, 'chlamdb/pan_genome.html', my_locals(locals()))


def blastnr_overview(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 ' \
               ' on t1.biodatabase_id=t2.biodatabase_id where name="%s";'



    tree = server.adaptor.execute_and_fetchall(sql_tree)[0][0]


    header_list2 = ['effectiveT3',  'BPBAac', 'T3MM', 'T4SEpre_bpbAac', 'T4SEpre_psAac', 'chaperones', 'intesect']

    from chlamdb.phylo_tree_display import phylo_tree_bar
    from chlamdb.plots import hmm_heatmap
    from chlamdb.plots import pathway_heatmap
    from chlamdb.plots import module_heatmap

    set2taxon2value = {} #, column_names = hmm_heatmap.get_set_data(biodb, score_cutoff=20)

    sql_checkm_completeness = 'select taxon_id, completeness from custom_tables_checkm;'
    sql_checkm_n_duplicated = 'select taxon_id,n_total from custom_tables_checkm;'
    sql_n_contigs = 'select taxon_id,n_contigs from genomes_info t1 ' \
                              ' inner join bioentry t2 on t1.accession=t2.accession ' \
                              ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
                              ' where t3.name="%s" and t2.description not like "%%%%plasmid%%%%" ' \
                              ' group by taxon_id;' % (biodb)
    sql_n_no_CDS = 'select taxon_id,n_no_CDS from genomes_info t1 ' \
                              ' inner join bioentry t2 on t1.accession=t2.accession ' \
                              ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
                              ' where t3.name="%s" and t2.description not like "%%%%plasmid%%%%" ' \
                              ' group by taxon_id;' % (biodb)
    sql_n_no_BBH_chlamydiae = 'select taxon_id,n_no_BBH_chlamydiae from genomes_info t1 ' \
                              ' inner join bioentry t2 on t1.accession=t2.accession ' \
                              ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
                              ' where t3.name="%s" and t2.description not like "%%%%plasmid%%%%" ' \
                              ' group by taxon_id;' % (biodb)

    sql_n_no_BBH_chlamydiae = 'select taxon_id,n_no_BBH_chlamydiae from genomes_info t1 ' \
                              ' inner join bioentry t2 on t1.accession=t2.accession ' \
                              ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
                              ' where t3.name="%s" and t2.description not like "%%%%plasmid%%%%" ' \
                              ' group by taxon_id;' % (biodb)

    sql_GC = 'select taxon_id,GC from genomes_info t1 ' \
                              ' inner join bioentry t2 on t1.accession=t2.accession ' \
                              ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
                              ' where t3.name="%s" and t2.description not like "%%%%plasmid%%%%" ' \
                              ' group by taxon_id;' % (biodb)

    sql_genome_size = 'select taxon_id,genome_size/1000000 from genomes_info t1 ' \
                              ' inner join bioentry t2 on t1.accession=t2.accession ' \
                              ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id ' \
                              ' where t3.name="%s" and t2.description not like "%%%%plasmid%%%%" ' \
                              ' group by taxon_id;' % (biodb)

    taxon_id2completeness = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_checkm_completeness))
    taxon_id2duplicate = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_checkm_n_duplicated))
    taxon_id2n_no_CDS = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_n_no_CDS))
    taxon_id2n_contigs = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_n_contigs))
    taxon_id2n_contigs_no_BBH_chlamydiae = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_n_no_BBH_chlamydiae))
    taxon_id2GC = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_GC))
    taxon_id2genome_size = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_genome_size))

    taxon_id2GC['127'] = 50
    # , 'T6SSi' 'rinke_et_al_2013', 'dupont_et_al_2012', 'eisen_et_al_2013'
    my_sets = ['duplicates', 'contigs','contigs_no_CDS', 'contigs_no_BBH', 'ADP/ATP carrier protein']

    #
    set2taxon2value['duplicates'] = taxon_id2duplicate
    set2taxon2value['contigs'] = taxon_id2n_contigs
    set2taxon2value['contigs_no_CDS'] = taxon_id2n_no_CDS
    set2taxon2value['contigs_no_BBH'] = taxon_id2n_contigs_no_BBH_chlamydiae


    set2taxon2value_list_simple_barplot = {}
    set2taxon2value_list_simple_barplot['gc'] = taxon_id2GC
    set2taxon2value_list_simple_barplot['complet.'] = taxon_id2completeness
    set2taxon2value_list_simple_barplot['size'] = taxon_id2genome_size

    header3 = ['gc', 'size', 'complet.' ]

    sql_ntt_transporters = 'select * from comparative_tables_interpro where id="IPR004667";'
    ntt_data = list(server.adaptor.execute_and_fetchall(sql_ntt_transporters)[0])
    sql_headers = 'show columns from  comparative_tables_interpro'
    ntt_taxid = [i[0] for i in server.adaptor.execute_and_fetchall(sql_headers)]

    taxon_id2count_ntt = dict(zip(ntt_taxid[1:], ntt_data[1:]))

    set2taxon2value['ADP/ATP carrier protein'] = taxon_id2count_ntt

    flagellum_data = pathway_heatmap.pathway_list2profile_dico(biodb, ["map02040"])
    peptidoglycan_data = pathway_heatmap.pathway_list2profile_dico(biodb, ["map00550"])
    purines_data = pathway_heatmap.pathway_list2profile_dico(biodb, ["map00230"])
    pyrim_data = pathway_heatmap.pathway_list2profile_dico(biodb, ["map00240"])

    set2taxon2value['Flagellum'] = flagellum_data[1]['Flagellar assembly']
    set2taxon2value['Peptidoglycan biosynthesis'] = peptidoglycan_data[1]['Peptidoglycan biosynthesis']
    set2taxon2value['Purine metabolism'] = purines_data[1]['Purine metabolism']
    set2taxon2value['Pyrimidine metabolism'] = pyrim_data[1]['Pyrimidine metabolism']

    module_list, code2taxon2count_modules = module_heatmap.module_list2profile_dico(biodb, ["M00001","M00004","M00009"])

    set2taxon2value['gylcolysis'] = code2taxon2count_modules['Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate [PATH:map01200 map00010]']
    set2taxon2value['TCA'] = code2taxon2count_modules['Citrate cycle (TCA cycle, Krebs cycle) [PATH:map01200 map00020]']
    set2taxon2value['PPP'] = code2taxon2count_modules['Pentose phosphate pathway (Pentose phosphate cycle) [PATH:map01200 map00030]']

    my_sets+= ["T3SS",
               "T4SS",
               'Flagellum',
                  "chemosensory",
                  'Peptidoglycan biosynthesis',
                  "Menaquinone",
                  "Menaquinone/futalosine pathway",
                  "NADH-quinone oxidoreductase",
                  "NQR",
                  "NDH-2",
                  "SDH",
                  "Cytochrome bc1 complex",
                  "Cbb3-type cytochrome c oxidase",
                  "Cytochrome bo(3) ubiquinol oxidase",
                  "Cytochrome bd-I ubiquinol oxidase",
                  "F-type ATPase",
                  "V-type ATPase",
                   ]

    my_sets+=['TCA',
              'gylcolysis',
              'PPP', "glycogen metabolism",'Purine metabolism','Pyrimidine metabolism']

    ######################
    # taxonomy project: locus int
    ######################
    #my_sets.append('ntt')
    try:
        # group by orthogroup ==> even if one prot has multiple homologs, count only one
        sql_complex = 'select C.taxon_id, count(*) as n from (select distinct A.category, A.orthogroup, B.taxon_id ' \
              ' from (select distinct t1.category,t2.orthogroup from custom_tables_annot_table t1 ' \
              ' inner join orthology_detail t2 on t1.locus_tag=t2.locus_tag ' \
              ' where category="%s") A inner join orthology_detail B ' \
              ' on A.orthogroup=B.orthogroup) C group by C.category,C.taxon_id;'
        # not group by orthogroup
        sql_simple = 'select C.taxon_id, count(*) as n from (select A.category, A.orthogroup, B.taxon_id ' \
              ' from (select distinct t1.category,t2.orthogroup from custom_tables_annot_table t1 ' \
              ' inner join orthology_detail t2 on t1.locus_tag=t2.locus_tag ' \
              ' where category="%s") A inner join orthology_detail B ' \
              ' on A.orthogroup=B.orthogroup) C group by C.category,C.taxon_id;'

        c_list = ["chemosensory",
                  "T4SS",
                  "T3SS",
                  "F-type ATPase",
                  "V-type ATPase",
                  "glycogen metabolism",
                  "Menaquinone",
                  "Menaquinone/futalosine pathway",
                  "NADH-quinone oxidoreductase",
                  "SDH",
                  "NQR",
                  "Cytochrome bc1 complex",
                  "Cbb3-type cytochrome c oxidase",
                  "Cytochrome bo(3) ubiquinol oxidase",
                  "Cytochrome bd-I ubiquinol oxidase"
                   ]
        for category in c_list:
            sql = sql_complex % ("2017_06_29b_motile_chlamydiae", biodb, category, biodb)
            #sql_ftype = sql_complex % (biodb, biodb, biodb)
            #sql_T3SS = sql_complex % (biodb, biodb, biodb)
            #sql_T4SS = sql_complex % (biodb, biodb, biodb)
            #sql_chemo = sql_complex % (biodb, biodb, biodb)
            #sql_futa = sql_complex % (biodb, biodb, biodb)
            #sql_mena = sql_complex % (biodb, biodb, biodb)

            #taxon_id2euo = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_euo,))
            #taxon_id2vtype = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_vtype,))
            #taxon_id2ftype = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_ftype,))
            #taxon_id2T3SS = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_T3SS,))
            #taxon_id2T4SS = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_T4SS,))
            #taxon_id2chemo = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_chemo,))
            #taxon_id2futa = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_futa,))
            taxon_id2counts = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

            set2taxon2value[category] = taxon_id2counts
            #my_sets.append(category)
            #set2taxon2value['euo'] = taxon_id2euo
            #set2taxon2value['V-type ATPase'] = taxon_id2vtype
            #set2taxon2value['F-type ATPase'] = taxon_id2ftype
            #set2taxon2value['T3SS'] = taxon_id2T3SS
            #set2taxon2value['T4SS'] = taxon_id2T4SS
            #set2taxon2value['chemosensory'] = taxon_id2chemo
            #set2taxon2value['Menaquinone/futalosine pathway'] = taxon_id2futa
            #set2taxon2value['Menaquinone'] = taxon_id2mena

        #my_sets+=['chemosensory', 'V-type ATPase', 'F-type ATPase', 'T4SS', 'T3SS', 'Menaquinone/futalosine pathway',
        #          'Menaquinone']

        s_list = ["NDH-2"]#

        s_list = ["NDH-2",
                  "uhpC",
                  "CPAF",
                  "CopN",
                  "NUE",
                  "pkn5",
                  "mip",
                  "Lda3",
                  "Tarp",
                  "TepP",
                  "CT_847",
                  "CT_868",
                  "CT_867",
                  "CT_082",
                  "capN",
                  "CT_695",
                  "CT_620",
                  "CT_621",
                  "CT_622",
                  "CT_711",
                  "CT_694",
                  "CT_163",
                  "CT_365",
                  "CT_610",
                  "CT_157",
                  "CT_119",
                  "CT_115",
                  "CT_118",
                  "CT_228",
                  "CT_229",
                  "CT_813",
                  "CT_850",
                  "CPn0585",
                  "CPn0517",
                  "euo",
                  "OmcB",
                  "OmpA",
                  "OmpA1",
                  "pmp1",
                  "pmp2",
                  "Ctad1",
                  "FtsH",
                  "FtsI/Pbp3",
                  "FtsY",
                  "FtsK",
                  "FtsQ",
                  "FtsL",
                  "FtsW",
                  "MreB",
                  "RodZ",
                  "Ddl",
                  "DacB",
                  "DacC",
                  "DapA",
                  "DapB",
                  "DapL(AspC3)",
                  "DapF",
                  "MraY",
                  "MurJ (MviN)",
                  "Pbp2",
                  "AmiA",
                  "AmiB",
                  "NlpD",
                  "Pal",
                  "TolQ",
                  "TolR",
                  "TolA",
                  "TolA",
                  "TolB"]

        for locus in s_list:
            if locus == 'ntt':
                continue
            sql_locus = sql_simple % ("2017_06_29b_motile_chlamydiae",locus)
            taxid2locus_count = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_locus,))
            set2taxon2value[locus] = taxid2locus_count
            my_sets.append(locus)
    except:
        pass




    #######################
    # metabolism aa, cofactors, nucleotides

    sql = 'select taxon_id,count(*) from ' \
          ' (select distinct t1.taxon_id,t1.ko_id,pathway_category from enzyme_locus2ko t1 ' \
          ' inner join enzyme_pathway2ko_v1 t2 on t1.ko_id=t2.ko_id ' \
          ' inner join enzyme_kegg_pathway t3 on t2.pathway_id=t3.pathway_id) A ' \
          ' where pathway_category="%s"  group by A.taxon_id,pathway_category;'

    #taxon_id2nucleo = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql % (biodb,"1.4 Nucleotide metabolism",)))
    #set2taxon2value['Nucleotide metabolism'] = taxon_id2nucleo
    #my_sets.append('Nucleotide metabolism')
    taxon_id2nucleo = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql % ("1.8 Metabolism of cofactors and vitamins",)))
    set2taxon2value['Cofactors and vitamins'] = taxon_id2nucleo
    my_sets.append('Cofactors and vitamins')
    taxon_id2nucleo = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql % ("1.5 Amino acid metabolism",)))
    set2taxon2value['AA metabolism'] = taxon_id2nucleo
    my_sets.append('AA metabolism')
    #taxon_id2nucleo = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql % ("1.10 Biosynthesis of other secondary metabolites",)))
    #set2taxon2value['1.10 Biosynthesis of other secondary metabolites'] = taxon_id2nucleo
    #my_sets.append('1.10 Biosynthesis of other secondary metabolites')
    #taxon_id2nucleo = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql % ("1.1 Carbohydrate metabolism",)))
    #set2taxon2value['Carbohydrate metabolism'] = taxon_id2nucleo
    #my_sets.append('Carbohydrate metabolism')
    taxon_id2nucleo = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql % ("1.3 Lipid metabolism",)))
    set2taxon2value['Lipid metabolism'] = taxon_id2nucleo
    my_sets.append('Lipid metabolism')


    sql = 'select taxon_id, count(*) from COG_seqfeature_id2best_COG_hit t1 ' \
          ' inner join COG_cog_names_2014 t2 on t1.hit_cog_id=t2.COG_id ' \
          ' inner join COG_cog_id2cog_category ta on t1.hit_cog_id=ta.COG_id' \
          ' inner join COG_code2category tb on ta.category_id=tb.category_id' \
          ' inner join annotation_seqfeature_id2locus tc on t1.seqfeature_id=tc.seqfeature_id' \
          ' where tb.code="X" group by taxon_id;'

    taxon_id2mobile_elements = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
    set2taxon2value['modile_elements'] = taxon_id2mobile_elements
    my_sets.append('modile_elements')




    set2taxon2value_new = {}
    for set in my_sets:
        set2taxon2value_new[set] = {}
        for taxon in set2taxon2value[set]:

            try:
                value = set2taxon2value[set][taxon]
                #if value < 4:
                #    value = 0
                set2taxon2value_new[set][taxon] = value
            except:
                set2taxon2value_new[set][taxon] = 0

    sql = 'select taxon_id, n_no_hits, n_less_100_hits, n_100_hits from blastnr_count_n_blast order by n_no_hits;'

    # blast_hits_taxonomy_overview
    sql2 = 'select * from blastnr_BBH_taxo_hit_number_1;'
    sql3 = 'select * from blastnr_BBH_taxo_hit_number_2;'


    sql4 = 'select t2.taxon_id, sum(n_CDS) from genomes_info t1 ' \
           ' inner join bioentry t2 on t1.ACCESSION=t2.ACCESSION ' \
           ' inner join biodatabase t3 on t2.biodatabase_id=t3.biodatabase_id where t3.name="%s" group by taxon_id;' % (biodb)


    taxon_id2blast_count = server.adaptor.execute_and_fetchall(sql,)
    taxon_id2BBH_1_taxonomy = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))
    taxon_id2BBH_2_taxonomy = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql3,))
    taxon_id2proteome_size = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql4,))




    taxon_id2values = {}
    for row in taxon_id2blast_count:
        taxo1 = list(taxon_id2BBH_1_taxonomy[row[0]])
        taxo2 = list(taxon_id2BBH_2_taxonomy[row[0]])
        taxo1.append(int(taxon_id2proteome_size[str(row[0])])-sum(taxo1))
        taxo2.append(int(taxon_id2proteome_size[str(row[0])])-sum(taxo2))
        taxo_filter = [taxo1[0], sum([taxo1[1], taxo1[2], taxo1[3], taxo1[4],taxo1[5]]), taxo1[6]]
        taxo_filter2 = [taxo2[0], sum([taxo2[1], taxo2[2], taxo2[3], taxo2[4],taxo2[5]]), taxo2[6]]
        taxon_id2values[row[0]] = [list(reversed(row[1:])), taxo1, taxo2]
    header_list = ['nr_n_hits', "nr_tax.1","nr_tax.2"]


    sql = 'select taxon_id, count(*) as n, best_hit_phylum from blastnr_blastnr_majority_phylum' \
          ' group by taxon_id,best_hit_phylum order by taxon_id,n DESC;'

    sql2 = 'select taxon_id, count(*) as n, majority_phylum from blastnr_blastnr_majority_phylum ' \
           ' group by taxon_id,majority_phylum order by taxon_id,n DESC;'

    data = server.adaptor.execute_and_fetchall(sql,)
    data2 = server.adaptor.execute_and_fetchall(sql2,)


    taxon2most_freq_phylum = {}
    for row in data:
        if row[0] not in taxon2most_freq_phylum:
            taxon2most_freq_phylum[row[0]] = row[2]
    taxon_match = []
    for row in data2:
        if row[0] not in taxon_match:
            taxon2most_freq_phylum[row[0]] += "/%s" % row[2]
            taxon_match.append(row[0])

    tree1, style1 = phylo_tree_bar.plot_tree_stacked_barplot(tree,
                                                    taxon_id2values,
                                                    header_list,
                                                    taxon2set2value_heatmap=set2taxon2value_new,
                                                    header_list2=my_sets,
                                                    biodb=biodb,
                                                    general_max=False,
                                                    taxon2label=taxon2most_freq_phylum,
                                                    set2taxon2value_list_simple_barplot=set2taxon2value_list_simple_barplot,
                                                    header_list3=header3) # taxon2most_freq_phylum

    # col = '#fc8d59' # col = '#91bfdb' '#99d594'
    path1 = settings.BASE_DIR + '/assets/temp/interpro_tree2.svg'
    asset_path1 = '/temp/interpro_tree2.svg'
    tree1.render(path1, dpi=600, tree_style=style1)
    return render(request, 'chlamdb/blastnr_overview.html', my_locals(locals()))



def blastnr_top_non_phylum(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    blastnr_form_class = make_blastnr_best_non_top_phylum_form(biodb)

    if request.method == 'POST':
        form = blastnr_form_class(request.POST)

        if form.is_valid():
            import pandas as pd

            accessions = form.cleaned_data['accessions']
            selection = form.cleaned_data['selection']

            taxon_filter = ','.join(accessions)

            if selection == 'all':
                sql = 'select t4.locus_tag,t5.product,t5.gene,t1.hit_number,t1.percent_identity,t3.kingdom,t3.class,' \
                      ' t3.order,t3.family,t3.species,t1.subject_accession,t1.subject_title from blastnr_blastnr_best_non_self_phylum t1' \
                      ' inner join blastnr_blastnr_taxonomy t3 on t1.subject_taxon_id=t3.taxon_id ' \
                      ' inner join custom_tables_locus2seqfeature_id t4 on t1.seqfeature_id=t4.seqfeature_id ' \
                      ' inner join orthology_detail t5 on t4.locus_tag=t5.locus_tag  ' \
                      'where t1.superkingdom="Eukaryota"' \
                      ' and query_taxon_id in (%s)' % (taxon_filter)

            if selection == 'specific':
                sql = 'select t4.locus_tag,t5.product,t5.gene,t1.hit_number,t1.percent_identity,t3.kingdom,t3.class,' \
                      't3.order,t3.family,t3.species,t1.subject_accession,t1.subject_title from blastnr_blastnr_best_non_self_phylum t1 ' \
                      ' inner join custom_tables_seqfeature_id2n_species t2 on t1.seqfeature_id=t2.seqfeature_id ' \
                      ' inner join blastnr_blastnr_taxonomy t3 on t1.subject_taxon_id=t3.taxon_id ' \
                      ' inner join custom_tables_locus2seqfeature_id t4 on t1.seqfeature_id=t4.seqfeature_id ' \
                      ' inner join orthology_detail t5 on t4.locus_tag=t5.locus_tag  ' \
                      ' where t1.superkingdom="Eukaryota" and n_species=1  and t1.query_taxon_id in (%s);' % (taxon_filter)

            data = server.adaptor.execute_and_fetchall(sql,)
            envoi = True

    else:  
        form = blastnr_form_class()
    return render(request, 'chlamdb/blastnr_locus_list.html', my_locals(locals()))

def blastnr_barchart(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    blastnr_form_class = make_blastnr_form(biodb)

    if request.method == 'POST':
        form = blastnr_form_class(request.POST)

        if form.is_valid():
            import pandas as pd

            target_accessions = form.cleaned_data['accession']
            rank = form.cleaned_data['rank']
            counttype = form.cleaned_data['type']
            top_n = form.cleaned_data['top_number']

            biodb_id_sql = 'select biodatabase_id from biodatabase where name="%s"'
            biodb_id = server.adaptor.execute_and_fetchall(biodb_id_sql,)[0][0]
            sql_accession = 'select bioentry_id,description from bioentry where biodatabase_id=%s and bioentry_id in (%s)' % (biodb_id,
                                                                                                                              '"'+'","'.join(target_accessions)+'"')
            taxon2description = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql_accession,))

            data_all_accessions = []
            for accession in target_accessions:

                if counttype == 'Majority':
                    sql = 'select seqfeature_id,%s, count(*) as n from blastnr_blastnr A ' \
                          ' inner join blastnr_blastnr_taxonomy B on A.subject_taxid=B.taxon_id where query_bioentry_id=%s and hit_number<=%s' \
                          ' group by seqfeature_id,%s order by seqfeature_id, n DESC' % (rank,
                                                                                         accession,
                                                                                         top_n,
                                                                                         rank)

                    data = server.adaptor.execute_and_fetchall(sql,)
                    category2count = {}
                    query_locus_list = []
                    for n, i in enumerate(data):
                        # KEEP ONY the first match (highest count ordered with mysql)
                        if i[0] not in query_locus_list:
                            if i[1] not in category2count:
                                category2count[i[1]] = 1
                            else:
                                category2count[i[1]] += 1
                            query_locus_list.append(i[0])
                            try:
                                if data[n+1][0] == data[n][0]:
                                    if data[n+1][2] == data[n][2]:
                                        print("Idem!:",data[n+1], data[n])
                            except:
                                pass

                    data = zip(category2count.keys(), category2count.values())

                elif counttype == 'BBH':
                    sql = ' select %s, count(*) as n from (select * from blastnr_blastnr ' \
                          ' where query_bioentry_id=%s and hit_number=1) A inner join blastnr_blastnr_taxonomy B on A.subject_taxid=B.taxon_id ' \
                          ' group by %s;' % (rank,
                                             accession,
                                             rank)
                    data = server.adaptor.execute_and_fetchall(sql,)
                else:
                    raise 'invalide type'

                for i in data:
                    data_all_accessions.append((accession,) + i)

            taxon_map = 'var taxon2description = {'
            for i in taxon2description:
                taxon_map += '"%s":"%s",' % (i, taxon2description[i])
            taxon_map = taxon_map[0:-1] + '};'

            taxon2category2count = {}
            all_categories = []
            for line in data_all_accessions:
                if line[0] not in taxon2category2count:
                    taxon2category2count[line[0]] = {}
                    taxon2category2count[line[0]][line[1]] = line[2]
                else:
                    taxon2category2count[line[0]][line[1]] = line[2]
                if line[1] not in all_categories:
                    all_categories.append(line[1])
            labels_template = '[\n' \
                              '%s\n' \
                              ']\n'

            serie_template = '[%s\n' \
                             ']\n'

            one_serie_template = '{\n' \
                                 'label: "%s",\n' \
                                 'values: [%s]\n' \
                                 '},\n'

            # count number of hits for each category and order based on the number of hits
            category2count = {}
            for taxon in taxon2category2count:
                for category in all_categories:
                    try:
                        if category not in category2count:
                            category2count[category] = int(taxon2category2count[taxon][category])
                        else:
                            category2count[category] += int(taxon2category2count[taxon][category])
                    except:
                        pass
            for key in category2count:
                print("%s\t%s" % (key, category2count[key]))
            data = pd.DataFrame({'category': list(category2count.keys()),
                    'count': list(category2count.values()) })
            data_sort = data.sort_values("count", ascending=0)

            all_series_templates = []
            for taxon in taxon2category2count:
                one_category_list = []
                for category in data_sort['category']:
                    try:
                        one_category_list.append(taxon2category2count[taxon][category])
                    except:
                        one_category_list.append(0)
                one_category_list = [str(i) for i in one_category_list]
                all_series_templates.append(one_serie_template % (taxon, ','.join(one_category_list)))

            series = serie_template % ''.join(all_series_templates)
            cat_list = [str(i) for i in data_sort['category']]
            labels = labels_template % ('"'+'","'.join(cat_list) + '"')


            circos_url = '?h=' + ('&h=').join(target_accessions) + '&t=%s&n=%s' % (counttype, top_n)


            envoi = True
    else:  
        form = blastnr_form_class()
    return render(request, 'chlamdb/blastnr_best_barplot.html', my_locals(locals()))


def format_gene(gene):
    if pd.isna(gene):
        return "-"
    else:
        return gene


def blastswissprot(request, locus_tag):
    print("FOO")
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

    if locus_tag is None:
        return render(request, 'chlamdb/blastswiss.html', my_locals({"valid_id": False}))
    try:
        seqid, _ = db.get_seqid(locus_tag=locus_tag)
    except:
        return render(request, 'chlamdb/blastswiss.html', my_locals({"valid_id": False}))

    swissprot_homologs = db.get_swissprot_homologs([seqid])
    transl = db.get_translation(seqid)

    # NOTE: to fix, some queries have a > 100% coverage. Need to use query start and query end.
    swissprot_homologs["perc_cover"] = (swissprot_homologs.match_len*100/len(transl)).round(2)
    blast_data = []
    for num, data in swissprot_homologs.iterrows():
        line = [num+1, data.accession, data.evalue, data.bitscore, data.perc_id,
                data.gaps, data.perc_cover, data.pe, data.gene, data.definition,
                data.organism, data.taxid]
        blast_data.append(line)
    ctx = {
        "valid_id": not swissprot_homologs.empty,
        "locus_tag": locus_tag,
        "blast_data": blast_data
    }
    return render(request, 'chlamdb/blastswiss.html', my_locals(ctx))


def blastnr(request, locus_tag):
    biodb = settings.BIODB

    print ('blastnr hits %s -- %s' % (biodb, locus_tag))

    server, db = manipulate_biosqldb.load_db(biodb)

    if request.method == 'GET': 

        server, db = manipulate_biosqldb.load_db(biodb)

        sql = 'select accession, organism from orthology_detail where locus_tag="%s"' % (biodb, locus_tag)
        data = server.adaptor.execute_and_fetchall(sql,)[0]
        accession = data[0]
        organism = data[1]

        server, db = manipulate_biosqldb.load_db(biodb)
        columns = 'hit_number, subject_accession, subject_kingdom, subject_scientific_name, ' \
                  ' subject_taxid, subject_title, evalue, bit_score, percent_identity, gaps, length'
        sql = 'select hit_number,subject_accession,superkingdom,subject_scientific_name,subject_taxid,subject_title,evalue,bit_score,percent_identity,gaps, length,' \
              ' B.phylum, B.order, B.family from (select %s from custom_tables_locus2seqfeature_id t1 inner join blastnr_blastnr t2' \
              ' on t1.seqfeature_id=t2.seqfeature_id where locus_tag="%s" order by hit_number) A ' \
              ' inner join blastnr_blastnr_taxonomy B on A.subject_taxid=B.taxon_id;' % (columns, locus_tag)

        blast_data = list(server.adaptor.execute_and_fetchall(sql))

        blast_data = list(server.adaptor.execute_and_fetchall(sql))

        if len(blast_data) > 0:
            valid_id = True
            blast_query_locus = blast_data[0][1]
            blast_query_protein_id = blast_data[0][2]


        return render(request, 'chlamdb/blastnr.html', my_locals(locals()))


    return render(request, 'chlamdb/blastnr.html', my_locals(locals()))


def plot_region(request):
    raise RuntimeException("To be re-implemented in d3.js")


def orthogroups(request):
    if request.method == 'POST':
        form = BiodatabaseForm(request.POST)
        if form.is_valid():
            biodb = form.cleaned_data['biodatabase']
            groups = "orthogroup_size_distrib_%s.svg" % biodb
            envoi = True
    else:
        form = BiodatabaseForm()

    return render(request, 'chlamdb/orthogroups.html', my_locals(locals()))


def sitemap(request):
    from io import StringIO
    map = '''<?xml version="1.0" encoding="UTF-8"?>
<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">
<url><loc>http://chlamdb.ch/</loc></url><url><loc>http://chlamdb.ch/about</loc></url>
</urlset>'''

    strio = StringIO()
    strio.write(map)
    response = HttpResponse(content_type='text/xml')
    #response.content_type = "text/xml"
    #response['Content-Disposition'] = 'attachment; filename="sitemap.xml"'
    response.write(strio.getvalue())
    return response


def circos_main(request):
    biodb = settings.BIODB

    if request.method == 'POST':
        
        reference_taxon = request.POST["reference_taxid"]
        include_taxids = eval(request.POST["include_taxids"])
        exclude_taxids = eval(request.POST["exclude_taxids"])
        og_list = eval(request.POST["og_list"])
        
        target_taxons = include_taxids + exclude_taxids
        
        
        print("reference_taxon", reference_taxon)
        print("target_taxons", target_taxons)
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
                sql_order = 'select taxon_2 from comparative_tables_core_orthogroups_identity_msa where taxon_1=%s order by identity desc;' % (reference_taxon)
                ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order)]
                target_taxons = ordered_taxons[0:10]
            except:
                sql_order = 'select taxon_2 from comparative_tables_shared_orthogroups where taxon_1=%s order by n_shared_orthogroups DESC;' % (reference_taxon)

                ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order)]
                target_taxons = ordered_taxons[0:10]
        else:
            target_taxons = [int(i) for i in request.GET.getlist('t')]
        highlight = request.GET.getlist('h')

        #sql = 'select locus_tag,traduction from orthology_detail_k_cosson_05_16 where orthogroup in (%s) and accession="NC_016845"' % ('"'+'","'.join(highlight)+'"')

        task = run_circos_main.delay(reference_taxon, target_taxons, highlight)
        task_id = task.id
        envoi_circos = True
        envoi_region = True

        #return HttpResponse(json.dumps({'task_id': task.id}), content_type='application/json')

    return render(request, 'chlamdb/circos_main.html', my_locals(locals()))


def circos_blastnr(request):
    biodb = settings.BIODB
    from chlamdb.plots import gbk2circos
    circos_form_class = make_genome_selection_form(biodb)
    server, db = manipulate_biosqldb.load_db(biodb)

    if request.method == 'POST':

        form = circos_form_class(request.POST)

        if form.is_valid():

            reference_taxon = form.cleaned_data['genome']

            description2accession_dict = manipulate_biosqldb.description2accession_dict(server, biodb)

            reference_accessions = manipulate_biosqldb.taxon_id2accessions(server, reference_taxon, biodb)

            record_list = []
            for accession in reference_accessions:

                biorecord = cache.get(biodb + "_" + accession)

                if not biorecord:
                    new_record = db.lookup(accession=accession)
                    biorecord = SeqRecord(Seq(new_record.seq.data, new_record.seq.alphabet),
                                                             id=new_record.id, name=new_record.name,
                                                             description=new_record.description,
                                                             dbxrefs =new_record.dbxrefs,
                                                             features=new_record.features,
                                                             annotations=new_record.annotations)
                    record_id = biorecord.id.split(".")[0]
                    cache.set(biodb + "_" + record_id, biorecord)
                    record_list.append(biorecord)
                else:
                    record_list.append(biorecord)


            ref_name = ''
            for i in reference_accessions:
                ref_name += i
            circos_file = "circos/%s.svg" % ref_name
            from chlamdb.biosqldb import circos
            from chlamdb.biosqldb import shell_command
            import ete3



            outtput_dir = settings.BASE_DIR + "/assets/circos/"

            myplot = circos.CircosAccession2blastnr_plot(server,
                                                         biodb,
                                                         record_list,
                                                         outtput_dir,
                                                         taxon_list=[],
                                                         highlight_BBH=True)

            original_map_file = settings.BASE_DIR + "/assets/circos/%s.html" % ref_name

            with open(original_map_file, "r") as f:
                map_string = ''.join([line for line in f.readlines()])

            circos_html = '<!DOCTYPE html>\n' \
                          ' <html>\n' \
                          ' <body>\n' \
                          ' %s\n' \
                          ' <img src="%s.svg" usemap="#%s">' \
                          ' </body>\n' \
                          ' </html>\n' % (map_string, ref_name, ref_name)

            circos_new_file = '/assets/circos/circos_clic_%s.html' % ref_name

            with open(settings.BASE_DIR + circos_new_file, "w") as f:
                f.write(circos_html)

            #target_map_file = settings.BASE_DIR + "/templates/circos/%s.html" % ref_name
            circos_svg_file = "circos/%s.svg" % ref_name
            circos_png_file = "circos/%s.png" % ref_name
            original_map_file_svg = settings.BASE_DIR + "/assets/circos/%s.svg" % ref_name
            #target_map_file_svg = settings.BASE_DIR + "/templates/circos/%s.svg" % ref_name
            map_file = "circos/%s.html" % ref_name
            svg_file = "circos/%s.svg" % ref_name
            #a, b, c = shell_command.shell_command("mv %s %s" % (original_map_file, target_map_file))
            #a, b, c = shell_command.shell_command("cp %s %s" % (original_map_file_svg, target_map_file_svg))
            map_name = ref_name
            envoi_circos = True
            envoi_region = True
    else:
        form = circos_form_class()
    return render(request, 'chlamdb/circos_blastnr.html', my_locals(locals()))



def get_circos_data(reference_taxon, target_taxons, highlight_og=False):

    from metagenlab_libs import circosjs
    
    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb_path)
    
    #task = run_circos.delay(reference_taxon, target_taxons)
    
    if highlight_og:
        df_genes = db.get_genes_from_og(highlight_og, taxon_ids=[reference_taxon], terms=["locus_tag"])
        locus_list = df_genes["locus_tag"].to_list()
    else:
        locus_list = []
    
    # "bioentry_id", "accession" ,"length"
    df_bioentry = db.get_bioentry_list(reference_taxon, min_bioentry_length=1000)
    
    # "bioentry_id", "seqfeature_id", "start_pos", "end_pos", "strand"
    df_feature_location = db.get_features_location(reference_taxon, ["CDS", "rRNA", "tRNA"]).set_index(["seqfeature_id"])
    
    # retrieve n_orthologs of list of seqids

    seq_og   = db.get_og_count(df_feature_location.index.to_list(), search_on="seqid")
    count_all_genomes = db.get_og_count(seq_og["orthogroup"].to_list(), search_on="orthogroup")
    orthogroup2count_all = count_all_genomes[count_all_genomes > 0].count(axis=1)
    homologs_count = df_feature_location.loc[df_feature_location.term_name=='CDS'].join(seq_og).reset_index().set_index("orthogroup").merge(orthogroup2count_all.rename('value'), left_index=True, right_index=True)[["bioentry_id", "start_pos", "end_pos", "value"]]

    df_identity = db.get_identity_closest_homolog(reference_taxon, target_taxons).set_index(["target_taxid"])
    
    c = circosjs.CircosJs()
    
    c.add_contigs_data(df_bioentry)
    
    # sort taxons by number of homologs (from mot similar to most dissmilar)
    target_taxon_n_homologs = df_identity.groupby(["target_taxid"]).count()["seqfeature_id_1"].sort_values(ascending=False)
    locus2seqfeature_id = db.get_hsh_locus_to_seqfeature_id()
    seqfeature_id2locus_tag = {value:key for key, value in locus2seqfeature_id.items()}
    heatmap_data_list = {}
    # "bioentry_id", "seqfeature_id", "start_pos", "end_pos", "strand"
    # "seqfeature_id_1", "seqfeature_id_2", "identity", "target_taxid"
    # join on seqfeature id
    df_feature_location["gene"] = df_feature_location["qualifier_value_gene"].fillna("-")
    df_feature_location["gene_product"] = df_feature_location["qualifier_value_product"].fillna("-")
    df_feature_location["locus_tag"] = df_feature_location["qualifier_value_locus_tag"].fillna("-")
    
    df_feature_location["color"]  = "grey"
    df_feature_location["color"][df_feature_location["term_name"] == 'tRNA']  = "magenta"
    df_feature_location["color"][df_feature_location["term_name"] == 'rRNA']  = "magenta"
    
    df_feature_location = df_feature_location.rename(columns={"locus_tag":"locus_ref"})
    #= [seqfeature_id2locus_tag[seqfeature_id] for seqfeature_id in df_feature_location.index]
    minus_strand = df_feature_location.set_index("strand").loc[-1]
    plus_strand = df_feature_location.set_index("strand").loc[1]
    c.add_histogram_track(homologs_count, "n_genomes", radius_diff=0.1, sep= 0.005, outer=True)
    c.add_gene_track(minus_strand, "minus", sep=0, radius_diff=0.04, highlight_list=locus_list)
    c.add_gene_track(plus_strand, "plus", sep=0, radius_diff=0.04, highlight_list=locus_list)

    # iterate ordered list of target taxids, add track to circos    
    for n, target_taxon in enumerate(target_taxon_n_homologs.index):
        df_combined = df_feature_location.join(df_identity.loc[target_taxon].reset_index().set_index("seqfeature_id_1")).reset_index()
        df_combined.identity = df_combined.identity.fillna(0).astype(int)
        df_combined.bioentry_id = df_combined.bioentry_id.astype(str)
        
        # only keep the highest identity for each seqfeature id
        df_combined = df_combined.sort_values('identity', ascending=False).drop_duplicates('index').sort_index()
        df_combined["locus_tag"] = [seqfeature_id2locus_tag[seqfeature_id] if not pd.isna(seqfeature_id) else None for seqfeature_id in df_combined["seqfeature_id_2"]]
        # comp is a custom scale with missing orthologs coloured in light blue
        if n == 0:
            sep = 0.03
        else:
            sep= 0.01
        c.add_heatmap_track(df_combined, f"target_{target_taxon}", color="comp", sep=sep, radius_diff=0.04)

    c.add_line_track(df_bioentry, "GC_content", radius_diff= 0.12, fillcolor="green")
    
    js_code = c.get_js_code()
    return js_code


def circos(request):
    
    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb_path)

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
    #local_vars["form"] = form
    
    return render(request, 'chlamdb/circos.html', local_vars)


def alignment(request, input_fasta):

    handle = open(input_fasta, "rU")
    for record in SeqIO.parse(handle, "fasta"):
        pass
    return render(request, 'chlamdb/alignment.html', my_locals(locals()))


def blast_profile(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    if request.method == 'POST': 

        form = BlastProfileForm(request.POST, request.FILES)

        if form.is_valid():  
            from tempfile import NamedTemporaryFile
            from io import StringIO
            from chlamdb.phylo_tree_display import biosqldb_plot_blast_hits_phylo
            from chlamdb.biosqldb import biosql_own_sql_tables
            from ete3 import Tree, TreeStyle

            fasta_file = request.FILES['fasta_file']
            fasta_string = StringIO(request.FILES['fasta_file'].read().decode("UTF-8"))
            fasta_rec = [i for i in SeqIO.parse(fasta_string, 'fasta')]

            try:
                ordered_labels_all = [i.id.split('|')[1] for i in fasta_rec]
            except:
                ordered_labels_all = [i.id for i in fasta_rec]

            blast_type = form.cleaned_data['blast']

            #my_record = SeqIO.read(request.FILES['fasta_file'].open(), 'fasta')


            tree, style1, tree2, style2, tree3, style3, locus2taxon2locus_closest = biosqldb_plot_blast_hits_phylo.plot_BBH_phylo(fasta_rec, 
                                                                                                                                  biodb,
                                                                                                                                  settings.BASE_DIR + '/assets/')


            sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

            tt = server.adaptor.execute_and_fetchall(sql_tree)[0][0]

            t1 = Tree(tt)

            R = t1.get_midpoint_outgroup()
            t1.set_outgroup(R)
            t1.ladderize()
            ordered_taxons_all = [str(i.name) for i in t1.iter_leaves()]

            sql = 'select taxon_id, t2.description from biodatabase t1 inner join bioentry t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                  ' where t1.name="%s" and t2.description not like "%%%%plasmid%%%%"' % (biodb)

            taxon2genome = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))
            taxon_list = taxon2genome.keys()


            ordered_taxons_keep=[]
            for taxon in ordered_taxons_all:
                if taxon in taxon_list:
                    ordered_taxons_keep.append(taxon)

            ordered_labels_keep=[]
            for label in ordered_labels_all:
                if label in locus2taxon2locus_closest:
                    ordered_labels_keep.append(label)

            path = settings.BASE_DIR + '/assets/temp/profile_tree.svg'
            path2 = settings.BASE_DIR + '/assets/temp/profile_tree2.svg'
            path3 = settings.BASE_DIR + '/assets/temp/profile_tree3.svg'
            asset_path = '/temp/profile_tree.svg'
            asset_path2 = '/temp/profile_tree2.svg'
            asset_path3 = '/temp/profile_tree3.svg'

            tree3.render(path3, dpi=500, tree_style=style3)
            tree.render(path2, dpi=500, tree_style=style1)
            tree2.render(path, dpi=500, tree_style=style2)


            envoi = True

    else:  
        form = BlastProfileForm()

    return render(request, 'chlamdb/blast_profile.html', my_locals(locals()))


def blast(request):
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)

    blast_form_class = make_blast_form(biodb)

    if request.method == 'POST': 

        form = blast_form_class(request.POST)

        if form.is_valid():  
            from Bio.Blast.Applications import NcbiblastpCommandline
            from Bio.Blast.Applications import NcbiblastnCommandline
            from Bio.Blast.Applications import NcbitblastnCommandline
            from Bio.Blast.Applications import NcbiblastxCommandline
            from Bio.Seq import reverse_complement, translate
            from tempfile import NamedTemporaryFile
            from io import StringIO
            from Bio.Blast import NCBIXML
            from Bio.Alphabet import IUPAC
            from Bio.Alphabet import _verify_alphabet
            import os
            from chlamdb.biosqldb import shell_command
            import re


            input_sequence = form.cleaned_data['blast_input']

            target_accession = form.cleaned_data['target']

            blast_type = form.cleaned_data['blast']


            unknown_format = False

            if '>' in input_sequence:
                my_record = [i for i in SeqIO.parse(StringIO(input_sequence), 'fasta', alphabet=IUPAC.ambiguous_dna)]
                seq = my_record[0]
                seq.alphabet = IUPAC.ambiguous_dna
                if _verify_alphabet(seq):
                    my_record = [i for i in SeqIO.parse(StringIO(input_sequence), 'fasta', alphabet=IUPAC.ambiguous_dna)]
                else:
                    seq.alphabet = IUPAC.extended_protein
                    if _verify_alphabet(seq):
                        my_record = [i for i in SeqIO.parse(StringIO(input_sequence), 'fasta', alphabet=IUPAC.extended_protein)]
                    else:
                        unknown_format = True
            else:
                input_sequence = input_sequence.rstrip(os.linesep)
                seq = Seq(input_sequence)
                seq.alphabet = IUPAC.ambiguous_dna
                if _verify_alphabet(seq):
                    my_record = [SeqRecord(seq, id="INPUT", description="INPUT")]
                else:
                    seq.alphabet = IUPAC.extended_protein
                    if _verify_alphabet(seq):
                        my_record = [SeqRecord(seq, id="INPUT", description="INPUT")]
                    else:
                        unknown_format = True
            if not unknown_format:
                if my_record[0].seq.alphabet == IUPAC.extended_protein and blast_type in ["blastn_ffn", "blastn_fna", "blastx"]:
                    wrong_format = True
                elif my_record[0].seq.alphabet == IUPAC.ambiguous_dna and blast_type in ["blastp", "tblastn"]:
                    wrong_format = True
                else:
                    query_file = NamedTemporaryFile(mode='w')
                    SeqIO.write(my_record, query_file, "fasta")
                    query_file.flush()

                    if blast_type=='blastn_ffn':
                        blastType = 'locus'
                        blastdb = settings.BASE_DIR + "/assets/%s/ffn/%s.ffn" % (biodb, target_accession)
                        blast_cline = NcbiblastnCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)
                    if blast_type=='blastn_fna':
                        blastType = 'genome'
                        blastdb = settings.BASE_DIR + "/assets/%s/fna/%s.fna" % (biodb, target_accession)
                        blast_cline = NcbiblastnCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)
                    if blast_type=='blastp':
                        blastType = 'locus'
                        blastdb = settings.BASE_DIR + "/assets/%s/faa/%s.faa" % (biodb, target_accession)
                        blast_cline = NcbiblastpCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)
                    if blast_type=='tblastn':
                        blastType = 'genome'
                        blastdb = settings.BASE_DIR + "/assets/%s/fna/%s.fna" % (biodb, target_accession)
                        blast_cline = NcbitblastnCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)
                        blast_cline2 = NcbitblastnCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=5)
                    if blast_type=='blastx':
                        blastType = 'locus'
                        blastdb = settings.BASE_DIR + "/assets/%s/faa/%s.faa" % (biodb, target_accession)
                        blast_cline = NcbiblastxCommandline(query=query_file.name, db=blastdb, evalue=10, outfmt=0)

                    blast_stdout, blast_stderr = blast_cline()

                    if blast_type=='tblastn':
                        from Bio.SeqUtils import six_frame_translations

                        blast_stdout2, blast_stderr2 = blast_cline2()

                        blast_records = NCBIXML.parse(StringIO(blast_stdout2))
                        all_data = []
                        best_hit_list = []
                        for record in blast_records:
                            for n, alignment in enumerate(record.alignments):
                                accession = alignment.title.split(' ')[1]

                                sql = 'select description from bioentry where accession="%s" ' % accession

                                description = server.adaptor.execute_and_fetchall(sql,)[0][0]

                                for n2, hsp in enumerate(alignment.hsps):
                                    if n == 0 and n2 == 0:
                                        best_hit_list.append([record.query, hsp.sbjct_start, hsp.sbjct_end])
                                
                                    start = hsp.sbjct_start
                                    end = hsp.sbjct_end
                                    if start > end:
                                        start = hsp.sbjct_end
                                        end = hsp.sbjct_start
                                    length = end-start
                                    leng = end-start
                                    seq = manipulate_biosqldb.location2sequence(server, accession, biodb, start, leng)
                                    
                                    anti = reverse_complement(seq)
                                    comp = anti[::-1]
                                    length = len(seq)
                                    frames = {}
                                    for i in range(0, 3):
                                        fragment_length = 3 * ((length-i) // 3)
                                        tem1 = translate(seq[i:i+fragment_length], 1)
                                        frames[i+1] = '<span style="color: #181407;">%s</span><span style="color: #bb60d5;">%s</span><span style="color: #181407;">%s</span>' % (tem1[0:100], tem1[100:len(tem1)-99], tem1[len(tem1)-99:])
                                        tmp2 = translate(anti[i:i+fragment_length], 1)[::-1]
                                        frames[-(i+1)] = tmp2
                                    all_data.append([accession, start, end, length, frames[1], frames[2], frames[3], frames[-1], frames[-2], frames[-3], description, seq])
                        if len(best_hit_list) > 0:
                            fig_list = []
                            for best_hit in best_hit_list:
                                accession = best_hit[0]
                                best_hit_start = best_hit[1]
                                best_hit_end = best_hit[2]
                                temp_location = os.path.join(settings.BASE_DIR, "assets/temp/")
                                temp_file = NamedTemporaryFile(delete=False, dir=temp_location, suffix=".svg")
                                name = 'temp/' + os.path.basename(temp_file.name)
                                fig_list.append([accession, name])
                                orthogroup_list = mysqldb_plot_genomic_feature.location2plot(db,
                                                                                            biodb,
                                                                                            target_accession,
                                                                                            temp_file.name,
                                                                                            best_hit_start-15000,
                                                                                            best_hit_end+15000,
                                                                                            cache,
                                                                                            color_locus_list = [],
                                                                                            region_highlight=[best_hit_start, best_hit_end])


                    no_match = re.compile('.* No hits found .*', re.DOTALL)

                    if no_match.match(blast_stdout):
                        print ("no blast hit")
                        blast_no_hits = blast_stdout
                    elif len(blast_stderr) != 0:
                        print ("blast error")
                        blast_err = blast_stderr
                    else:
                        from Bio import SearchIO
                        blast_record = [i for i  in SearchIO.parse(StringIO(blast_stdout), 'blast-text')]
                        all_locus_tag = []
                        for query in blast_record:
                            for hit in query:
                                locus_tag = hit.id
                                all_locus_tag.append(locus_tag)

                        locus_filter = '"' + '","'.join(all_locus_tag) + '"'
                        sql = 'select locus_tag, product from orthology_detail_%s where locus_tag in (%s)' % (biodb, locus_filter)
                        locus2product = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql,))

                        rand_id = id_generator(6)

                        blast_file_l = settings.BASE_DIR + '/assets/temp/%s.xml' % rand_id
                        f = open(blast_file_l, 'w')
                        f.write(blast_stdout)
                        f.close()

                        asset_blast_path = '/assets/temp/%s.xml' % rand_id
                        js_out = True

                envoi = True

    else:  
        form = blast_form_class()

    return render(request, 'chlamdb/blast.html', my_locals(locals()))


def get_record_from_memory(biodb, cache_obj, record_key, accession):

        biorecord = cache_obj.get(record_key)
        if not biorecord:
            new_record = biodb.lookup(accession=accession)
            biorecord = SeqRecord(Seq(new_record.seq.data, new_record.seq.alphabet),
                                                             id=new_record.id, name=new_record.name,
                                                             description=new_record.description,
                                                             dbxrefs =new_record.dbxrefs,
                                                             features=new_record.features,
                                                             annotations=new_record.annotations)
            cache_obj.set(record_key, biorecord)
        return biorecord



def circos2genomes(request):
    biodb = settings.BIODB
    from chlamdb.biosqldb import circos
    from chlamdb.biosqldb import shell_command
    server, db = manipulate_biosqldb.load_db(biodb)
    circos2genomes_form_class = make_circos2genomes_form(biodb)



    if request.method == 'POST': 
        #make_circos2genomes_form

        form = circos2genomes_form_class(request.POST)
        if form.is_valid():  
            valid_id = True
            server, db = manipulate_biosqldb.load_db(biodb)
            reference_taxon = form.cleaned_data['reference_genome']
            query_taxon = form.cleaned_data['query_genome']
            import re
            protein_locus_list = re.sub(" ", "", form.cleaned_data['locus_list'])
            protein_locus_list = list(protein_locus_list.split(","))

            #reference_taxon = manipulate_biosqldb.description2taxon_id(server, reference_genome, biodb)
            #query_taxon = manipulate_biosqldb.description2taxon_id(server, query_genome, biodb)

            reference_accessions = manipulate_biosqldb.taxon_id2accessions(server, reference_taxon, biodb)
            query_accessions = manipulate_biosqldb.taxon_id2accessions(server, query_taxon, biodb)

            reference_records = []
            for accession in reference_accessions:
                reference_records.append(get_record_from_memory(db, cache, biodb + "_" + accession, accession))

            query_records = []
            for accession in query_accessions:

                query_records.append(get_record_from_memory(db, cache, biodb + "_" + accession, accession))

            orthogroup_list = []
            if len(protein_locus_list[0]) > 0:
                for protein in protein_locus_list:

                    sql = 'select orthogroup from orthology_detail where protein_id="%s" or locus_tag="%s"' % (protein, 
                                                                                                               protein)

                    try:
                        protein_group = server.adaptor.execute_and_fetchall(sql,)[0][0]
                        orthogroup_list.append(protein_group)
                    except IndexError:
                        valid_id = False

            if valid_id:

                #accession2description = manipulate_biosqldb.accession2description_dict(server, biodb)
                taxon_id2description = manipulate_biosqldb.taxon_id2genome_description(server, biodb)
                reference_name = taxon_id2description[reference_taxon]
                query_name = taxon_id2description[query_taxon]

                accession2taxon_id = manipulate_biosqldb.accession2taxon_id(server, biodb)
                taxon_id_reference = accession2taxon_id[reference_accessions[0]]
                taxon_id_query = accession2taxon_id[query_accessions[0]]

                reference_n_orthogroups = manipulate_biosqldb.get_genome_number_of_orthogroups(server, biodb, taxon_id_reference)
                reference_n_proteins = manipulate_biosqldb.get_genome_number_of_proteins(server, biodb, taxon_id_reference)

                query_n_orthogroups = manipulate_biosqldb.get_genome_number_of_orthogroups(server, biodb, taxon_id_query)
                query_n_proteins = manipulate_biosqldb.get_genome_number_of_proteins(server, biodb, taxon_id_query)

                n_shared_orthogroups = manipulate_biosqldb.get_number_of_shared_orthogroups(server, biodb, taxon_id_reference, taxon_id_query)

                from chlamdb.biosqldb import circos

                path = settings.BASE_DIR + "/assets/circos"



                biplot = circos.CircosAccession2biplot(server, db, biodb, reference_records, query_records,
                                                       orthogroup_list, path)

                reference_file = "circos/%s" % biplot.reference_circos
                query_file = "circos/%s" % biplot.query_circos

            envoi = True

    else:  
        form = circos2genomes_form_class()

    return render(request, 'chlamdb/circos2genomes.html', my_locals(locals()))


def circos2genomes_main(request):
    return render(request, 'chlamdb/circos2genomes_main.html', my_locals(locals()))


def plot_heatmap(request, type):
    import seaborn as sns
    from datetime import datetime

    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb)
    form_class = make_venn_from(db)

    if request.method != "POST":
        form_venn = form_class()
        return render(request, 'chlamdb/plot_heatmap.html', my_locals(locals()))

    form_venn = form_class(request.POST)
    if not form_venn.is_valid():
        return render(request, 'chlamdb/plot_heatmap.html', my_locals(locals()))

    taxon_ids = form_venn.get_taxids()
    if type=="COG":
        mat = db.get_cog_hits(taxon_ids, indexing="taxid", search_on="taxid")
    elif type=="orthology":
        mat = db.get_og_count(taxon_ids)
    elif type == "ko":
        mat = db.get_ko_hits(taxon_ids)
    elif type == "Pfam":
        mat = db.get_pfam_hits(taxon_ids)
    else:
        form_venn = form_class()
        return render(request, 'chlamdb/plot_heatmap.html', my_locals(locals()))

    target2description = db.get_genomes_description().description.to_dict()
    mat.columns = (target2description[i] for i in mat.columns.values)
    cur_time = datetime.now().strftime("%H%M%S")

    filename = f"heatmap_{cur_time}.png"
    path = settings.BASE_DIR + '/assets/temp/' + filename
    asset_path = '/temp/' + filename
    cm = sns.clustermap(mat)
    cm.savefig(path)
    envoi_heatmap = True
    return render(request, 'chlamdb/plot_heatmap.html', my_locals(locals()))


def format_pathway(path_id, to_url=False, taxid=None):
    base_string = f"map{path_id:05d}"
    if not to_url:
        return base_string
    if taxid is None:
        to_page = f"/KEGG_map_ko/{base_string}"
    else:
        to_page = f"/KEGG_mapp_ko/{base_string}/{taxid}"
    return f"<a href=\"{to_page}\">{base_string}</a>"


def priam_kegg(request):
    db = db_utils.DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    priam_form_class = make_priam_form(db)
    if request.method != "POST":
        form = priam_form_class()
        return render(request, 'chlamdb/priam_kegg.html', my_locals(locals()))

    form = priam_form_class(request.POST)
    if not form.is_valid():
        form = priam_form_class()
        return render(request, 'chlamdb/priam_kegg.html', my_locals(locals()))

    taxid, _ = form.get_genome()
    ko_hits = db.get_ko_hits([taxid], search_on="taxid")
    kos = ko_hits.index.tolist()
    ko_pathways = db.get_ko_pathways(kos)
    header = ["Pathway", "Description", "Number of KO"]
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

    ctx = {"envoi":True, "data":data, "header": header}
    return render(request, 'chlamdb/priam_kegg.html', my_locals(ctx))


def locus_list2circos(request, target_taxon):
    biodb = settings.BIODB
    from chlamdb.plots import gbk2circos
    from chlamdb.biosqldb import circos
    import re

    server, db = manipulate_biosqldb.load_db(biodb)

    locus_list = [i for i in request.GET.getlist('l')]

    # 1. check if locus_list is part of the target genome
    # 2. if not part of the target genome, get orthogroup id

    locus_filter = '"'+'","'.join(locus_list)+'"'
    sql = 'select locus_tag, product from orthology_detail where taxon_id=%s and locus_tag in (%s)' % (target_taxon,
                                                                                                       locus_filter)
    data = server.adaptor.execute_and_fetchall(sql,)
    locus_target_genome = [i[0] for i in data]

    locus2label = {}
    for row in data:
        locus2label[row[0]] = row[0] #+ "-%s" % re.sub(" ","-",row[1])

    locus_other_genomes = []
    for locus in locus_list:
        if locus not in locus_target_genome:
            locus_other_genomes.append(locus)
    locus_other_genomes_filter = '"' + '","'.join(locus_other_genomes) + '"'
    sql2 = 'select B.locus_tag, A.orthogroup, B.product from (select orthogroup from orthology_detail where locus_tag in (%s) group by orthogroup) A' \
           ' inner join orthology_detail B on A.orthogroup=B.orthogroup where taxon_id=%s' % (locus_other_genomes_filter,
                                                                                              target_taxon)

    data = server.adaptor.execute_and_fetchall(sql2,)

    for row in data:
        locus2label[row[0]] = row[1] #+ "-%s" % re.sub(" ","-",row[2])

    reference_accessions = manipulate_biosqldb.taxon_id2accessions(server, target_taxon, biodb) # ["NC_009648"] NZ_CP009208 NC_016845

    record_list = []
    for accession in reference_accessions:

        biorecord = cache.get(biodb + "_" + accession)

        if not biorecord:
            new_record = db.lookup(accession=accession)
            biorecord = SeqRecord(Seq(new_record.seq.data, new_record.seq.alphabet),
                                                     id=new_record.id, name=new_record.name,
                                                     description=new_record.description,
                                                     dbxrefs =new_record.dbxrefs,
                                                     features=new_record.features,
                                                     annotations=new_record.annotations)
            record_id = biorecord.id.split(".")[0]
            cache.set(biodb + "_" + record_id, biorecord)
            record_list.append(biorecord)
        else:
            record_list.append(biorecord)

    ref_name = ('').join(reference_accessions)

    circos_file = "circos/%s.svg" % ref_name

    draft_data = []
    for biorecord in record_list:
        draft_data.append(gbk2circos.circos_fasta_draft_misc_features(biorecord))

    home_dir = os.path.dirname(__file__)

    temp_location = os.path.join(home_dir, "../assets/circos/")

    myplot = circos.CircosAccession2multiplot(server,
                                              db,
                                              biodb,
                                              record_list,
                                              [],
                                              locus_highlight=locus2label.keys(),
                                              out_directory=temp_location,
                                              draft_fasta=draft_data,
                                              href="/chlamdb/locusx/",
                                              ordered_taxons=[],
                                              locus2label=locus2label,
                                              show_homologs=False,
                                              radius=0.3)

    original_map_file = settings.BASE_DIR + "/assets/circos/%s.html" % ref_name
    with open(original_map_file, "r") as f:
        map_string = ''.join([line for line in f.readlines()])

    circos_html = '<!DOCTYPE html>\n' \
                  ' <html>\n' \
                  ' <body>\n' \
                  ' %s\n' \
                  ' <img src="%s.svg" usemap="#%s">' \
                  ' </body>\n' \
                  ' </html>\n' % (map_string, ref_name, ref_name)


    circos_new_file = '/assets/circos/circos_clic.html'

    with open(settings.BASE_DIR + circos_new_file, "w") as f:
        f.write(circos_html)

    original_map_file_svg = settings.BASE_DIR + "/assets/circos/%s.svg" % ref_name
    map_file = "circos/%s.html" % ref_name
    svg_file = "circos/%s.svg" % ref_name
    map_name = ref_name

    envoi_circos = True


    return render(request, 'chlamdb/locus2circos.html', my_locals(locals()))


def hmm2circos(request):
    from chlamdb.phylo_tree_display import ete_motifs
    biodb = settings.BIODB
    server, db = manipulate_biosqldb.load_db(biodb)
    hmm_form = hmm_sets_form_circos(biodb)

    if request.method == 'POST': 
        form = hmm_form(request.POST)
        if form.is_valid():
            from chlamdb.plots import hmm_heatmap
            from chlamdb.plots import gbk2circos
            from chlamdb.biosqldb import circos

            sql_biodb_id = 'select biodatabase_id from biodatabase where name="%s"' % biodb

            database_id = server.adaptor.execute_and_fetchall(sql_biodb_id,)[0][0]

            hmm_set = form.cleaned_data['hmm_set']
            reference_taxon = form.cleaned_data['genome']
            score_cutoff = form.cleaned_data['score_cutoff']
            query_coverage_cutoff = form.cleaned_data['query_coverage_cutoff']


            sql = 'select locus_tag from hmm.hmm_sets t1 ' \
                  ' inner join hmm.hmm_sets_entry t2 on t1.set_id=t2.set_id ' \
                  ' inner join hmm_hmm_hits_annotated_genome t3 on t2.hmm_id=t3.hmm_id' \
                  ' inner join custom_tables_locus2seqfeature_id t4 on t3.seqfeature_id=t4.seqfeature_id ' \
                  ' where t1.name="%s" and t3.taxon_id=%s and bitscore>=%s ' \
                  ' and query_coverage>=%s order by bitscore;' % (hmm_set,
                                                                  reference_taxon,
                                                                  score_cutoff,
                                                                  query_coverage_cutoff)

            target_locus_list = [i[0] for i in server.adaptor.execute_and_fetchall(sql,)]


            sql2 = 'select locus_tag,t5.name  from hmm.hmm_sets t1 ' \
                  ' inner join hmm.hmm_sets_entry t2 on t1.set_id=t2.set_id ' \
                  ' inner join hmm_hmm_hits_annotated_genome t3 on t2.hmm_id=t3.hmm_id' \
                  ' inner join custom_tables_locus2seqfeature_id t4 on t3.seqfeature_id=t4.seqfeature_id ' \
                  ' inner join hmm.hmm_profiles t5 on t2.hmm_id=t5.hmm_id' \
                  ' where t1.name="%s" and t3.taxon_id=%s and bitscore>=%s ' \
                   ' and query_coverage>=%s order by bitscore;' % (hmm_set, 
                                                                    reference_taxon, 
                                                                    score_cutoff, 
                                                                    query_coverage_cutoff)

            locus2label = manipulate_biosqldb.to_dict(server.adaptor.execute_and_fetchall(sql2,))

            description2accession_dict = manipulate_biosqldb.description2accession_dict(server, biodb)

            reference_accessions = manipulate_biosqldb.taxon_id2accessions(server, reference_taxon, biodb) # ["NC_009648"] NZ_CP009208 NC_016845


            record_list = []
            for accession in reference_accessions:
                biorecord = cache.get(biodb + "_" + accession)

                if not biorecord:
                    new_record = db.lookup(accession=accession)
                    biorecord = SeqRecord(Seq(new_record.seq.data, new_record.seq.alphabet),
                                                             id=new_record.id, name=new_record.name,
                                                             description=new_record.description,
                                                             dbxrefs =new_record.dbxrefs,
                                                             features=new_record.features,
                                                             annotations=new_record.annotations)
                    record_id = biorecord.id.split(".")[0]
                    cache.set(biodb + "_" + record_id, biorecord)
                    record_list.append(biorecord)
                else:
                    record_list.append(biorecord)

            ref_name = ('').join(reference_accessions)

            circos_file = "circos/%s.svg" % ref_name

            draft_data = []
            for biorecord in record_list:
                draft_data.append(gbk2circos.circos_fasta_draft_misc_features(biorecord))

            home_dir = os.path.dirname(__file__)

            temp_location = os.path.join(home_dir, "../assets/circos/")

            #sql_tree = 'select tree from reference_phylogeny as t1 inner join biodatabase as t2 on t1.biodatabase_id=t2.biodatabase_id where name="%s";' % biodb

            sql_order1 = 'select A.taxon_1 from (select taxon_1,median_identity from comparative_tables_shared_og_av_id where taxon_2=%s ' \
                        ' union select taxon_2,median_identity from comparative_tables_shared_og_av_id ' \
                        ' where taxon_1=%s order by median_identity DESC) A;' % (reference_taxon, 
                                                                                 reference_taxon)
            try:
                sql_order = 'select taxon_2 from comparative_tables_core_orthogroups_identity_msa where taxon_1=%s order by identity desc;' % (reference_taxon)

                ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order)]
            except:
                sql_order2 = 'select taxon_2 from comparative_tables_shared_og_av_id where taxon_1=%s order by median_identity desc;' % (reference_taxon)

                ordered_taxons = [i[0] for i in server.adaptor.execute_and_fetchall(sql_order1)]

            myplot = circos.CircosAccession2multiplot(server,
                                                      db,
                                                      biodb,
                                                      record_list,
                                                      [],
                                                      locus_highlight=target_locus_list,
                                                      out_directory=temp_location,
                                                      draft_fasta=draft_data,
                                                      href="/chlamdb/locusx/",
                                                      ordered_taxons = ordered_taxons,
                                                      locus2label=locus2label,
                                                      show_homologs=False)

            original_map_file = settings.BASE_DIR + "/assets/circos/%s.html" % ref_name
            with open(original_map_file, "r") as f:
                map_string = ''.join([line for line in f.readlines()])

            circos_html = '<!DOCTYPE html>\n' \
                          ' <html>\n' \
                          ' <body>\n' \
                          ' %s\n' \
                          ' <img src="%s.svg" usemap="#%s">' \
                          ' </body>\n' \
                          ' </html>\n' % (map_string, ref_name, ref_name)


            circos_new_file = '/assets/circos/circos_clic.html'

            with open(settings.BASE_DIR + circos_new_file, "w") as f:
                f.write(circos_html)

            original_map_file_svg = settings.BASE_DIR + "/assets/circos/%s.svg" % ref_name
            map_file = "circos/%s.html" % ref_name
            svg_file = "circos/%s.svg" % ref_name
            map_name = ref_name

            envoi_circos = True

    else:  
        form = hmm_form()

    return render(request, 'chlamdb/hmm2circos.html', my_locals(locals()))


def blast_sets(request):
    biodb = settings.BIODB
    from chlamdb.phylo_tree_display import ete_motifs

    server, db = manipulate_biosqldb.load_db(biodb)
    sets_form = blast_sets_form(biodb)

    if request.method == 'POST': 
        form = sets_form(request.POST)
        if form.is_valid():
            from chlamdb.plots import blast_heatmap
            #if request.method == 'POST': 

            hmm_sets = form.cleaned_data['blast_set']
            score_cutoff = form.cleaned_data['score_cutoff']
            query_coverage_cutoff = form.cleaned_data['query_coverage_cutoff']
            hit_coverage_cutoff = form.cleaned_data['hit_coverage_cutoff']

            counts = request.POST['counts']

            if counts == 'detailed':
                blast_set_filter = '"'+'","'.join(hmm_sets)+'"'

                gene2taxon2score, gene_list = blast_heatmap.get_multiple_set_profiles(biodb,
                                                                              hmm_sets,
                                                                              'bitscore',
                                                                              score_cutoff,
                                                                              query_coverage_cutoff,
                                                                              hit_coverage_cutoff)
                if hmm_sets == 'flagellum':
                    gene_list = ['Flg_fliE','Flg_flgC', 'Flg_flgB','Flg_sctJ_FLG','Flg_sctN_FLG','Flg_sctQ_FLG','Flg_sctR_FLG',
                                 'Flg_sctS_FLG','Flg_sctT_FLG','Flg_sctU_FLG','Flg_sctV_FLG']

                sql = 'select accession, t1.name from blast.blast_sets t1 inner join blast.blast_sets_entry t2 on t1.set_id=t2.set_id ' \
                      ' inner join blast.blast_db t3 on t2.seq_id=t3.seq_id where t1.name in (%s) order by t1.name;' % (blast_set_filter)

                ordered_gene_list = []
                for row in server.adaptor.execute_and_fetchall(sql):
                    ordered_gene_list.append('%s (%s)' % (row[0], row[1]))

                tree1, style1 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                                     ordered_gene_list,
                                                                     gene2taxon2score,
                                                                     identity_scale=True,
                                                                     show_labels=True,
                                                                     column_scale=True)

                path1 = settings.BASE_DIR + '/assets/temp/ortho_tree1.svg'
                asset_path1 = '/temp/ortho_tree1.svg'
                tree1.render(path1, dpi=800, tree_style=style1)


                gene2taxon2score2, gene_list2 = blast_heatmap.get_multiple_set_profiles(biodb,
                                                                              hmm_sets,
                                                                              'query_coverage',
                                                                              score_cutoff,
                                                                              query_coverage_cutoff,
                                                                              hit_coverage_cutoff)


                tree2, style2 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                            gene_list2,
                                                            gene2taxon2score2,
                                                            identity_scale=True,
                                                            show_labels=True,
                                                            column_scale=True)

                path2 = settings.BASE_DIR + '/assets/temp/ortho_tree2.svg'
                asset_path2 = '/temp/ortho_tree2.svg'
                tree2.render(path2, dpi=800, tree_style=style2)

                gene2taxon2score3, gene_list3 = blast_heatmap.get_multiple_set_profiles(biodb,
                                                                              hmm_sets,
                                                                              'locus_tag',
                                                                              score_cutoff,
                                                                              query_coverage_cutoff,
                                                                              hit_coverage_cutoff)


                tree3, style3 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                            gene_list3,
                                                            gene2taxon2score3,
                                                            identity_scale=False,
                                                            show_labels=True,
                                                            column_scale=False)

                path3 = settings.BASE_DIR + '/assets/temp/ortho_tree3.svg'
                asset_path3 = '/temp/ortho_tree3.svg'
                tree3.render(path3, dpi=800, tree_style=style3)

                gene2taxon2score4, gene_list4 = blast_heatmap.get_multiple_set_profiles(biodb,
                                                                              hmm_sets,
                                                                              'identity',
                                                                              score_cutoff,
                                                                              query_coverage_cutoff,
                                                                              hit_coverage_cutoff)


                tree4, style4 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                            gene_list4,
                                                            gene2taxon2score4,
                                                            identity_scale=True,
                                                            show_labels=True,
                                                            column_scale=True)

                path4 = settings.BASE_DIR + '/assets/temp/ortho_tree4.svg'
                asset_path4 = '/temp/ortho_tree4.svg'
                tree4.render(path4, dpi=800, tree_style=style4)


                gene2taxon2score5, gene_list5 = blast_heatmap.get_multiple_set_profiles(biodb,
                                                                              hmm_sets,
                                                                              'evalue',
                                                                              score_cutoff,
                                                                              query_coverage_cutoff,
                                                                              hit_coverage_cutoff)


                tree5, style5 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                            gene_list5,
                                                            gene2taxon2score5,
                                                            identity_scale=False,
                                                            show_labels=True,
                                                            column_scale=True,
                                                                     as_float=True)

                path5 = settings.BASE_DIR + '/assets/temp/ortho_tree5.svg'
                asset_path5 = '/temp/ortho_tree5.svg'
                tree5.render(path5, dpi=800, tree_style=style5)

                envoi = True
            else:
                from ete3 import Tree
                sql_tree = 'select tree from reference_phylogeny t1 inner join biodatabase t2 on t1.biodatabase_id=t2.biodatabase_id ' \
                           ' where t2.name="%s";' % biodb
                           
                server, db = manipulate_biosqldb.load_db(biodb)
                tree = server.adaptor.execute_and_fetchall(sql_tree,)[0][0]

                #acc_list = ['NC_010655',u'NC_013720',u'NZ_CP006571', u'NZ_AWUS01000000', u'NZ_APJW00000000', u'NC_002620', u'NZ_CCJF00000000', u'NZ_AYKJ01000000', u'NC_000117', u'LNES01000000', u'LJUH01000000', u'NC_004552', u'NC_003361', u'NC_007899', u'NC_015408', u'NC_000922', u'NC_015470', u'NZ_CCEJ000000000', u'CWGJ01000001', u'NZ_JSDQ00000000', u'NZ_BASK00000000', u'NZ_JRXI00000000', u'NZ_BAWW00000000', u'NZ_ACZE00000000', u'NC_015702', u'NZ_BBPT00000000', u'NZ_JSAN00000000', u'NC_005861', u'FCNU01000001', u'NZ_LN879502', u'NZ_BASL00000000', u'Rht', u'CCSC01000000', u'NC_015713', u'NC_014225']
                #filter = '"' + '","'.join(acc_list) + '"'
                #sql = 'select taxon_id from bioentry where biodatabase_id=102 and accession in (%s)' % filter
                #taxon_list = [str(i[0]) for i in server.adaptor.execute_and_fetchall(sql,)]
                t1 = Tree(tree)
                R = t1.get_midpoint_outgroup()
                t1.set_outgroup(R)
                #try:
                #    t2 = t1.prune(taxon_list)
                #except:
                #    pass
                t1.ladderize()

                set2taxon2count = blast_heatmap.get_multiple_set_counts(biodb,
                                                                        hmm_sets,
                                                                        score_cutoff,
                                                                        query_coverage_cutoff,
                                                                        hit_coverage_cutoff)

                tree1, style1 = ete_motifs.multiple_profiles_heatmap(biodb,
                                                                     hmm_sets,
                                                                     set2taxon2count,
                                                                     tree=t1,
                                                                     identity_scale=False,
                                                                     show_labels=True,
                                                                     column_scale=True,
                                                                     as_float=False,
                                                                     rotate=True)
                style1.rotation = 90
                path1 = settings.BASE_DIR + '/assets/temp/ortho_tree1.svg'
                asset_path1 = '/temp/ortho_tree1.svg'
                tree1.render(path1, dpi=800, tree_style=style1)
                envoi = True





    else:  
        form = sets_form()

    return render(request, 'chlamdb/blast_sets_profiles.html', my_locals(locals()))


def kegg_module_subcat(request):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)

    module_overview_form = make_module_overview_form(db, True)
    if request.method != "POST":
        form = module_overview_form()
        return render(request, 'chlamdb/module_subcat.html', my_locals(locals()))

    form = module_overview_form(request.POST)
    if not form.is_valid():
        # TODO: add error message
        form = module_overview_form()
        return render(request, 'chlamdb/module_subcat.html', my_locals(locals()))

    category        = form.cleaned_data["category"]
    leaf_to_name    = db.get_genomes_description()
    ko_count_subcat = db.get_ko_count_cat(subcategory=category)

    grouped_count = ko_count_subcat.groupby(["taxon_id", "module_id"]).sum()

    unique_module_ids = ko_count_subcat.index.get_level_values("module_id").unique().tolist()
    if len(unique_module_ids) == 0:
        # add error message : no module found
        envoi = True
        return render(request, 'chlamdb/module_subcat.html', my_locals(locals()))

    module_infos = db.get_modules_info(unique_module_ids)
    expression_tree = {}
    for module_id, descr, definition, *other in module_infos:
        parser = ModuleParser(definition)
        expression_tree[module_id] = parser.parse()

    labels = [format_module(val[0]) for val in module_infos]

    grouped_count = grouped_count.unstack(level=1, fill_value=0)
    grouped_count.columns = [col for col in grouped_count["count"].columns]
    ref_tree = db.get_reference_phylogeny()
    e_tree = EteTree.default_tree(ref_tree)
    for module, counts in grouped_count.iteritems():
        values = counts.to_dict()
        header = format_module(module)
        expr_tree = expression_tree[module]
        n_missing = {}
        for bioentry, count in counts.iteritems():
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
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)
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
    for module, counts in grouped_count.iteritems():
        values = counts.to_dict()
        header = format_module(module)
        expr_tree = expression_tree[module]
        n_missing = {}
        for bioentry, count in counts.iteritems():
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


def module_comparison(request):
    biodb = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db(biodb, settings.BIODB_CONF)
    comp_metabo_form = make_metabo_from(db)

    if request.method!="POST":
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
    envoi_comp=True
    return render(request, 'chlamdb/module_comp.html', my_locals(locals()))


def pfam_comparison(request):
    db = db_utils.DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    comp_metabo_form = make_metabo_from(db)
    if request.method != 'POST': 
        form = comp_metabo_form(request.POST)
        return render(request, 'chlamdb/pfam_comp.html', my_locals(locals()))

    form = comp_metabo_form(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/pfam_comp.html', my_locals(locals()))

    try:
        all_targets = form.get_choices()
    except:
        # TODO: add error message
        return render(request, 'chlamdb/pfam_comp.html', my_locals(locals()))

    genomes = db.get_genomes_description().description.to_dict()

    pfam_hits = db.get_pfam_hits(ids=all_targets)
    pfam_defs = db.get_pfam_def(pfam_hits.index.tolist(), add_ttl_count=True)
    data = []
    for pfam, entry_data in pfam_hits.iterrows():
        entry_infos = pfam_defs.loc[pfam]
        base_infos = [format_pfam(pfam, to_url=True), entry_infos["def"], entry_infos.ttl_cnt]
        data.append(base_infos+entry_data.values.tolist())
    ctx = {
        "envoi_comp": True,
        "taxon_list": pfam_hits.columns.tolist(),
        "pfam2data" : data,
        "taxon_id2description": genomes
    }
    return render(request, 'chlamdb/pfam_comp.html', my_locals(ctx))


def orthogroup_comparison(request):
    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb_path)

    comp_metabo_form = make_metabo_from(db)
    if request.method != 'POST': 
        form = comp_metabo_form(request.POST)
        return render(request, 'chlamdb/ortho_comp.html', my_locals(locals()))

    form = comp_metabo_form(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/ortho_comp.html', my_locals(locals()))

    try:
        all_targets = form.get_choices()
    except:
        # TODO: add error message
        return render(request, 'chlamdb/ortho_comp.html', my_locals(locals()))

    genomes = db.get_genomes_description().description.to_dict()
    og_count = db.get_og_count(all_targets)
    og_count.columns = [genomes[int(col)] for col in og_count.columns]
    annotations = db.get_genes_from_og(orthogroups=og_count.index.tolist(),
            taxon_ids=all_targets, terms=["product"])

    products = annotations.groupby("orthogroup")["product"].apply(list)
    n_orthogroups = len(og_count.index)

    og_data = []
    for og, items in og_count.iterrows():
        piece = [format_orthogroup(og)]
        piece.append(items.tolist())
        if og in products.index:
            piece.append(format_lst_to_html(products.loc[og]))
        else:
            piece.append("-")
        og_data.append(piece)

    genomes_list = og_count.columns.tolist()
    envoi_comp = True
    return render(request, 'chlamdb/ortho_comp.html', my_locals(locals()))


def ko_comparison(request):
    biodb_path = settings.BIODB_DB_PATH
    db = db_utils.DB.load_db_from_name(biodb_path)

    comp_metabo_form = make_metabo_from(db)

    if request.method != "POST":
        form = comp_metabo_form()
        return render(request, 'chlamdb/ko_comp.html', my_locals(locals()))

    form = comp_metabo_form(request.POST)
    if not form.is_valid():
        return render(request, 'chlamdb/ko_comp.html', my_locals(locals()))

    include = form.get_choices()
    mat_include = db.get_ko_count(include).unstack(level=0, fill_value=0)
    mat_include.columns = [col for col in mat_include["count"].columns.values]

    ko2annot = db.get_ko_desc(mat_include.index.tolist())
    df_ttl   = db.get_ko_count(mat_include.index.tolist(), search_on="ko_id")
    ko2total_count = df_ttl.groupby("KO").sum()["count"].to_dict()
    ko2counts = mat_include.to_dict()
    ko2counts = {}
    ko2_print = {}
    for key, values in mat_include.iterrows():
        ko2counts[key] = values.values.tolist()
        ko2_print[key] = format_ko(key)

    hsh_gen_desc = db.get_genomes_description().description.to_dict()
    taxon_list = [hsh_gen_desc[int(col)] for col in mat_include.columns.values]
    n_ko = len(mat_include.index.tolist())
    envoi_comp = True
    return render(request, 'chlamdb/ko_comp.html', my_locals(locals()))


def genomic_locus_tag_infos(request):
    db = db_utils.DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
    locus_tag = request.POST.get("locus_tag", None)
    
    if locus_tag is None:
        return JsonResponse({"error_msg": "no locus_tag"})
    seqid, is_pseudo = db.get_seqid(locus_tag=locus_tag, feature_type=False)
    prot_infos = db.get_proteins_info([seqid], as_df=True)
    entry = prot_infos.loc[seqid]
    data = {"product": entry["product"]}
    return JsonResponse(data)

