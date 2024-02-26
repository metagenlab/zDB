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
from chlamdb.forms import (make_blast_form, make_circos_form, make_metabo_from,
                           make_module_overview_form,
                           make_pathway_overview_form, make_plot_form,
                           make_single_genome_form, make_venn_from)
from django.conf import settings
from django.http import JsonResponse
from django.shortcuts import render
from django.views import View
from ete3 import SeqMotifFace, StackedBarFace, TextFace, Tree, TreeStyle
from lib import search_bar as sb
from lib.db_utils import DB, NoPhylogenyException
from lib.ete_phylo import (Column, EteTree, KOAndCompleteness,
                           ModuleCompletenessColumn, SimpleColorColumn)
from lib.KO_module import ModuleParser
from reportlab.lib import colors

from views.mixins import CogViewMixin, ComparisonViewMixin, KoViewMixin
from views.utils import (format_amr, format_cog, format_hmm_url, format_ko,
                         format_ko_modules, format_ko_path, format_locus,
                         format_orthogroup, format_pfam,
                         format_refseqid_to_ncbi, my_locals, optional2status,
                         page2title, to_s)


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

    hsh_files = db.get_filenames_to_taxon_id()
    number_of_files = len(hsh_files)

    number_ort = db.get_n_orthogroups()
    versions = db.get_versions_table()
    return render(request, 'chlamdb/home.html', my_locals(locals()))


class ComparisonIndexView(ComparisonViewMixin, View):

    @property
    def view_name(self):
        return f"index_comp_{self.object_type}"

    def get(self, request):
        context = my_locals({
            "page_title": self.page_title,
            "compared_obj_name": self.object_name_plural,
            "comp_type": self.object_type,
            "boxes": self.available_views,
            })
        return render(request, 'chlamdb/index_comp.html', context)


def get_genomes_data(db):
    genomes_data = db.get_genomes_infos()
    genomes_descr = db.get_genomes_description()

    genomes_data = genomes_data.join(genomes_descr)

    genomes_data.gc = genomes_data.gc.apply(lambda x: round(100 * x))
    genomes_data.coding_density = genomes_data.coding_density.apply(lambda x: round(100 * x))
    genomes_data.length = genomes_data.length.apply(lambda x: round(x / pow(10, 6), 2))
    return genomes_data


def genomes(request):
    biodb_path = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb_path)
    page_title = page2title["genomes"]
    genomes_data = get_genomes_data(db)

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
        "GC %",
        "N proteins",
        "N contigs",
        "Size (Mbp)",
        "Coding %",
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

    col_titles = {"closest_seq": "Closest Sequence",
                  "hmm_id": "HMM"}

    amr_hits["closest_seq"] = amr_hits["closest_seq"].map(format_refseqid_to_ncbi)
    amr_hits["gene"] = amr_hits["gene"].apply(format_amr, to_url=True)
    amr_hits["hmm_id"] = amr_hits["hmm_id"].apply(format_hmm_url)
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


def format_module(mod_id, base=None, to_url=False):
    if base is None:
        formated = f"M{mod_id:05d}"
    else:
        formated = base

    if to_url:
        return f"<a href=/KEGG_module_map/{formated}>{formated}</a>"
    return formated


class CogPhyloHeatmap(CogViewMixin, View):

    @property
    def view_name(self):
        return f"{self.object_type}_phylo_heatmap"

    def get_context(self, **kwargs):
        context = {
            "page_title": self.page_title,
        }
        context.update(kwargs)
        return my_locals(context)

    def post(self, request, *args, **kwargs):
        return render(request, 'chlamdb/cog_phylo_heatmap.html', self.get_context())

    def get(self, request, frequency, *args, **kwargs):
        freq = frequency != "False"

        tree = self.db.get_reference_phylogeny()
        descr = self.db.get_genomes_description()
        all_taxids = descr.index.tolist()

        all_cog_hits = self.db.get_cog_hits(all_taxids, search_on="taxid")
        all_cog_funcs = self.db.get_cog_summaries(all_cog_hits.index.unique().tolist(),
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
            col = SimpleColorColumn.fromSeries(func_count,
                                               header=detailed_func + "(" + func + ")", color_gradient=True)
            e_tree.add_column(col)

        freq = frequency
        path = settings.BASE_DIR + f"/assets/temp/COG_tree_{freq}.svg"
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


def js_bioentries_to_description(hsh):
    taxon_map = 'var taxon2description = { '
    mid = ",".join(f"{to_s(bioentry)}: {to_s(description)}" for bioentry, description in hsh.items())
    return taxon_map + mid + "};"


class KoBarchart(KoViewMixin, View):

    @property
    def view_name(self):
        return f"{self.object_type}_barchart"

    def get_context(self, **kwargs):
        context = {
            "page_title": self.page_title,
            "form": self.form,
        }
        context.update(kwargs)
        return my_locals(context)

    def get(self, request, *args, **kwargs):
        venn_form_class = make_venn_from(self.db)
        self.form = venn_form_class()
        return render(request, 'chlamdb/ko_barplot.html', self.get_context())

    def post(self, request, *args, **kwargs):
        venn_form_class = make_venn_from(self.db)
        self.form = venn_form_class(request.POST)
        if not self.form.is_valid():
            self.form = venn_form_class(request.POST)
            return render(request, 'chlamdb/ko_barplot.html', self.get_context())

        taxids = self.form.get_taxids()
        taxon2description = self.db.get_genomes_description().description.to_dict()

        ko_counts = self.db.get_ko_count(taxids, keep_seqids=True, as_multi=False)
        ko_ids = ko_counts.KO.unique()
        ko_module_ids = self.db.get_ko_modules(ko_ids.tolist(), as_pandas=True, compact=True)
        ko_modules_info = self.db.get_modules_info(ko_module_ids["module_id"].unique().tolist(), as_pandas=True)

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
        self.form = venn_form_class()
        context = self.get_context(
            envoi=True, series=series, labels=labels, taxids=taxids)
        return render(request, 'chlamdb/ko_barplot.html', context)


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


class CogBarchart(CogViewMixin, View):

    @property
    def view_name(self):
        return f"{self.object_type}_barchart"

    def get_context(self, **kwargs):
        context = {
            "page_title": self.page_title,
            "form": self.form,
        }
        context.update(kwargs)
        return my_locals(context)

    def get(self, request, *args, **kwargs):
        venn_form_class = make_venn_from(self.db)
        self.form = venn_form_class()
        return render(request, 'chlamdb/cog_barplot.html', self.get_context())

    def post(self, request, *args, **kwargs):
        venn_form_class = make_venn_from(self.db)
        self.form = venn_form_class(request.POST)
        if not self.form.is_valid():
            # add error message
            return render(request, 'chlamdb/cog_barplot.html', self.get_context())

        target_bioentries = self.form.get_taxids()

        hsh_counts = self.db.get_cog_counts_per_category(target_bioentries)
        taxon2description = self.db.get_genomes_description().description.to_dict()
        category_dico = self.db.get_cog_code_description()

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
        context = self.get_context(
            envoi=True, series=series, labels=labels, taxids=taxids,
            category_map=category_map, taxon_map=taxon_map)
        return render(request, 'chlamdb/cog_barplot.html', context)


class PanGenome(ComparisonViewMixin, View):

    @property
    def view_name(self):
        return f"pan_genome_{self.object_type}"

    def get_context(self, **kwargs):
        context = {
            "page_title": self.page_title,
            "form": self.form,
            "object_type": self.object_type,
            "object_name_plural": self.object_name_plural
        }
        context.update(kwargs)
        return my_locals(context)

    def get(self, request, *args, **kwargs):
        venn_form_class = make_venn_from(self.db, plasmid=False)
        self.form = venn_form_class()
        return render(request, 'chlamdb/pan_genome.html', self.get_context())

    def post(self, request, *args, **kwargs):
        venn_form_class = make_venn_from(self.db)
        self.form = venn_form_class(request.POST)
        if not self.form.is_valid():
            # add error message
            self.form = venn_form_class()
            return render(request, 'chlamdb/pan_genome.html', self.get_context())

        taxids = self.form.get_taxids()
        df_hits = self.get_hit_counts(taxids, search_on="taxid")

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

        context = self.get_context(
            envoi=True, js_data_acc=js_data_acc,
            js_data_count=js_data_count, js_data_core=js_data_core
            )
        return render(request, 'chlamdb/pan_genome.html', context)


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


class PlotHeatmap(ComparisonViewMixin, View):

    @property
    def view_name(self):
        return f"plot_heatmap_{self.object_type}"

    def get_context(self, **kwargs):
        context = {
            "page_title": self.page_title,
            "object_type": self.object_type,
            "form": self.form,
        }
        context.update(kwargs)
        return my_locals(context)

    def get(self, request, *args, **kwargs):
        venn_form_class = make_venn_from(self.db)
        self.form = venn_form_class()
        return render(request, 'chlamdb/plot_heatmap.html', self.get_context())

    def post(self, request, *args, **kwargs):
        venn_form_class = make_venn_from(self.db)
        self.form = venn_form_class(request.POST)
        if not self.form.is_valid():
            # add error message
            return render(request, 'chlamdb/plot_heatmap.html', self.get_context())

        import plotly.graph_objects as go
        import scipy.cluster.hierarchy as shc
        from scipy.cluster import hierarchy

        taxon_ids = self.form.get_taxids()

        if len(taxon_ids) <= 1:
            error_message = "Please select at least two genomes"
            error_title = "Wrong input"
            ctx = self.get_context(error=True, error_title=error_title,
                                   error_message=error_message)
            return render(request, 'chlamdb/plot_heatmap.html', ctx)

        mat = self.get_hit_counts(taxon_ids, search_on="taxid")
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
        context = self.get_context(envoi_heatmap=True, html_plot=html_plot)
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


def faq(request):
    a = 2
    return render(request, 'chlamdb/FAQ.html', my_locals(locals()))


def phylogeny(request):
    page_title = page2title["phylogeny_intro"]

    biodb_path = settings.BIODB_DB_PATH
    db = DB.load_db_from_name(biodb_path)

    genomes_data = get_genomes_data(db)

    asset_path = "/temp/species_tree.svg"
    path = settings.BASE_DIR + '/assets/temp/species_tree.svg'

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
        e_tree.add_column(SimpleColorColumn.fromSeries(data, header=None, use_col=False))
        stack = StackedBarColumn(data.to_dict(), colours=[col, "white"],
                                 relative=is_relative, header=header,
                                 header_params=header_params, face_params=stacked_face_params)
        e_tree.add_column(stack)

    e_tree.rename_leaves(genomes_data.description.to_dict())
    e_tree.render(path, dpi=500)
    return render(request, 'chlamdb/phylogeny_intro.html', my_locals(locals()))
