import pandas as pd
from django.conf import settings
from lib.db_utils import DB


def safe_replace(string, search_string, replace_string):
    if string:
        return string.replace(search_string, replace_string)
    return string


title2page = {
    'Antimicrobial Resistance Gene': ["fam_amr"],
    'COG Ortholog': ['fam_cog'],
    'Comparisons: Antimicrobial Resistance': [
        'amr_comparison', 'index_comp_amr', 'entry_list_amr', 'extract_amr',
        'plot_heatmap_amr', 'venn_amr'],
    'Comparisons: Clusters of Orthologous groups (COGs)': [
        'cog_barchart', 'index_comp_cog', 'cog_phylo_heatmap',
        'cog_comparison', 'entry_list_cog', 'extract_cog', 'heatmap_COG',
        'pan_genome_cog', 'plot_heatmap_cog', 'venn_cog'],
    'Comparisons: Kegg Orthologs (KO)': [
        'entry_list_ko', 'extract_ko', 'heatmap_ko', 'index_comp_ko',
        'ko_comparison', 'ko_barchart', 'pan_genome_ko', 'plot_heatmap_ko',
        'venn_ko'],
    'Comparisons: PFAM domains': [
        'entry_list_pfam', 'extract_pfam', 'heatmap_pfam', 'pan_genome_pfam',
        'pfam_comparison', 'venn_pfam', 'index_comp_pfam', 'plot_heatmap_pfam'],
    'Comparisons: orthologous groups': [
        'extract_orthogroup', 'heatmap_orthogroup', 'index_comp_orthogroup',
        'orthogroup_comparison', 'pan_genome_orthogroup', 'plot_heatmap_orthogroup',
        'venn_orthogroup'],
    'Genome alignments: Circos plot': ['circos'],
    'Genome alignments: Plot region': ['plot_region'],
    'Genomes: table of contents': ['extract_contigs', 'genomes_intro'],
    'Homology search: Blast': ['blast'],
    'Kegg Ortholog': ['fam_ko'],
    'Kegg metabolic pathways': ['kegg_genomes'],
    'Kegg module': ['KEGG_module_map'],
    'Metabolism: kegg based': [
        'KEGG_mapp_ko', 'kegg', 'kegg_genomes_modules', 'kegg_module',
        'kegg_module_subcat', 'module_comparison'],
    'Pfam domain': ['fam_pfam'],
    'Phylogeny': ['phylogeny_intro'],
    'Table of content: Genomes': ['genomes'],
}

page2title = {}
for value, keys in title2page.items():
    page2title.update({key: value for key in keys})


# could also be extended to cache the results of frequent queries
# (e.g. taxid -> organism name) to avoid db queries
with DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF) as db:
    hsh_config = db.get_config_table(ret_mandatory=True)
    optional2status = {name: value for name,
                       (mandatory, value) in hsh_config.items() if not mandatory}
    missing_mandatory = [name for name, (mandatory, value) in hsh_config.items()
                         if mandatory and not value]


def my_locals(local_dico):
    local_dico["optional2status"] = optional2status
    local_dico["missing_mandatory"] = missing_mandatory
    return local_dico


def to_s(f):
    return "\"" + str(f) + "\""


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


def format_cog(cog_id, as_url=False, base=None):
    if base is None:
        base = f"COG{int(cog_id):04d}"
    if as_url is False:
        return base
    return f"<a href=\"/fam_cog/{base}\">{base}</a>"


def format_cog_url(cog_id):
    return format_cog(cog_id, as_url=True)


def format_ko(ko_id, as_url=False, base=None):
    if base is None:
        base = f"K{int(ko_id):05d}"
    if not as_url:
        return base
    return f"<a href=\"/fam_ko/{base}\">{base}</a>"


def format_ko_url(ko_id):
    return format_ko(ko_id, as_url=True)


def format_amr(gene, to_url=False):
    if not to_url:
        return gene
    return f"<a href=\"/fam_amr/{gene}\">{gene}</a>"


def format_hmm_url(hmm_id):
    if hmm_id:
        hmm_id = hmm_id.rsplit(".", 1)[0]
        return f"<a href=\"https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/{hmm_id}\">{hmm_id}</a>"  # noqa
    return hmm_id


def format_pfam(pfam_id, base=None, to_url=False):
    if base is None:
        fmt_entry = f"PF{pfam_id:04d}"
    else:
        fmt_entry = base
    if to_url:
        return f"<a href=/fam_pfam/{fmt_entry}>{fmt_entry}</a>"
    return fmt_entry


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
