
from django.conf import settings
from django.shortcuts import render
from django.views import View
from ete3 import Tree
from lib.db_utils import DB
from lib.ete_phylo import EteTree, SimpleColorColumn

from views.mixins import AmrViewMixin
from views.utils import (format_cog, format_ko, format_ko_module,
                         format_ko_path, format_orthogroup, my_locals,
                         page2title)


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


class FamAmrView(View, AmrViewMixin):

    template = 'chlamdb/fam.html'
    accessors = ["seq_name", "scope", "type", "class", "subclass", "hmm_id"]

    @property
    def view_name(self):
        return f"fam_{self.object_type}"

    def get(self, request, entry_id, *args, **kwargs):
        # Get hits for that entry:
        hit_counts = self.get_hit_counts(
            [entry_id], indexing="seqid", search_on=self.object_type,
            keep_taxid=True)

        if len(hit_counts) == 0:
            return render(request, self.template,
                          {"msg": f"No entry for {self.format_entry(entry_id)}"})

        seqids = hit_counts.index.tolist()

        orthogroups = self.db.get_og_count(seqids, search_on="seqid",
                                           keep_taxid=True)
        infos = self.get_hit_descriptions([entry_id])
        infos = infos.iloc[0]
        all_locus_data, group_count = get_all_prot_infos(
            self.db, seqids, orthogroups)

        hit_counts = hit_counts.groupby(["taxid"]).count()
        fam = self.format_entry(entry_id)
        e_tree = tab_gen_profile_tree(
            self.db, getattr(hit_counts, self.object_type),
            self.format_entry(entry_id), orthogroups)
        asset_path = f"/temp/fam_tree_{entry_id}.svg"
        path = settings.BASE_DIR + "/assets/" + asset_path
        e_tree.render(path, dpi=500)

        info = {self.colname_to_header[key]: infos[key]
                for key in self.accessors if infos[key]}

        context = {
            "page_title": self.page_title,
            "type": self.object_type,
            "fam": fam,
            "info": info,
            "all_locus_data": all_locus_data,
            "group_count": group_count,
            "asset_path": asset_path,
            "object_name": self.object_name,
            "object_name_plural": self.object_name_plural,
            "object_name_singular_or_plural": self.object_name_singular_or_plural,
        }
        return render(request, self.template, my_locals(context))


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
