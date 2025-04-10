from django.conf import settings
from django.shortcuts import render
from django.views import View
from ete3 import Tree
from lib.ete_phylo import EteTree
from lib.ete_phylo import SimpleColorColumn
from views.mixins import AmrViewMixin
from views.mixins import CogViewMixin
from views.mixins import KoViewMixin
from views.mixins import PfamViewMixin
from views.mixins import VfViewMixin
from views.utils import format_ko_module
from views.utils import format_ko_path
from views.utils import format_orthogroup


class FamCogColorFunc:
    def __init__(self, og, red_color):
        self.og = og
        self.red_color = red_color

    def get_color(self, taxid):
        if (self.og, taxid) in self.red_color:
            return "#FA5858"
        else:
            return EteTree.GREEN


class FamBaseView(View):
    template = "chlamdb/fam.html"

    @property
    def view_name(self):
        return f"fam_{self.object_type}"

    def get(self, request, entry_id, *args, **kwargs):
        context = self.prepare_context(request, entry_id, *args, **kwargs)
        return render(request, self.template, context)

    def get_orthogroups(self, seqids):
        return self.db.get_og_count(seqids, search_on="seqid", keep_taxid=True)

    def get_profile_tree(self, main_series, header, intersect):
        """
        Generate the tree from the profiles tab in the pfam/ko/cog pages:
        -ref_tree: the phylogenetic tree
        - main_series: the cog/ko/pfam count per taxid
        - header: the header of the main_series in the tree
        -intersect: a dataframe containing the seqid, taxid and orthogroups of the pfam/cog/ko hits
        """
        ref_tree = self.db.get_reference_phylogeny()
        ref_names = self.db.get_genomes_description().description.to_dict()

        tree = Tree(ref_tree)
        R = tree.get_midpoint_outgroup()
        if R is not None:
            tree.set_outgroup(R)
        tree.ladderize()
        e_tree = EteTree(tree)
        e_tree.rename_leaves(ref_names)

        e_tree.add_column(SimpleColorColumn.fromSeries(main_series, header=header))
        self.add_additional_columns(e_tree, intersect)
        return e_tree

    def add_additional_columns(self, e_tree, intersect):
        # the (group, taxid) in this dataframe are those that should be colored in red
        # in the profile (correspondance between a cog entry and an orthogroup)
        unique_og = intersect.orthogroup.unique().tolist()
        red_color = set(tuple(entry) for entry in intersect.to_numpy())
        df_og_count = self.db.get_og_count(list(unique_og), search_on="orthogroup").T
        for og in df_og_count:
            og_serie = df_og_count[og]
            color_chooser = FamCogColorFunc(og, red_color)
            col_column = SimpleColorColumn(
                og_serie.to_dict(),
                header=format_orthogroup(og),
                col_func=color_chooser.get_color,
            )
            e_tree.add_column(col_column)

    def get_all_prot_infos(self, seqids, orthogroups):
        hsh_gene_locs = self.db.get_gene_loc(seqids)
        hsh_prot_infos = self.db.get_proteins_info(seqids)
        hsh_organisms = self.db.get_organism(seqids)
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
            data = (
                index,
                fmt_orthogroup,
                locus,
                prot_id,
                start,
                end,
                strand,
                gene,
                product,
                organism,
            )
            all_locus_data.append(data)
        return all_locus_data, group_count

    def prepare_context(self, request, entry_id, *args, **kwargs):
        # Get hits for that entry:
        hit_counts = self.get_hit_counts(
            [entry_id], indexing="seqid", search_on=self.object_type, keep_taxid=True
        )

        if len(hit_counts) == 0:
            return render(
                request,
                self.template,
                {"msg": f"No entry for {self.format_entry(entry_id)}"},
            )

        if hit_counts.index.name == "seqid":
            seqids = hit_counts.index.tolist()
        else:
            # Pfam hits are not indexed with seqid...
            seqids = hit_counts.seqid.unique().tolist()

        orthogroups = self.get_orthogroups(seqids)
        infos = self.get_hit_descriptions(
            [entry_id], columns=self.accessors, extended_data=False
        )
        infos = infos.iloc[0]
        all_locus_data, group_count = self.get_all_prot_infos(seqids, orthogroups)

        hit_counts = hit_counts.groupby(["taxid"]).count()
        fam = self.format_entry(entry_id)
        e_tree = self.get_profile_tree(
            getattr(hit_counts, self.object_column),
            self.format_entry(entry_id),
            orthogroups,
        )
        asset_path = f"/temp/fam_tree_{entry_id}.svg"
        path = settings.ASSET_ROOT + asset_path
        e_tree.render(path, dpi=500)

        info = {
            self.colname_to_header(key): infos[key]
            for key in self.accessors
            if infos[key]
        }

        context = self.get_context(
            fam=fam,
            info=info,
            all_locus_data=all_locus_data,
            group_count=group_count,
            asset_path=asset_path,
            object_name_singular_or_plural=self.object_name_singular_or_plural,
        )
        return context


class FamAmrView(FamBaseView, AmrViewMixin):
    accessors = ["seq_name", "scope", "type", "class", "subclass", "hmm_id"]


class FamVfView(FamBaseView, VfViewMixin):
    accessors = [
        "prot_name",
        "vfid",
        "category",
        "gb_accession",
        "characteristics",
        "structure",
        "function",
        "mechanism",
    ]


class FamCogView(FamBaseView, CogViewMixin):
    accessors = ["description", "function_descr"]

    def get(self, request, entry_id, *args, **kwargs):
        entry_id = int(entry_id[3:])
        return super(FamCogView, self).get(request, entry_id, *args, **kwargs)


class FamKoView(FamBaseView, KoViewMixin):
    accessors = ["ko", "description"]

    def get(self, request, entry_id, *args, **kwargs):
        entry_id = int(entry_id[len("K") :])
        return super(FamKoView, self).get(request, entry_id, *args, **kwargs)

    def prepare_context(self, request, entry_id, *args, **kwargs):
        pathways = self.db.get_ko_pathways([entry_id])
        modules = self.db.get_ko_modules([entry_id])
        modules_id = [
            mod_id for key, values in modules.items() for mod_id, desc in values
        ]
        modules_data = self.db.get_modules_info(modules_id)
        pathway_data = format_ko_path(pathways, entry_id, as_list=True)
        module_data = [
            (format_ko_module(mod_id), cat, mod_desc)
            for mod_id, mod_desc, mod_def, path, cat in modules_data
        ]

        context = super(FamKoView, self).prepare_context(
            request, entry_id, *args, **kwargs
        )
        context["pathway_data"] = pathway_data
        context["module_data"] = module_data
        return context


class FamPfamView(FamBaseView, PfamViewMixin):
    accessors = ["def"]

    def get(self, request, entry_id, *args, **kwargs):
        entry_id = int(entry_id[len("PF") :])
        return super(FamPfamView, self).get(request, entry_id, *args, **kwargs)
