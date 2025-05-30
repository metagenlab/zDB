import pandas as pd
from django.conf import settings
from django.shortcuts import render
from django.views import View
from ete3 import Tree
from lib.ete_phylo import EteTree
from lib.ete_phylo import SimpleColorColumn
from views.mixins import AmrViewMixin
from views.mixins import CogViewMixin
from views.mixins import GiViewMixin
from views.mixins import KoViewMixin
from views.mixins import OrthogroupViewMixin
from views.mixins import PfamViewMixin
from views.mixins import VfViewMixin
from views.utils import ResultTab
from views.utils import TabularResultTab
from views.utils import format_ko_module
from views.utils import format_ko_path
from views.utils import format_locus
from views.utils import format_orthogroup
from views.utils import locusx_genomic_region
from views.utils import optional2status
from views.views import prepare_genomic_regions


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

    table_headers = [
        "Orthogroup",
        "Locus",
        "Start",
        "Stop",
        "S.",
        "Gene",
        "Product",
        "Organism",
    ]

    table_accessors = table_headers
    tabular_result_tab_header = "Protein list"
    hit_count_indexing = "seqid"

    external_link_tag = (
        '<a href="{}" target="_blank"> {} <i class="fas fa-external-link-alt"></i></a>'
    )

    @property
    def profile_tab_help_text(self):
        return (
            f"phylogenetic tree annotated with"
            f"<br>- the presence of the {self.object_name} of interest within all "
            f"the genomes of the database (first column)"
            f"<br>- the size of the orthogroup(s) in which the reported {self.object_name} has been "
            f"clustered."
            f'<br>Orthogroups are shown in <font size="2" color="red"> red if they are annotated with that '
            f'particular {self.object_name} in a given genome </font>, in <font size="2" color="green">'
            f"green otherwise, i.e. if there is a discrepency between orthogroup clustering and {self.object_name}"
            f" prediction</font>."
            f"<br>Green homologs (same orthogroup) <strong>are not</strong> positive hit(s) for the "
            f"considered {self.object_name}."
            f"<br>Variations within orthogroups may be due to the clustering of multi domain proteins"
            f" or because of erroneous homolog clustering or {self.object_name} prediction."
        )

    @property
    def help_text(self):
        text = (
            f"{'Three' if len(self.tabs) == 3 else 'Four'} outputs have been generated:"
        )
        for tab in self.tabs:
            text += f"<br> <b>{tab.title}</b>: {tab.help_text}"
        return text

    @property
    def view_name(self):
        return f"fam_{self.object_type}"

    def get(self, request, entry_id, *args, **kwargs):
        context = self.prepare_context(request, entry_id, *args, **kwargs)
        return render(request, self.template, context)

    def get_orthogroups(self, seqids):
        return self.db.get_og_count(seqids, search_on="seqid", keep_taxid=True)

    def get_profile_tree(self, main_series, header):
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
        return e_tree

    def add_additional_columns(self, e_tree):
        # the (group, taxid) in this dataframe are those that should be colored in red
        # in the profile (correspondance between a cog entry and an orthogroup)
        unique_og = self.orthogroups.orthogroup.unique().tolist()
        red_color = set(tuple(entry) for entry in self.orthogroups.to_numpy())
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

    def get_table(self, seqids):
        hsh_gene_locs = self.db.get_gene_loc(seqids)
        hsh_prot_infos = self.db.get_proteins_info(seqids)
        hsh_organisms = self.db.get_organism(seqids)
        all_locus_data = []

        for seqid in seqids:
            # NOTE: all seqids are attributed an orthogroup, the case where
            # seqid is not in orthogroups should therefore not arise.
            og = self.orthogroups.loc[seqid].orthogroup
            fmt_orthogroup = format_orthogroup(og, to_url=True)
            strand, start, end = hsh_gene_locs[seqid]
            organism = hsh_organisms[seqid]
            locus, prot_id, gene, product = hsh_prot_infos[seqid]
            if gene is None:
                gene = ""
            data = (
                fmt_orthogroup,
                locus,
                start,
                end,
                strand,
                gene,
                product,
                organism,
            )
            all_locus_data.append(data)

        table_data = pd.DataFrame(all_locus_data, columns=self.table_headers)
        table_data["Locus"] = table_data["Locus"].apply(format_locus)
        return table_data, self.table_headers, self.table_accessors

    def get_associated_entries(self, table_data):
        return table_data["Orthogroup"].unique()

    def maybe_add_external_link(self, info):
        pass

    def prepare_context(self, request, entry_id, *args, **kwargs):
        # Get hits for that entry:
        hit_counts = self.get_hit_counts(
            [entry_id],
            indexing=self.hit_count_indexing,
            search_on=self.object_type,
            keep_taxid=True,
        )

        if len(hit_counts) == 0:
            return render(
                request,
                self.template,
                {"msg": f"No entry for {self.format_entry(entry_id)}"},
            )

        if hit_counts.index.name in ["seqid", "gi"]:
            seqids = hit_counts.index.tolist()
        else:
            # Pfam hits are not indexed with seqid...
            seqids = hit_counts.seqid.unique().tolist()

        self.orthogroups = self.get_orthogroups(seqids)
        infos = self.get_hit_descriptions(
            [entry_id], columns=self.accessors, extended_data=False
        )
        infos = infos.iloc[0]

        hit_counts = hit_counts.groupby(["taxid"]).count()
        self.fam = self.format_entry(entry_id)
        e_tree = self.get_profile_tree(
            getattr(hit_counts, self.object_column),
            self.format_entry(entry_id),
        )

        self.add_additional_columns(e_tree)
        self.asset_path = f"/temp/fam_tree_{entry_id}.svg"
        path = settings.ASSET_ROOT + self.asset_path
        e_tree.render(path, dpi=500)

        info = {}
        self.maybe_add_external_link(info)
        info.update(
            {
                self.colname_to_header(key): infos[key]
                for key in self.accessors
                if infos[key]
            }
        )

        table_data, table_headers, table_accessors = self.get_table(seqids)
        self.table_size = len(table_data)
        self.tabs = self.result_tabs(table_data, table_headers, table_accessors)
        context = self.get_context(
            fam=self.fam,
            info=info,
            table_size=self.table_size,
            group_count=self.get_associated_entries(table_data),
            object_name_singular_or_plural=self.object_name_singular_or_plural,
            help_text=self.help_text,
            result_tabs=self.tabs,
        )
        return context

    def result_tabs(self, table_data, table_headers, table_accessors):
        return [
            ResultTab(
                "general",
                "General",
                "chlamdb/fam_general_tab.html",
                help_text=f"this tab contains the description, occurence in the database "
                f"and other information related to the selected {self.object_name} {self.fam}",
            ),
            TabularResultTab(
                "distribution",
                self.tabular_result_tab_header,
                table_headers=table_headers,
                table_data=table_data,
                table_data_accessors=table_accessors,
                selectable=True,
                help_text=f"lists the {self.table_size} occurences of the {self.object_name} within the database. "
                "The table reports information on each occurence.",
            ),
            ResultTab(
                "profile",
                "Profile",
                "chlamdb/result_asset.html",
                asset_path=self.asset_path,
                help_text=self.profile_tab_help_text,
            ),
        ]


class FamAmrView(FamBaseView, AmrViewMixin):
    accessors = ["seq_name", "scope", "type", "class", "subclass", "hmm_id"]

    def maybe_add_external_link(self, info):
        info["Gene"] = self.external_link_tag.format(
            f"https://www.ncbi.nlm.nih.gov/pathogens/refgene/#virulence_genotypes:{self.fam}",
            self.fam,
        )


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

    def maybe_add_external_link(self, info):
        info["External link"] = self.external_link_tag.format(
            f"http://www.genome.jp/dbget-bin/www_bget?{self.fam}",
            self.fam,
        )

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

    def maybe_add_external_link(self, info):
        info["External link"] = self.external_link_tag.format(
            f"https://www.ebi.ac.uk/interpro/entry/pfam/{self.fam}",
            self.fam,
        )

    def get(self, request, entry_id, *args, **kwargs):
        entry_id = int(entry_id[len("PF") :])
        return super(FamPfamView, self).get(request, entry_id, *args, **kwargs)


class FamGiClusterView(FamBaseView, GiViewMixin):
    accessors = ["cluster_id", "length"]
    hit_count_indexing = "gi"
    tabular_result_tab_header = "GIs list"

    @property
    def profile_tab_help_text(self):
        return (
            f"phylogenetic tree annotated with"
            f"<br>- the presence of the {self.object_name} of interest within all "
            f"the genomes of the database (first column)"
            f"<br>- the size of the orthogroup(s) associated with the {self.object_name} in each genome."
            f'<br>Orthogroups are shown in <font size="2" color="red"> red if they are found withinin the '
            f'that particular {self.object_name} in a given genome </font>, in <font size="2" color="green">'
            f"green otherwise, i.e. if the orthogroup is present in the genome but not in that "
            f"{self.object_name}</font>."
            f"<br>Green homologs (same orthogroup) <strong>are not</strong> found inside the "
            f"considered {self.object_name}."
            f"<br>Variations within orthogroups may be due to the clustering of multi domain proteins"
            f" or because of erroneous homolog clustering or {self.object_name} prediction."
        )

    def get(self, request, entry_id, *args, **kwargs):
        entry_id = int(entry_id[3:])
        return super(FamGiClusterView, self).get(request, entry_id, *args, **kwargs)

    def get_orthogroups(self, gis_ids):
        self.gi_descriptions = self.get_gi_descriptions(gis_ids)
        self.gics = self.gi_descriptions[
            ["bioentry.bioentry_id", "start_pos", "end_pos"]
        ]
        seqids = set()
        for gis_id, row in self.gics.iterrows():
            bioentry = int(row["bioentry.bioentry_id"])
            start_pos = int(row["start_pos"])
            end_pos = int(row["end_pos"])
            df_seqids = self.db.get_features_location(bioentry, search_on="bioentry_id")
            df_seqids = df_seqids[
                (df_seqids.end_pos > start_pos) & (df_seqids.start_pos < end_pos)
            ]
            seqids.update(df_seqids.seqfeature_id)
        self.seqids = list(seqids)
        return super(FamGiClusterView, self).get_orthogroups(self.seqids)

    def get_table(self, gis_ids):
        table_data = self.gi_descriptions.drop(
            columns=["bioentry.bioentry_id", "cluster_id"]
        )
        return (
            table_data,
            [self.colname_to_header(colname) for colname in table_data.columns],
            table_data.columns,
        )

    def get_associated_entries(self, table_data):
        return table_data["gis_id"].unique()

    def add_genomic_region(self, tabs):
        genomic_regions = []
        for gis_id, row in self.gics.iterrows():
            genomic_regions.append(
                list(
                    locusx_genomic_region(
                        self.db,
                        bioentry=int(row["bioentry.bioentry_id"]),
                        window_start=int(row["start_pos"]),
                        window_stop=int(row["end_pos"]),
                    )
                )
            )
        all_regions, connections, all_identities = prepare_genomic_regions(
            self.db, genomic_regions, allow_flips=True
        )
        to_highlight = {}
        if optional2status.get("vf", False):
            vfs = self.db.vf.get_hits_from_seqids(self.seqids, columns=("seqid",))
            to_highlight.update(
                {
                    el: "purple"
                    for el in self.db.get_proteins_info(
                        vfs.seqid.to_list(), to_return=["locus_tag"], as_df=True
                    ).get("locus_tag", [])
                }
            )
        if optional2status.get("amr", False):
            amrs = self.db.get_amr_hits_from_seqids(self.seqids, columns=("seqid",))
            to_highlight.update(
                {
                    el: "magenta"
                    for el in self.db.get_proteins_info(
                        amrs.seqid.to_list(), to_return=["locus_tag"], as_df=True
                    ).get("locus_tag", [])
                }
            )

        tabs[0].show_genomic_region = True
        tabs[0].genomic_regions = "[" + "\n,".join(all_regions) + "]"
        tabs[0].connections = "[" + ",".join(connections) + "]"
        tabs[0].to_highlight = to_highlight
        if len(all_identities) > 0:
            tabs[0].max_ident = max(all_identities)
            tabs[0].min_ident = min(all_identities)
        tabs[0].description = (
            "This plot shows the genomic regions of the selected GIs."
            "<br>Genes in magenta are AMR genes, while purple correspond to VF genes."
        )

    def add_orthogroup_table(self, tabs):
        og_mixin = OrthogroupViewMixin()
        og_table = og_mixin.get_hit_descriptions(list(set(self.orthogroups.orthogroup)))
        tabs.insert(
            -1,
            TabularResultTab(
                "orthogroups",
                "Orthogroup list",
                table_data=og_table,
                table_data_accessors=og_table.columns,
                table_headers=og_mixin.table_headers,
                selectable=True,
                help_text=f"lists the {len(og_table)} orthogroups associated with the {self.object_name}",
            ),
        )

    def result_tabs(self, table_data, table_headers, table_accessors):
        tabs = super(FamGiClusterView, self).result_tabs(
            table_data, table_headers, table_accessors
        )
        self.add_genomic_region(tabs)
        self.add_orthogroup_table(tabs)
        return tabs
