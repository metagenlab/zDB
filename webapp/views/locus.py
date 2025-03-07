import collections
import json
from io import StringIO
from wsgiref.util import FileWrapper

import matplotlib.colors as mpl_col
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
from django.conf import settings
from django.http import HttpResponse
from django.shortcuts import render
from django.views import View
from ete3 import SeqMotifFace
from ete3 import TextFace
from ete3 import Tree
from ete3 import TreeStyle
from lib.db_utils import DB
from lib.db_utils import NoPhylogenyException
from lib.ete_phylo import Column
from lib.ete_phylo import EteTree
from lib.ete_phylo import SimpleColorColumn
from views.mixins import AmrViewMixin
from views.mixins import CogViewMixin
from views.mixins import KoViewMixin
from views.mixins import PfamViewMixin
from views.mixins import VfViewMixin
from views.object_type_metadata import my_locals
from views.utils import DataTableConfig
from views.utils import format_gene
from views.utils import format_genome
from views.utils import format_locus
from views.utils import format_orthogroup
from views.utils import format_pfam
from views.utils import format_refseqid_to_ncbi
from views.utils import format_swissprot_entry
from views.utils import format_taxid_to_ncbi
from views.utils import genomic_region_df_to_js
from views.utils import locusx_genomic_region
from views.utils import make_div
from views.utils import optional2status


class DownloadSequences(View):
    def post(self, request, *args, **kwargs):
        db = DB.load_db(settings.BIODB_DB_PATH, settings.BIODB_CONF)
        data = json.loads(self.request.body)
        loci = data.get("loci", [])
        dna = data.get("dna", False)
        sequences = []
        for locus in loci:
            seqid, is_pseudogene = db.get_seqid(locus_tag=locus)
            hsh_infos = db.get_proteins_info(
                [seqid], to_return=["gene"], inc_non_CDS=True, inc_pseudo=True
            )
            organism = db.get_organism([seqid])[seqid]
            description = f"[organism={organism}]"
            gene = hsh_infos.get(seqid)
            if gene:
                description += f" {gene[0]}"
            if dna:
                l_flank, seq, r_flank = get_sequence(db, seqid)
            else:
                seq = db.get_translation(seqid)
            sequences.append(SeqRecord(seq, id=locus, description=description))
        out = StringIO()
        SeqIO.write(sequences, out, "fasta")
        out.seek(0)
        response = HttpResponse(
            FileWrapper(out), content_type="application/octet-stream; charset=utf-8"
        )
        response["Content-Disposition"] = 'attachment; filename="bla.fasta"'
        return response


def tab_general(db, seqid):
    hsh_infos = db.get_proteins_info(
        [seqid],
        to_return=["locus_tag", "gene", "product"],
        inc_non_CDS=True,
        inc_pseudo=True,
    )
    hsh_organism = db.get_organism([seqid], taxid_and_name=True)
    gene_loc = db.get_gene_loc([seqid], as_hash=False)

    gene_pos = []
    nucl_length = 0
    for index, row in gene_loc.iterrows():
        gene_pos.append((row.start, row.end, row.strand))
        nucl_length += row.end - row.start + 1

    locus_tag, gene, product = hsh_infos[seqid]
    if pd.isna(gene):
        gene = "-"
    taxid_and_organism = hsh_organism[seqid]
    return {
        "locus_tag": locus_tag,
        "organism": format_genome(taxid_and_organism),
        "gene_pos": gene_pos,
        "gene": gene,
        "nucl_length": nucl_length,
        "prot": product,
    }


def tab_contig(db, seqid):
    bioentry, accession, contig_size, _ = db.get_bioentry_list(seqid, search_on="seqid")
    qualifiers = db.get_bioentry_qualifiers(bioentry).set_index("term")["value"]
    return {
        "contig_size": contig_size,
        "contig_topology": qualifiers["topology"],
        "contig_is_plasmid": qualifiers["plasmid"] == "1",
        "contig_accession": accession,
    }


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

    e_tree.add_column(
        SimpleColorColumn.fromSeries(count.loc[group], header="Number of homologs")
    )
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
        col = SimpleColorColumn.fromSeries(
            max_identity.identity,
            color_gradient=True,
            header="Identity",
            default_val="-",
        )
        e_tree.add_column(col)

    e_tree.rename_leaves(leaf_to_name)

    dpi = 1200
    asset_path = f"/temp/og_conservation{group}.svg"
    path = settings.ASSET_ROOT + asset_path
    e_tree.render(path, dpi=dpi)
    return {"asset_path": asset_path}


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
    orga_set = set()
    for seqid, data in infos.iterrows():
        organism = hsh_organism[seqid]
        locus = format_locus(data.locus_tag, to_url=False)
        locus_fmt = format_locus(data.locus_tag, to_url=True)
        entry = [locus, locus_fmt, organism, format_gene(data.gene), data["product"]]
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

    n_genomes = (
        len(orga_set) if ref_seqid is not None else len(set(hsh_organism.values()))
    )
    return {
        "orthogroup": orthogroup_title,
        "n_genomes": "1 genome" if n_genomes == 1 else f"{n_genomes} genomes",
        "headers": headers,
        "homologues": homologues,
    }


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
        annots = db.get_genes_from_og(
            orthogroups=[og_id], terms=["locus_tag", "length"]
        )
        pfams = db.get_pfam_hits_info(annots.index.tolist())
        unique_pfams = pfams.pfam.unique()
        color_palette = (
            mpl_col.to_hex(col) for col in sns.color_palette(None, len(unique_pfams))
        )
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
        seqid_to_locus = db.get_proteins_info(
            seqids, to_return=["locus_tag"], as_df=True
        )
        all_infos = identity_matrix.join(seqid_to_locus).set_index("locus_tag").round(1)
        all_infos.loc[seqid_to_locus.loc[compare_to].locus_tag] = 100.0

        ident_col = SimpleColorColumn.fromSeries(
            all_infos.identity,
            color_gradient=True,
            header="Identity",
            default_val="-",
            is_str_index=True,
        )

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
    path = settings.ASSET_ROOT + asset_path
    e_tree.render(path, dpi=1200)

    algn_file = f"/alignments/{og_filename}"
    return {"og_phylogeny": asset_path, "root": root, "og_alignment": algn_file}


class AnnotationTableBase:
    def __init__(self, include_occurences=False, from_taxid=None):
        self.include_occurences = include_occurences
        self.from_taxid = from_taxid

    @property
    def table_data_accessors(self):
        accessors = super(AnnotationTableBase, self).table_data_accessors.copy()
        if self.include_occurences:
            accessors.insert(1, "occurences")
        return accessors

    def get_results(self, seqids):
        hits = self.get_hit_counts(seqids, search_on="seqid", indexing="seqid")
        if hits.empty:
            return {}

        descriptions = self.get_hit_descriptions(
            hits[self.object_column].to_list(),
            taxid_for_pathway_formatting=self.from_taxid,
        )

        if self.include_occurences:
            counts = hits[self.object_column].value_counts()
            counts.name = "occurences"
            descriptions = descriptions.merge(counts, left_index=True, right_index=True)

        return {
            "table_headers": self.table_headers,
            "table_data": descriptions,
            "table_data_accessors": self.table_data_accessors,
            "data_table_config": DataTableConfig(
                export_buttons=False, display_as_datatable=False
            ),
            "title": f"{self.object_name} Annotation(s)",
        }


class KoAnnotationTable(AnnotationTableBase, KoViewMixin):
    pass


class CogAnnotationTable(AnnotationTableBase, CogViewMixin):
    @property
    def table_data_accessors(self):
        accessors = super(CogAnnotationTable, self).table_data_accessors
        accessors.insert(2, "function")
        return accessors


class PfamAnnotationTable(AnnotationTableBase, PfamViewMixin):
    pass


class AmrAnnotationTable(AnnotationTableBase, AmrViewMixin):
    pass


class VfAnnotationTable(AnnotationTableBase, VfViewMixin):
    pass


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
            entry = [
                start,
                end,
                "[]",
                None,
                8,
                "black",
                self.pfam_cmap[pfam],
                fmt_entry,
            ]
            pfam_entries.append(entry)
        return SeqMotifFace(dummy_seq, motifs=pfam_entries, seq_format="line")


def og_tab_get_swissprot_homologs(db, annotations):
    homologs = db.get_swissprot_homologs(
        annotations.index.tolist(), indexing="accession"
    )
    summary = homologs.groupby("definition").count()
    return {
        "table_data": summary,
        "table_data_accessors": ["name", "accession"],
        "table_headers": ["Annotation", "Number of occurrences"],
        "data_table_config": DataTableConfig(table_id="swissprot", display_index=False),
        "title": f"Annotations of the {len(homologs)} SwissProt hits",
    }


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
        data = (
            "["
            + ",".join(f"{{x: {start}, y: {end}}}" for start, end in zip(starts, ends))
            + "]"
        )
        feature = (
            f'{{data: {data},  name: "{name}",   color: "#0F8292",  type : "rect",}}'
        )
        pfam_def = pfam_defs_df["def"].loc[pfam]
        pfam_defs.append((format_pfam(pfam, to_url=True), pfam_def))
        feature_viewer_fet.append(feature)
    return {
        "pfam_domains": "[" + ",".join(feature_viewer_fet) + "]",
        "pfam_def": pfam_defs,
    }


def locus_tab_swissprot_hits(db, seqid):
    swissprot_homologs = db.get_swissprot_homologs([seqid])
    header = [
        "Swissprot accession",
        "Eval",
        "Score",
        "ID (%)",
        "N gaps",
        "Alignment length",
        "Annot score",
        "Gene",
        "Description",
        "Organism",
    ]
    swissprot_homologs["ncbi"] = swissprot_homologs[["organism", "taxid"]].apply(
        format_taxid_to_ncbi, axis=1
    )
    swissprot_homologs["accession"] = swissprot_homologs["accession"].apply(
        format_swissprot_entry
    )

    accessors = [
        "accession",
        "evalue",
        "bitscore",
        "perc_id",
        "gaps",
        "match_len",
        "pe",
        "gene",
        "definition",
        "ncbi",
    ]
    return {
        "table_data": swissprot_homologs,
        "table_data_accessors": accessors,
        "table_headers": header,
        "data_table_config": DataTableConfig(table_id="swissprot", display_index=False),
        "n_swissprot_hits": len(swissprot_homologs),
    }


def tab_get_refseq_homologs(db, seqid):
    refseq_hits = db.get_refseq_hits([seqid]).set_index("match_id")
    refseq_hits_infos = db.get_refseq_matches_info(refseq_hits.index.tolist())
    all_infos = refseq_hits.join(refseq_hits_infos)

    header = [
        "Refseq accession",
        "Evalue",
        "Score",
        "ID(%)",
        "# gaps",
        "Len",
        "Description",
        "Organism",
    ]
    entries = []
    for match_id, data in all_infos.iterrows():
        to_ncbi = format_refseqid_to_ncbi(data.accession)
        entries.append(
            (
                to_ncbi,
                data.evalue,
                data.bitscore,
                data.pident,
                data.gaps,
                data.length,
                data.description,
                data.organism,
            )
        )
    return {
        "n_refseq_homologs": len(refseq_hits),
        "refseq_headers": header,
        "blast_data": entries,
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
        leaf.add_face(
            TextFace(f"{leaf.name} | {orga_name}", fgcolor=color), 0, "branch-right"
        )

    asset_path = f"/temp/og_best_hit_phylogeny_{orthogroup}.svg"
    path = settings.ASSET_ROOT + asset_path
    ts = TreeStyle()
    ts.show_leaf_name = False
    ete_tree.render(path, tree_style=ts, dpi=1200)
    return {"best_hits_phylogeny": asset_path, "has_refseq_phylo": True}


def get_sequence(db, seqid, flanking=0):
    loc = db.get_gene_loc([seqid], as_hash=False)
    bioentry, accession, length, seq = db.get_bioentry_list(seqid, search_on="seqid")

    if len(loc) == 2:
        # Need to handle the special case where a gene is overlapping both ends
        # of a circular contig.
        # This code assumes that the gene overlaps both ends of a circular
        # contig and won't work otherwise
        loc0 = loc.loc[0]
        loc1 = loc.loc[1]
        if loc0.strand != loc1.strand:
            raise Exception("Unsupported case of fragment gene on different strands")

        _, strand, start, stop = (int(i) for i in loc0.tolist())
        if start == 1:
            fet1 = SeqFeature(
                FeatureLocation(start - 1, stop + flanking, strand=strand)
            )
            fet0 = SeqFeature(
                FeatureLocation(
                    int(loc1.start - flanking - 1), int(loc1.end), strand=strand
                )
            )
        else:
            fet0 = SeqFeature(
                FeatureLocation(start - 1 - flanking, stop, strand=strand)
            )
            fet1 = SeqFeature(
                FeatureLocation(
                    int(loc1.start) - 1, int(loc1.end) + flanking, strand=strand
                )
            )
        extracted0 = fet0.extract(seq)
        extracted1 = fet1.extract(seq)
        extracted = extracted0 + extracted1
        red_start = flanking
        red_stop = len(extracted0) + (
            fet1.location.end - fet1.location.start - flanking
        )
    elif len(loc) == 1:
        _, strand, start, stop = (int(i) for i in loc.loc[0].tolist())
        start -= 1
        if start < flanking:
            start_w_flank = 0
            red_start = start
        else:
            start_w_flank = start - flanking
            red_start = flanking

        if stop + flanking > len(seq):
            stop_w_flank = len(seq) - 1
        else:
            stop_w_flank = stop + flanking
        red_stop = red_start + stop - start
        fet = SeqFeature(FeatureLocation(start_w_flank, stop_w_flank, strand=strand))
        extracted = fet.extract(seq)
    else:
        raise Exception("Unsupported case of fragmented gene")
    return (extracted[0:red_start], extracted[red_start:red_stop], extracted[red_stop:])


class ViewBase(View):
    locus = None
    seqid = None
    _db = None

    @property
    def db(self):
        if self._db is None:
            biodb = settings.BIODB_DB_PATH
            self._db = DB.load_db(biodb, settings.BIODB_CONF)
        return self._db

    def render_invalid(self, request, **kwargs):
        return render(request, self.template, my_locals({"valid_id": False, **kwargs}))

    def get(self, request, context):
        self.og_counts = self.db.get_og_count([self.og_id], search_on="orthogroup")

        self.og_annot = self.db.get_genes_from_og(
            orthogroups=[self.og_id], terms=["locus_tag", "gene", "product", "length"]
        )

        self.all_org = self.db.get_organism(self.og_annot.index.tolist())

        if self.show_homology_info:
            context.update(
                tab_og_conservation_tree(self.db, self.og_id, compare_to=self.seqid)
            )
            context.update(
                tab_homologs(
                    self.db,
                    self.og_annot,
                    self.all_org,
                    ref_seqid=self.seqid,
                    og=self.og_id,
                )
            )
            try:
                context.update(
                    tab_og_phylogeny(self.db, self.og_id, compare_to=self.seqid)
                )
            except NoPhylogenyException:
                pass
        else:
            context["n_genomes"] = "1 genome"

        context["result_tables"] = []
        for annotation_table in self.annotation_tables:
            if not annotation_table.is_enabled:
                continue
            context["result_tables"].append(annotation_table.get_results(self.seqids))

        if optional2status.get("BBH_phylogenies", False):
            context.update(tab_og_best_hits(self.db, self.og_id, locus=self.locus))
        return context


class LocusX(ViewBase):
    template = "chlamdb/locus.html"

    @property
    def annotation_tables(self):
        taxid = self.db.get_organism(self.seqids, as_taxid=True)[self.seqid]
        return [
            CogAnnotationTable(),
            KoAnnotationTable(from_taxid=taxid),
            AmrAnnotationTable(),
            VfAnnotationTable(),
        ]

    @property
    def seqids(self):
        return [self.seqid]

    def get(self, request, locus=None, menu=True):
        context = {"menu": True, "valid_id": True}
        if locus is None:
            return self.render_invalid(request)
        self.locus = locus
        try:
            self.seqid, feature_type, is_pseudogene = self.db.get_seqid(
                locus_tag=locus, feature_type=True
            )
        except Exception:
            return self.render_invalid(request)

        context["page_title"] = f"Locus tag: {locus}"
        l_flank, coding, r_flank = get_sequence(self.db, self.seqid, flanking=50)
        context["seq"] = l_flank + "<font color='red'>" + coding + "</font>" + r_flank
        context["feature_type"] = feature_type

        window_size = 8000
        all_infos, wd_start, wd_end, contig_size, contig_topology = (
            locusx_genomic_region(self.db, self.seqid, window=window_size)
        )
        context["genomic_region"] = genomic_region_df_to_js(
            all_infos, wd_start, wd_end, contig_size, contig_topology
        )
        context["window_size"] = window_size * 2
        context.update(tab_general(self.db, self.seqid))
        context.update(tab_contig(self.db, self.seqid))

        if feature_type != "CDS" or is_pseudogene:
            if is_pseudogene:
                context["feature_type"] = "Pseudogene"
            return render(request, self.template, my_locals(context))
        else:
            context["sequence_type"] = feature_type

        og_inf = self.db.get_og_count(self.seqids, search_on="seqid")
        # need to convert from numpy64 to int
        self.og_id = int(og_inf.loc[self.seqid].orthogroup)

        super(LocusX, self).get(request, context)

        # a bit of an hack
        translation = self.db.get_translation(self.seqid)
        context["length"] = len(translation)

        n_homologues = self.og_counts.loc[self.og_id].sum() - 1
        og_size = n_homologues + 1
        og_num_genomes = len(set(self.all_org.values()))

        if optional2status.get("pfam", False):
            context.update(tab_get_pfam_annot(self.db, self.seqids))

        if optional2status.get("BLAST_swissprot", False):
            context["swissprot"] = locus_tab_swissprot_hits(self.db, self.seqid)

        if optional2status.get("BLAST_database", False):
            context.update(tab_get_refseq_homologs(self.db, self.seqid))

        context.update(
            {
                "n_homologues": n_homologues,
                "og_id": format_orthogroup(self.og_id, to_url=True),
                "og_size": og_size,
                "og_num_genomes": og_num_genomes,
                "translation": translation,
                "locus": locus,
            }
        )
        return render(request, self.template, my_locals(context))

    def show_homology_info(self):
        return self.n_homologues > 1


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
        fig1.layout.margin.update(
            {
                "l": 80,
                "r": 20,
                "b": 40,
                "t": 20,
                "pad": 10,
            }
        )
        html_plot_prot_length = make_div(fig1, div_id="distplot")
    else:
        return {
            "length_distrib": True,
            "single_length": True,
            "prot_length": lengths.iloc[0],
        }

    return {
        "length_distrib": True,
        "max_protein_length": max_protein_length,
        "std_protein_length": std_protein_length,
        "min_protein_length": min_protein_length,
        "mean_protein_length": mean_protein_length,
        "median_protein_length": median_protein_length,
        "html_plot_prot_length": html_plot_prot_length,
    }


def format_lst(lst):
    hsh_values = {}
    for item in lst:
        val = hsh_values.get(item, 0)
        hsh_values[item] = val + 1
    return hsh_values


class Orthogroup(ViewBase):
    template = "chlamdb/og.html"

    @property
    def annotation_tables(self):
        return [
            CogAnnotationTable(include_occurences=True),
            KoAnnotationTable(include_occurences=True),
            PfamAnnotationTable(include_occurences=True),
            AmrAnnotationTable(include_occurences=True),
            VfAnnotationTable(include_occurences=True),
        ]

    @property
    def seqids(self):
        return self.og_annot.index.tolist()

    def get(self, request, og):
        context = {"valid_id": True}

        tokens = og.split("_")
        try:
            self.og_id = int(tokens[1])
        except Exception:
            return self.render_invalid(request, menu=True)

        super(Orthogroup, self).get(request, context)

        if len(self.og_counts.index) == 0:
            return self.render_invalid(request)

        hsh_genes = format_lst(self.og_annot["gene"].tolist())
        hsh_products = format_lst(self.og_annot["product"].tolist())

        n_homologues = self.og_counts.loc[self.og_id].sum()

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

        if optional2status.get("BLAST_swissprot", False):
            context["result_tables"].append(
                og_tab_get_swissprot_homologs(self.db, self.og_annot)
            )

        context.update(tab_lengths(n_homologues, self.og_annot))
        context.update(
            {
                "n_homologues": n_homologues,
                "og": og,
                "menu": True,
                "gene_annotations": gene_annotations,
                "product_annotations": product_annotations,
            }
        )
        return render(request, self.template, my_locals(context))

    def show_homology_info(self):
        return True
