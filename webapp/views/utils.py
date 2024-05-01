import json
import re

import pandas as pd
from django.conf import settings
from django.core.serializers.json import DjangoJSONEncoder
from lib.db_utils import DB


def safe_replace(string, search_string, replace_string):
    if string:
        return string.replace(search_string, replace_string)
    return string


title2page = {
    'Antimicrobial Resistance Gene': ["fam_amr"],
    'COG Ortholog': ['fam_cog'],
    'Comparisons: Antimicrobial Resistance': ['amr_comparison'],
    'Comparisons: Clusters of Orthologous groups (COGs)': ['cog_comparison'],
    'Comparisons: Kegg Orthologs (KO)': ['ko_comparison'],
    'Comparisons: PFAM domains': ['pfam_comparison'],
    'Comparisons: orthologous groups': ['orthogroup_comparison'],
    'Comparisons: Virulence Factors': ['vf_comparison'],
    'Genome alignments: Circos plot': ['circos'],
    'Genome alignments: Plot region': ['plot_region'],
    'Genomes: table of contents': ['extract_contigs'],
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
    'Virulence Factor Gene': ["fam_vf"],
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
    optional2status["cog"] = optional2status["COG"]
    optional2status["ko"] = optional2status["KEGG"]
    optional2status["module"] = optional2status["KEGG"]
    optional2status["pathway"] = optional2status["KEGG"]
    optional2status["amr"] = optional2status["AMR"]
    optional2status["vf"] = optional2status["BLAST_vfdb"]

    missing_mandatory = [name for name, (mandatory, value) in hsh_config.items()
                         if mandatory and not value]


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


icon_external_link = '<i class="fas fa-external-link-alt"></i>'


def format_hmm_url(hmm_id):
    if hmm_id:
        hmm_id = hmm_id.rsplit(".", 1)[0]
        return f'<a href="https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/{hmm_id}"'\
               f' target="_blank">{hmm_id}&nbsp{icon_external_link}</a>'  # noqa
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


def format_refseqid_to_ncbi(seqid):
    return f"<a href=\"http://www.ncbi.nlm.nih.gov/protein/{seqid}\">{seqid}</a>"


def format_gene(gene):
    if pd.isna(gene):
        return "-"
    else:
        return gene


def format_taxid_to_ncbi(organism_and_taxid):
    organism, taxid = organism_and_taxid
    val = (
        f"""<a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}" target="_top">"""
        f"""{organism}</a>"""
    )
    return val


def format_swissprot_entry(entry_id):
    return f'<a href="https://www.uniprot.org/uniprot/{entry_id}">{entry_id}</a>'


class DataTableConfig():

    def __init__(self, table_id="results", ordering=True, paging=True,
                 export_buttons=True, colvis_button=False, display_index=False,
                 display_as_datatable=True):
        self.table_id = table_id
        self.ordering = ordering
        self.paging = paging
        self.export_buttons = export_buttons
        self.colvis_button = colvis_button
        self.display_index = display_index
        self.display_as_datatable = display_as_datatable
        if self.display_as_datatable:
            self.style = "margin-top: 3em;"
        else:
            self.style = ""

    @property
    def buttons(self):
        buttons = []
        if self.colvis_button:
            buttons.append({"extend": 'colvis', "columns": ':not(.noVis)'})
        if self.export_buttons:
            buttons.extend([
                {"extend": 'excel', "title": self.table_id},
                {"extend": 'csv', "title": self.table_id}
            ])
        return buttons

    @property
    def dom(self):
        if self.export_buttons or self.colvis_button:
            return 'lBfrtip'
        return 'lfrtip'

    def to_json(self):
        return json.dumps(
            {
                "paging": self.paging,
                "ordering": self.ordering,
                "info": False,
                "buttons": self.buttons,
                "dom": self.dom,
                "display_as_datatable": self.display_as_datatable,
            },
            cls=DjangoJSONEncoder)


class ResultTab():

    def __init__(self, tabid, title, template, show_badge=False, badge=None, **kwargs):
        self.id = tabid
        self.title = title
        self.template = template
        self.show_badge = show_badge
        self.badge = badge
        for key, val in kwargs.items():
            setattr(self, key, val)


class TabularResultTab(ResultTab):

    def __init__(self, tabid, title, template="chlamdb/result_table.html",
                 ordering=True, paging=True, export_buttons=True,
                 colvis_button=False, display_index=False,
                 show_badge=False, **kwargs):
        self.data_table_config = DataTableConfig(
            table_id=tabid, ordering=ordering, paging=paging,
            export_buttons=export_buttons, colvis_button=colvis_button,
            display_index=display_index)
        if show_badge:
            badge = len(kwargs["table_data"])
        else:
            badge = None
        super(TabularResultTab, self).__init__(
            tabid, title, template, show_badge=show_badge, badge=badge, **kwargs)


class EntryIdParser():

    og_re = re.compile("group_([0-9]*)")
    cog_re = re.compile("COG([0-9]{4})")
    pfam_re = re.compile("PF([0-9]{4,5})")
    ko_re = re.compile("K([0-9]{5})")
    vf_re = re.compile("VFG[0-9]{6}")

    def __init__(self, db):
        self.db = db

    def id_to_object_type(self, identifier):
        match = self.og_re.match(identifier)
        parsed_id = match and int(match.groups()[0])
        if parsed_id and self.db.check_orthogroup_entry_id(parsed_id):
            return "orthogroup", parsed_id

        match = self.cog_re.match(identifier)
        parsed_id = match and int(match.groups()[0])
        if parsed_id and self.db.check_cog_entry_id(parsed_id):
            return "cog", parsed_id

        match = self.pfam_re.match(identifier)
        parsed_id = match and int(match.groups()[0])
        if parsed_id and self.db.check_pfam_entry_id(parsed_id):
            return "pfam", parsed_id

        match = self.ko_re.match(identifier)
        parsed_id = match and int(match.groups()[0])
        if parsed_id and self.db.check_ko_entry_id(parsed_id):
            return "ko", parsed_id

        match = self.vf_re.match(identifier)
        if match and self.db.check_vf_entry_id(identifier):
            return "vf", identifier

        if self.db.check_amr_entry_id(identifier):
            return "amr", identifier

        return None


def locusx_genomic_region(db, seqid, window):
    hsh_loc = db.get_gene_loc([seqid])
    strand, start, end = hsh_loc[seqid]
    window_start, window_stop = start - window, start + window

    bioentry, _, contig_size, _ = db.get_bioentry_list(seqid, search_on="seqid")
    qualifiers = db.get_bioentry_qualifiers(bioentry)
    is_circular = "circular" in qualifiers["value"].values
    df_seqids = db.get_features_location(bioentry, search_on="bioentry_id")

    if 2*window >= contig_size:
        window_start = 0
        window_stop = contig_size
    elif window_start < 0 and not is_circular:
        window_start = 0
        df_seqids = df_seqids[df_seqids.start_pos < window_stop]
    elif window_stop > contig_size and not is_circular:
        window_stop = contig_size
        df_seqids = df_seqids[df_seqids.end_pos > window_start]
    elif window_start < 0:
        # circular contig
        diff = contig_size+window_start
        mask_circled = (df_seqids.end_pos >= diff)
        mask_same = (df_seqids.start_pos <= window_stop)

        df_seqids.loc[mask_same, "start_pos"] -= window_start
        df_seqids.loc[mask_same, "end_pos"] -= window_start
        df_seqids.loc[mask_circled, "start_pos"] -= diff
        df_seqids.loc[mask_circled, "end_pos"] -= diff
        df_seqids = df_seqids.loc[mask_same | mask_circled]
        window_stop -= window_start
        window_start = 0
        min_val = df_seqids["start_pos"].min()
        if min_val < 0:
            df_seqids["start_pos"] -= min_val
            df_seqids["end_pos"] -= min_val
            window_start -= min_val
            window_stop -= min_val

    elif window_stop > contig_size:
        # circular contig
        diff = window_stop-contig_size

        mask_same = (df_seqids.end_pos >= window_start)
        mask_circled = (df_seqids.start_pos <= diff)

        df_seqids.loc[mask_same, "start_pos"] -= diff
        df_seqids.loc[mask_same, "end_pos"] -= diff
        df_seqids.loc[mask_circled, "start_pos"] += (contig_size-diff)
        df_seqids.loc[mask_circled, "end_pos"] += (contig_size-diff)
        df_seqids = df_seqids.loc[mask_same | mask_circled]
        window_start -= diff
        window_stop = contig_size
    else:
        df_seqids = df_seqids[(df_seqids.end_pos>window_start)]
        df_seqids = df_seqids[(df_seqids.start_pos < window_stop)]

    if len(df_seqids) != len(df_seqids["seqfeature_id"].unique()):
        # This case may happen when a gene overlaps the break of a circular contig.
        # The location of this gene will be coded as join(...,...) in the gbk file
        # and stored as two separate genes with the same seqid in BioSQL.
        # If we want to display the whole contig as a continuous sequence, it is necessary
        # to detect this and manually merge this overlapping gene.
        grouped = df_seqids[["seqfeature_id", "strand", "end_pos",
                             "start_pos"]].groupby("seqfeature_id")
        start = grouped["start_pos"].min()
        end = grouped["end_pos"].max()
        strands = df_seqids[["seqfeature_id", "strand"]].drop_duplicates(
            "seqfeature_id")
        df_seqids = start.to_frame().join(end).join(
            strands.set_index("seqfeature_id"))
    else:
        df_seqids = df_seqids.set_index("seqfeature_id")

    # Some parts are redundant with get_features_location
    # those two function should be merged at some point
    infos = db.get_proteins_info(df_seqids.index.tolist(),
                                 to_return=["gene", "locus_tag", "product"],
                                 as_df=True, inc_non_CDS=True, inc_pseudo=True)
    cds_type = db.get_CDS_type(df_seqids.index.tolist())
    all_infos = infos.join(cds_type)
    all_infos = all_infos.join(df_seqids)
    return all_infos, window_start, window_stop


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


def make_div(figure_or_data, include_plotlyjs=False, show_link=False,
             div_id=None):
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


class AccessionFieldHandler():

    plasmid_prefix = "plasmid:"
    group_prefix = "group:"
    _db = None

    @classmethod
    def is_plasmid(cls, key):
        return key.startswith(cls.plasmid_prefix)

    def plasmid_key_to_id(self, key):
        return int(key.lstrip(self.plasmid_prefix))

    def plasmid_id_to_key(self, identifier):
        return f"{self.plasmid_prefix}{identifier}"

    @classmethod
    def is_group(cls, key):
        return key.startswith(cls.group_prefix)

    def group_key_to_id(self, key):
        return key.rsplit(self.group_prefix, 1)[-1]

    def group_id_to_key(self, identifier):
        return f"{self.group_prefix}{identifier}"

    @classmethod
    def is_taxid(cls, key):
        return not (cls.is_group(key) or cls.is_plasmid(key))

    @property
    def db(self):
        if self._db is None:
            biodb_path = settings.BIODB_DB_PATH
            self._db = DB.load_db_from_name(biodb_path)
        return self._db

    def get_choices(self, with_plasmids=True, exclude=[], exclude_taxids_in_groups=[]):
        result = self.db.get_genomes_description()
        result.set_index(result.index.astype(str), inplace=True)
        accession_choices = []
        for taxid, data in result.iterrows():
            accession_choices.append((taxid, data.description))
            if with_plasmids and data.has_plasmid:
                # Distinguish plasmids from taxons
                plasmid = self.plasmid_id_to_key(taxid)
                accession_choices.append((plasmid,
                                          f"{data.description} plasmid"))

        accession_choices.extend([(self.group_id_to_key(group[0]), group[0])
                                  for group in self.db.get_groups()])

        # We also exclude taxids contained in the excluded groups
        groups_to_exclude = [self.group_key_to_id(key) for key in exclude
                             if self.is_group(key)]
        in_groups = self.db.get_taxids_for_groups(groups_to_exclude)
        exclude = set(exclude).union({str(el) for el in in_groups})

        # And we exclude groups containing an excluded taxid
        taxids_to_exclude = list(filter(self.is_taxid, exclude))
        exclude = exclude.union(
            {self.group_id_to_key(groupid) for groupid in
             self.db.get_groups_containing_taxids(taxids_to_exclude)})

        # Finally we exclude taxids from groups in exclude_taxids_in_groups
        groups_to_exclude = [self.group_key_to_id(key)
                             for key in exclude_taxids_in_groups
                             if self.is_group(key)]
        in_groups = self.db.get_taxids_for_groups(groups_to_exclude)
        exclude = set(exclude).union({str(el) for el in in_groups})

        accession_choices = filter(lambda choice: choice[0] not in exclude,
                                   accession_choices)
        return tuple(accession_choices)

    def extract_choices(self, indices, include_plasmids):
        plasmids = []
        groups = []
        taxids = set()
        for key in indices:
            if self.is_plasmid(key):
                plasmids.append(self.plasmid_key_to_id(key))
            elif self.is_group(key):
                groups.append(self.group_key_to_id(key))
            else:
                taxids.add(int(key))

        if not include_plasmids:
            plasmids = None

        taxids.update(self.db.get_taxids_for_groups(groups))
        return list(taxids), plasmids
