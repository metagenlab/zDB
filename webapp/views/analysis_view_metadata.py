class AnalysisViewMetadata:
    _types_available_for = ["amr", "cog", "ko", "pfam", "orthogroup", "vf"]
    static_icon = True

    def __init__(self, object_type, object_name_plural):
        self.object_type = object_type
        self.object_name_plural = object_name_plural

    @classmethod
    def available_for(cls, object_type):
        return object_type in cls._types_available_for

    @property
    def description(self):
        return self._description.format(self.object_name_plural)

    @property
    def url(self):
        return self._url.format(self.object_type)


class EntryListMetadata(AnalysisViewMetadata):
    name = "index"
    _types_available_for = ["amr", "cog", "ko", "pfam", "vf"]
    title = "Index"
    _description = "Index of all {} identified in all genomes"
    _url = "/entry_list_{}"
    icon = "/img/icons8-index-50.png"


class ExtractionMetadata(AnalysisViewMetadata):
    name = "extract"
    title = "Extraction Form"
    _description = "List of {} present/absent from selected genomes"
    _url = "/extract_{}/"
    icon = "/img/icons8-compare-64.png"


class VennMetadata(AnalysisViewMetadata):
    name = "venn"
    title = "Venn Diagram"
    _description = "Venn diagrams of unique and shared {} among selected genomes"
    _url = "/venn_{}/"
    icon = "https://img.icons8.com/plasticine/100/000000/venn-diagram.png"
    static_icon = False


class TabularComparisonMetadata(AnalysisViewMetadata):
    name = "table_comp"
    title = "Presence/absence table"
    _description = "Presence/absence table of {} in selected genomes"
    _url = "/{}_comparison/"
    icon = "/img/icons8-spreadsheet-80.png"


class HeatmapMetadata(AnalysisViewMetadata):
    name = "heatmap"
    title = "Heatmap"
    _description = "Heatmap of {} presence/absence in selected genomes"
    _url = "/plot_heatmap/{}"
    icon = "https://img.icons8.com/plasticine/100/000000/heat-map.png"
    static_icon = False


class AccumulationRarefactionMetadata(AnalysisViewMetadata):
    name = "pan_genome"
    title = "Accumulation / rarefaction plot"
    _description = "Accumulation / rarefaction plot of {} in selected genomes"
    _url = "/pan_genome/{}/"
    icon = "https://img.icons8.com/plasticine/100/000000/line-chart.png"
    static_icon = False


class CategoriesBarchartMetadata(AnalysisViewMetadata):
    name = "bar"
    _types_available_for = ["cog", "ko"]
    title = "Categories barchart"
    _description = "Barcharts of {} categories in selected genomes"
    _url = "/{}_barchart/"
    icon = "/img/barchart-1.1s-200px.svg"


class CategoriesFreqHeatmapMetadata(AnalysisViewMetadata):
    name = "cat_freq"
    _types_available_for = ["cog"]
    title = "Categories freq. heatmap"
    _description = "Heatmap of frequencies of genes identifed in each COG category"
    _url = "/{}_phylo_heatmap/True/"
    icon = "/img/icons8-heat-map-100.phylo.svg"


class CategoriesCountHeatmapMetadata(AnalysisViewMetadata):
    name = "cat_count"
    _types_available_for = ["cog"]
    title = "Categories count heatmap"
    _description = "Heatmap of counts of genes identifed in each COG category"
    _url = "/{}_phylo_heatmap/False/"
    icon = "/img/icons8-heat-map-100.phylo.svg"


class GwasMetadata(AnalysisViewMetadata):
    name = "gwas"
    _types_available_for = ["amr", "cog", "ko", "pfam", "vf", "orthogroup"]
    title = "Genome Wide Association Study"
    _description = "Association of {} with user-defined phenotype"
    _url = "/gwas_{}/"
    icon = "/img/gwas.png"


analysis_views_metadata = [
    EntryListMetadata,
    ExtractionMetadata,
    VennMetadata,
    TabularComparisonMetadata,
    HeatmapMetadata,
    AccumulationRarefactionMetadata,
    CategoriesBarchartMetadata,
    CategoriesFreqHeatmapMetadata,
    CategoriesCountHeatmapMetadata,
    GwasMetadata,
]
