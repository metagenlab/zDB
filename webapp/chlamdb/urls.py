from django.urls import re_path
from django.views.generic import TemplateView
from django.views.generic.base import RedirectView
from views import autocomplete
from views import custom_plots
from views import entry_lists
from views import fam
from views import genomic_island
from views import groups
from views import gwas
from views import hits_extraction
from views import locus
from views import tabular_comparison
from views import venn
from views import views

favicon_view = RedirectView.as_view(url="/assets/favicon.ico", permanent=True)


urlpatterns = [
    re_path(
        r"^vf_comparison",
        tabular_comparison.VfComparisonView.as_view(),
        name="vf_comparison",
    ),  # noqa
    re_path(r"^venn_vf/$", venn.VennVfView.as_view(), name="venn_vf"),
    re_path(r"^venn_pfam/$", venn.VennPfamView.as_view(), name="venn_pfam"),
    re_path(
        r"^venn_orthogroup/$", venn.VennOrthogroupView.as_view(), name="venn_orthogroup"
    ),  # noqa
    re_path(r"^venn_ko/$", venn.VennKoView.as_view(), name="venn_ko"),
    re_path(r"^venn_gic/$", venn.VennGiView.as_view(), name="venn_gic"),
    re_path(r"^venn_cog/$", venn.VennCogView.as_view(), name="venn_cog"),
    re_path(r"^venn_amr/$", venn.VennAmrView.as_view(), name="venn_amr"),
    re_path(r"^search_suggest/.*$", views.search_suggest, name="search_suggest"),  # noqa
    re_path(
        r"^search_bar/([a-zA-Z0-9_\.\-]+)/([a-zA-Z0-9_\.\-]+)",
        views.search_bar,
        name="search_bar",
    ),  # noqa
    re_path(r"^search_bar$", views.search_bar, name="search_bar"),
    re_path(r"^plot_region/$", views.plot_region, name="plot_region"),
    re_path(
        r"^plot_heatmap/([a-zA-Z0-9_\-]+)",
        views.PlotHeatmap.as_view(),
        name="plot_heatmap",
    ),  # noqa
    re_path(r"^phylogeny", views.phylogeny, name="phylogeny"),
    re_path(
        r"^pfam_comparison",
        tabular_comparison.PfamComparisonView.as_view(),
        name="pfam_comparison",
    ),  # noqa
    re_path(
        r"^pan_genome/([a-zA-Z0-9_]+)", views.PanGenome.as_view(), name="pan_genome"
    ),  # noqa
    re_path(
        r"^orthogroup_comparison",
        tabular_comparison.OrthogroupComparisonView.as_view(),
        name="orthogroup_comparison",
    ),  # noqa
    re_path(
        r"^orthogroup/([a-zA-Z0-9_\.\-]+)",
        locus.Orthogroup.as_view(),
        name="orthogroup",
    ),  # noqa
    re_path(r"^module_comparison/$", views.module_comparison, name="module_comparison"),  # noqa
    re_path(
        r"^module_cat_info/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.\+-]+)$",
        views.module_cat_info,
        name="module_cat_info",
    ),  # noqa
    re_path(
        r"^locusx/([a-zA-Z0-9_\.\-]+)/([a-zA-Z0-9_\.\-]+)",
        locus.LocusX.as_view(),
        name="locusx",
    ),  # noqa
    re_path(r"^locusx/([a-zA-Z0-9_\.\-]+)", locus.LocusX.as_view(), name="locusx"),
    re_path(r"^locusx$", locus.LocusX.as_view(), name="locusx"),
    re_path(
        r"^ko_venn_subset/([a-zA-Z0-9_\.\+-]+)$",
        venn.VennKoSubsetView.as_view(),
        name="ko_venn_subset",
    ),  # noqa
    re_path(
        r"^ko_comparison",
        tabular_comparison.KoComparisonView.as_view(),
        name="ko_comparison",
    ),  # noqa
    re_path(r"^ko_barchart/$", views.KoBarchart.as_view(), name="ko_barchart"),
    re_path(
        r"^kegg_module_subcat$", views.kegg_module_subcat, name="kegg_module_subcat"
    ),  # noqa
    re_path(
        r"^KEGG_module_map/([a-zA-Z0-9_\.]+)$",
        views.KEGG_module_map,
        name="KEGG_module_map",
    ),  # noqa
    re_path(r"^kegg_module/$", views.kegg_module, name="kegg_module"),
    re_path(
        r"^KEGG_mapp_ko/([a-zA-Z0-9_\.]+)/([0-9]+)$",
        views.KEGG_mapp_ko,
        name="KEGG_mapp_ko",
    ),  # noqa
    re_path(
        r"^KEGG_mapp_ko/([a-zA-Z0-9_\.]+)$", views.KEGG_mapp_ko, name="KEGG_mapp_ko"
    ),  # noqa
    re_path(r"^KEGG_mapp_ko$", views.KEGG_mapp_ko, name="KEGG_mapp_ko"),
    re_path(
        r"^kegg_genomes_modules/$",
        views.kegg_genomes_modules,
        name="kegg_genomes_modules",
    ),  # noqa
    re_path(r"^kegg_genomes/$", views.kegg_genomes, name="kegg_genomes"),
    re_path(r"^kegg/$", views.kegg, name="kegg"),
    re_path(
        r"^index_comp/([a-zA-Z0-9_\.\-]+)",
        views.ComparisonIndexView.as_view(),
        name="index_comp",
    ),  # noqa
    re_path(r"^home/$", views.home, name="home"),
    re_path(r"^help", views.help, name="help"),
    re_path(r"^gwas_vf/", gwas.VfGwasView.as_view(), name="gwas_vf"),
    re_path(r"^gwas_pfam/", gwas.PfamGwasView.as_view(), name="gwas_pfam"),
    re_path(
        r"^gwas_orthogroup/", gwas.OrthogroupGwasView.as_view(), name="gwas_orthogroup"
    ),  # noqa
    re_path(r"^gwas_ko/", gwas.KoGwasView.as_view(), name="gwas_ko"),
    re_path(r"^gwas_cog/", gwas.CogGwasView.as_view(), name="gwas_cog"),
    re_path(r"^gwas_amr/", gwas.AmrGwasView.as_view(), name="gwas_amr"),
    re_path(
        r"^groups/add/$", groups.GroupAdd.as_view(), name=groups.GroupAdd.view_name
    ),  # noqa
    re_path(
        r"^groups/([a-zA-Z0-9_\.\(\)\-\'\s]+)/delete/$",
        groups.GroupDelete.as_view(),
        name=groups.GroupDelete.view_name,
    ),  # noqa
    re_path(
        r"^groups/([a-zA-Z0-9_\.\(\)\-\'\s]+)",
        groups.GroupDetails.as_view(),
        name=groups.GroupDetails.view_name,
    ),  # noqa
    re_path(
        r"^groups/$",
        groups.GroupsOverview.as_view(),
        name=groups.GroupsOverview.view_name,
    ),  # noqa
    re_path(
        r"^gic_comparison",
        tabular_comparison.GiComparisonView.as_view(),
        name="gic_comparison",
    ),  # noqa
    re_path(
        r"^get_cog/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.\%]+)$",
        views.get_cog,
        name="get_cog",
    ),  # noqa
    re_path(
        r"^genomic_island/([a-zA-Z0-9_\.]+)",
        genomic_island.GenomicIsland.as_view(),
        name=genomic_island.GenomicIsland.view_name,
    ),
    re_path(r"^genomes", views.Genomes.as_view(), name=views.Genomes.view_name),
    re_path(r"^favicon\.ico$", favicon_view),
    re_path(r"^FAQ", views.faq, name="FAQ"),
    re_path(r"^fam_vf/(VFG[0-9]+)$", fam.FamVfView.as_view(), name="fam_vf"),
    re_path(r"^fam_pfam/(PF[0-9]+)$", fam.FamPfamView.as_view(), name="fam_pfam"),
    re_path(r"^fam_ko/(K[0-9]+)$", fam.FamKoView.as_view(), name="fam_ko"),
    re_path(r"^fam_gic/(GIC[0-9]+)$", fam.FamGiClusterView.as_view(), name="fam_gic"),
    re_path(r"^fam_cog/(COG[0-9]+)$", fam.FamCogView.as_view(), name="fam_cog"),  # noqa
    re_path(
        r"^fam_amr/([a-zA-Z0-9_\.\(\)\-\']+)$", fam.FamAmrView.as_view(), name="fam_amr"
    ),  # noqa
    re_path(
        r"^extract_vf/$", hits_extraction.ExtractVfView.as_view(), name="extract_vf"
    ),  # noqa
    re_path(
        r"^extract_pfam/$",
        hits_extraction.ExtractPfamView.as_view(),
        name="extract_pfam",
    ),  # noqa
    re_path(
        r"^extract_orthogroup/$",
        hits_extraction.ExtractOrthogroupView.as_view(),
        name="extract_orthogroup",
    ),  # noqa
    re_path(
        r"^extract_ko/$", hits_extraction.ExtractKoView.as_view(), name="extract_ko"
    ),  # noqa
    re_path(
        r"^extract_gic/$", hits_extraction.ExtractGiView.as_view(), name="extract_gic"
    ),  # noqa
    re_path(
        r"^extract_contigs/([0-9]+)",
        views.ExtractContigs.as_view(),
        name="extract_contigs",
    ),  # noqa
    re_path(
        r"^extract_cog/$", hits_extraction.ExtractCogView.as_view(), name="extract_cog"
    ),  # noqa
    re_path(
        r"^extract_amr/$", hits_extraction.ExtractAmrView.as_view(), name="extract_amr"
    ),  # noqa
    re_path(
        r"^entry_list_vf$", entry_lists.VfEntryListView.as_view(), name="entry_list_vf"
    ),  # noqa
    re_path(
        r"^entry_list_pfam$",
        entry_lists.PfamEntryListView.as_view(),
        name="entry_list_pfam",
    ),  # noqa
    re_path(
        r"^entry_list_ko$", entry_lists.KoEntryListView.as_view(), name="entry_list_ko"
    ),  # noqa
    re_path(
        r"^entry_list_gic$",
        entry_lists.GiEntryListView.as_view(),
        name="entry_list_gic",
    ),  # noqa
    re_path(
        r"^entry_list_cog$",
        entry_lists.CogEntryListView.as_view(),
        name="entry_list_cog",
    ),  # noqa
    re_path(
        r"^entry_list_amr$",
        entry_lists.AmrEntryListView.as_view(),
        name="entry_list_amr",
    ),  # noqa
    re_path(
        r"^download_sequences$",
        locus.DownloadSequences.as_view(),
        name="download_sequences",
    ),  # noqa
    re_path(
        r"^custom_plots/$", custom_plots.CusomPlotsView.as_view(), name="custom_plots"
    ),  # noqa
    re_path(
        r"^cog_venn_subset/([A-Z])$",
        venn.VennCogSubsetView.as_view(),
        name="cog_venn_subset",
    ),  # noqa
    re_path(
        r"^cog_phylo_heatmap/([a-zA-Z0-9_\-]+)",
        views.CogPhyloHeatmap.as_view(),
        name="cog_phylo_heatmap",
    ),  # noqa
    re_path(
        r"^cog_comparison",
        tabular_comparison.CogComparisonView.as_view(),
        name="cog_comparison",
    ),  # noqa
    re_path(r"^cog_barchart/$", views.CogBarchart.as_view(), name="cog_barchart"),
    re_path(r"^circos_main/$", views.circos_main, name="circos_main"),
    re_path(r"^circos/$", views.circos, name="circos"),
    re_path(r"^blast/$", views.blast, name="blast"),
    re_path(
        r"^autocomplete_taxid/$",
        autocomplete.AutocompleteTaxid.as_view(),
        name="autocomplete_taxid",
    ),  # noqa
    re_path(
        r"^autocomplete_n_missing/$",
        autocomplete.AutocompleteNMissing.as_view(),
        name="autocomplete_n_missing",
    ),  # noqa
    re_path(
        r"^amr_comparison",
        tabular_comparison.AmrComparisonView.as_view(),
        name="amr_comparison",
    ),  # noqa
    re_path(r"^about$", views.about, name="about"),
    re_path(r"^.*$", views.home, name="home"),
    re_path(
        "^robots.txt$",
        TemplateView.as_view(template_name="robots.txt", content_type="text/plain"),
    ),  # noqa
    # re_path(r'^FAQ',TemplateView.as_view(template_name='FAQ.html')),
]
