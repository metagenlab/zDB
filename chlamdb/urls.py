from django.conf.urls import url
from . import views
from django.contrib.auth import logout
from django.contrib.sitemaps.views import sitemap
from django.views.generic import TemplateView
from django.urls import include, path
from django.conf import settings
from django.urls import reverse
from django.contrib.sitemaps import Sitemap
from django.views import static
from django.views.generic.base import RedirectView

favicon_view = RedirectView.as_view(url='/assets/favicon.ico', permanent=True)

class ViewSitemap(Sitemap):
    """Reverse 'static' views for XML sitemap."""

    def items(self):
        # Return list of url names for views to include in sitemap
        return ['home', 'about']

    def location(self, item):
        #site = Site(domain='chlamdb.ch', name='chlamdb.ch')
        return reverse(item)


sitemaps = {'views': ViewSitemap}

# url(r'^sitemap.xml$', sitemap, {'sitemaps': sitemaps}),

urlpatterns = [        
    url('^robots.txt$', TemplateView.as_view(template_name='robots.txt', content_type='text/plain')),
    url(r'^sitemap$', views.sitemap, name="sitemap"),
    url(r'^home/$', views.home, name="home"),
    url(r'^cog_barchart/$', views.cog_barchart, name="cog_barchart"),
    url(r'^pan_genome/([a-zA-Z0-9_]+)', views.pan_genome, name="pan_genome"),
    url(r'^COG_phylo_heatmap/([a-zA-Z0-9_\-]+)', views.COG_phylo_heatmap, name="COG_phylo_heatmap"),
    url(r'^module_barchart/$', views.module_barchart, name="module_barchart"),
    url(r'^plot_heatmap/([a-zA-Z0-9_\-]+)', views.plot_heatmap, name="plot_heatmap"),
    url(r'^get_cog/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.\%]+)$', views.get_cog, name="get_cog"),
    url(r'^module_cat_info/([a-zA-Z0-9_\.]+)/([a-zA-Z0-9_\.\+-]+)$', views.module_cat_info, name="module_cat_info"),
    url(r'^ko_venn_subset/([a-zA-Z0-9_\.\+-]+)$', views.ko_venn_subset, name="ko_venn_subset"),
    url(r'^hydropathy/([a-zA-Z0-9_\.\+-]+)$', views.hydropathy, name="hydropathy"),
    url(r'^hmm2circos/$', views.hmm2circos, name="hmm2circos"),
    url(r'^priam_kegg/$', views.priam_kegg, name="priam_kegg"),
    url(r'^locusx/([a-zA-Z0-9_\.\-]+)/([a-zA-Z0-9_\.\-]+)', views.locusx, name="locusx"),
    url(r'^search_bar$', views.search_bar, name="search_bar"),
    url(r'^search_bar/([a-zA-Z0-9_\.\-]+)/([a-zA-Z0-9_\.\-]+)', views.search_bar, name="search_bar"),
    url(r'^locusx/([a-zA-Z0-9_\.\-]+)', views.locusx, name="locusx"),
    url(r'^orthogroup/([a-zA-Z0-9_\.\-]+)', views.orthogroup, name="orthogroup"),
    url(r'^locusx$', views.locusx, name="locusx"),
    url(r'^blastswissprot/([a-zA-Z0-9_\.\-]+)$', views.blastswissprot, name="blastswissprot"),
    url(r'^plot_region/$', views.plot_region, name="plot_region"),
    url(r'^circos/$', views.circos, name="circos"),
    url(r'^circos_main/$', views.circos_main, name="circos_main"),
    url(r'^orthogroup_list_cog_barchart/$', views.orthogroup_list_cog_barchart, name="orthogroup_list_cog_barchart"),
    url(r'^orthogroup_list_cog_barchart/([a-zA-Z0-9_]+)/$', views.orthogroup_list_cog_barchart, name="orthogroup_list_cog_barchart"),
    url(r'^blast/$', views.blast, name="blast"),
    url(r'^extract_orthogroup/$', views.extract_orthogroup, name="extract_orthogroup"),
    url(r'^extract_orthogroup/([a-zA-Z0-9_]+)$$', views.extract_orthogroup, name="extract_orthogroup"),
    url(r'^venn_orthogroup/$', views.venn_orthogroup, name="venn_orthogroup"),
    url(r'^extract_cog/([a-zA-Z0-9_]+)$$', views.extract_cog, name="extract_cog"),
    url(r'^extract_cog/$', views.extract_cog, name="extract_cog"),
    url(r'^venn_cog/$', views.venn_cog, name="venn_cog"),
    url(r'^cog_venn_subset/([A-Z])$', views.cog_venn_subset, name="cog_venn_subset"),
    url(r'^venn_cog/([a-zA-Z0-9_]+)$$', views.venn_cog, name="venn_cog"),
    url(r'^venn_ko/$', views.venn_ko, name="venn_ko"),
    url(r'^extract_pfam/$', views.extract_pfam, name="extract_pfam"),
    url(r'^extract_pfam/([a-zA-Z0-9_]+)$', views.extract_pfam, name="extract_pfam"),
    url(r'^extract_ko/$', views.extract_ko, name="extract_ko"),
    url(r'^extract_ko/([a-zA-Z0-9_]+)$', views.extract_ko, name="extract_ko"),
    url(r'^venn_pfam/$', views.venn_pfam, name="venn_pfam"),
    url(r'^KEGG_mapp_ko/([a-zA-Z0-9_\.]+)$', views.KEGG_mapp_ko, name="KEGG_mapp_ko"),
    url(r'^KEGG_mapp_ko/([a-zA-Z0-9_\.]+)/([0-9]+)$', views.KEGG_mapp_ko, name="KEGG_mapp_ko"),
    url(r'^KEGG_module_map/([a-zA-Z0-9_\.]+)$', views.KEGG_module_map, name="KEGG_module_map"),
    url(r'^kegg_module_subcat$', views.kegg_module_subcat, name="kegg_module_subcat"),
    url(r'^kegg_module/$', views.kegg_module, name="kegg_module"),
    url(r'^module_comparison/$', views.module_comparison, name="module_comparison"),
    url(r'^pfam_comparison', views.pfam_comparison, name="pfam_comparison"),
    url(r'^ko_comparison', views.ko_comparison, name="ko_comparison"),
    url(r'^orthogroup_comparison', views.orthogroup_comparison, name="orthogroup_comparison"),
    url(r'^logout/$', logout, {'next_page': '/'},),
    url(r'^about$', views.about, name="about"),
    url(r'^help', views.help, name="help"),

    url(r'^fam_pfam/(PF[0-9]+)$', views.fam_pfam, name="fam_pfam"),
    url(r'^fam_cog/(COG[0-9]+)$', views.fam_cog, name="fam_cog"),
    url(r'^fam_ko/(K[0-9]+)$', views.fam_ko, name="fam_ko"),
    url(r'^get-task-info/', views.get_task_info),
    url(r'^docs/(?P<path>.*)$', static.serve, {'document_root': settings.DOCS_ROOT}),
    url(r'^favicon\.ico$', favicon_view),
    url(r'^FAQ', views.faq, name='FAQ'),
    url(r'^phylogeny_intro', views.phylogeny_intro, name='phylogeny_intro'),
    url(r'^genomes_intro', views.genomes_intro, name='genomes_intro'),
    url(r'^extract_contigs/([0-9]+)', views.extract_contigs, name='extract_contigs'),
    url(r'^extract_region', views.extract_region, name='extract_region'),
    url(r'^.*$', views.home, name="home"),
    #url(r'^FAQ',TemplateView.as_view(template_name='FAQ.html')),
]


