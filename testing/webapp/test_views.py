from chlamdb.urls import urlpatterns
from django.test import SimpleTestCase
from django.urls import resolve

broken_views = [
    '/orthogroup_list_cog_barchart/',
    '/orthogroup_list_cog_barchart/True/',
    '/circos_main/',
    '/genomes_intro/',
]

untested_patterns = {
    '^robots.txt$',
    '^search_suggest/.*$',
    '^search_bar/([a-zA-Z0-9_\\.\\-]+)/([a-zA-Z0-9_\\.\\-]+)',
    '^module_cat_info/([a-zA-Z0-9_\\.]+)/([a-zA-Z0-9_\\.\\+-]+)$',
    '^favicon\\.ico$',
    '^search_bar$',
    '^plot_region/$',
}

urls = [
    '/home/',
    '/cog_barchart/',
    '/pan_genome/COG',
    '/pan_genome/Pfam',
    '/pan_genome/ko',
    '/pan_genome/orthology',
    '/plot_heatmap/COG',
    '/plot_heatmap/Pfam',
    '/plot_heatmap/ko',
    '/plot_heatmap/orthology',
    '/COG_phylo_heatmap/True',
    '/COG_phylo_heatmap/False',
    '/module_barchart/',
    '/get_cog/3/L?h=1&h=2&h=3',
    '/ko_venn_subset/L?h=1&h=2&h=3',
    '/kegg/',
    '/locusx/CHUV_00025/True',
    '/index_comp/orthology',
    '/index_comp/COG',
    '/index_comp/Pfam',
    '/index_comp/ko',
    '/locusx/CHUV_00025',
    '/orthogroup/group_85',
    '/locusx',
    '/entry_list_pfam',
    '/entry_list_cog',
    '/entry_list_ko',
    '/circos/',
    '/blast/',
    '/extract_orthogroup/',
    '/venn_orthogroup/',
    '/extract_cog/',
    '/venn_cog/',
    '/cog_venn_subset/L?h=1&h=2',
    '/venn_cog/True',
    '/venn_ko/',
    '/extract_pfam/',
    '/extract_pfam/taxon_id',
    '/extract_ko/',
    '/venn_pfam/',
    '/KEGG_mapp_ko',
    '/KEGG_mapp_ko/map01100',
    '/KEGG_mapp_ko/map00230/1',
    '/KEGG_module_map/M00035',
    '/kegg_module_subcat',
    '/kegg_module/',
    '/module_comparison/',
    '/pfam_comparison',
    '/cog_comparison',
    '/ko_comparison',
    '/orthogroup_comparison',
    '/amr_comparison',
    '/about',
    '/help',
    '/fam_pfam/PF10423',
    '/fam_cog/COG0775',
    '/fam_ko/K01241',
    '/FAQ',
    '/phylogeny',
    '/extract_contigs/1',
    '/genomes',
    '/kegg_genomes/',
    '/kegg_genomes_modules/',
    '/invalid/url'
]


class TestViewsAreHealthy(SimpleTestCase):

    def test_no_broken_views(self):
        for url in urls:
            resp = self.client.get(url)
            self.assertEqual(200, resp.status_code, f"{url} is broken")

    def test_broken_views(self):
        for url in broken_views:
            self.assertRaises(Exception, self.client.get, url)

    def test_all_urlpatterns_are_tested(self):
        all_patterns = set(el.pattern._regex for el in urlpatterns)
        tested_patterns = set()
        for url in urls:
            tested_patterns.add(resolve(url.rsplit("?")[0]).route)

        for url in broken_views:
            tested_patterns.add(resolve(url.rsplit("?")[0]).route)

        covered_patterns = tested_patterns | untested_patterns
        self.assertFalse(
            all_patterns - covered_patterns,
            "Some patterns are not covered in the tests: please add them to "
            "one of broken_views, untested_patterns or urls")

    def test_all_views_render_valid_html(self):
        for url in urls:
            resp = self.client.get(url)
            try:
                self.assertContains(resp, "", html=True)
            except Exception as exc:
                print(f"\n\nInvalid html for {url}")
                print(f"Templates {[el.name for el in resp.templates]}\n\n")
                raise exc


class TestViewsContent(SimpleTestCase):

    def assertNoPlot(self, resp):
        self.assertNotIn("envoi", resp.context.keys())
        self.assertNotContains(resp, '<div class="panel panel-success"')
        self.assertNotContains(resp, "Help to interpret the results")

    def assertPlot(self, resp):
        self.assertIn("envoi", resp.context.keys())
        self.assertContains(resp, '<div class="panel panel-success"')
        self.assertContains(resp, "Help to interpret the results")

    def assertTitle(self, resp, title):
        self.assertContains(
            resp, f'<p class="home-title">{title}</p>', html=True)

    def test_home_page(self):
        resp = self.client.get("/home")
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/home.html')
        self.assertEqual(3, resp.context["number_of_files"])
        self.assertContains(resp, "A comparative genomics database", html=True)
        self.assertContains(resp, '<a href="/genomes" ><b>Genomes</b></a>', html=True)
        self.assertContains(resp, '<a href="/phylogeny"><b>Phylogeny</b></a>', html=True)
        self.assertContains(resp, '<a class="link_boxes" href="/blast/"><span class="link"></span> Blast </a>', html=True)
        self.assertContains(resp, '<a class="link_boxes" href="/index_comp/orthology"><span class="link"></span>Orthogroups</a>', html=True)
        self.assertContains(resp, '<a class="link_boxes" href="/index_comp/COG"><span class="link"></span>COGs</a>', html=True)
        self.assertContains(resp, '<a class="link_boxes" href="/index_comp/ko"><span class="link"></span>Kegg Orthologs</a>', html=True)
        self.assertContains(resp, '<a class="link_boxes" href="/index_comp/Pfam"><span class="link"></span>Pfam domains</a>', html=True)
        self.assertContains(resp, '<a class="link_boxes" href="/circos/"><span class="link"></span>Plot region</a>', html=True)
        self.assertContains(resp, '<a class="link_boxes" href="/circos/"><span class="link"></span>Circos</a>', html=True)
        self.assertContains(resp, '<a class="link_boxes" href="/kegg/"><span class="link"></span>Kegg based</a>', html=True)

    def test_cog_barchart(self):
        resp = self.client.get("/cog_barchart/")
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/cog_barplot.html')
        self.assertTitle(resp, "Comparisons: Clusters of Orthologous groups (COGs)")
        self.assertNoPlot(resp)
        self.assertContains(resp, "Distribution of COGs within COG categories")

        resp = self.client.post("/cog_barchart/", data={"targets": ["0", "1"]})
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/cog_barplot.html')
        self.assertTitle(resp, "Comparisons: Clusters of Orthologous groups (COGs)")
        self.assertPlot(resp)
        self.assertContains(resp, "Distribution of COGs within COG categories")

    def test_pan_genome_cog(self):
        resp = self.client.get("/pan_genome/COG")
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/pan_genome.html')
        self.assertEqual("COG", resp.context["type"])
        self.assertTitle(resp, "Comparisons: Clusters of Orthologous groups (COGs)")
        self.assertNoPlot(resp)

        resp = self.client.post("/pan_genome/COG", data={"targets": ["0", "1"]})
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/pan_genome.html')
        self.assertEqual("COG", resp.context["type"])
        self.assertTitle(resp, "Comparisons: Clusters of Orthologous groups (COGs)")
        self.assertPlot(resp)
