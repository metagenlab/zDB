import os
import re
from unittest import skip

from chlamdb.urls import urlpatterns
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import SimpleTestCase
from django.urls import resolve

dump_html = False
html_dumps_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "html_dumps")

untested_patterns = {
    '^favicon\\.ico$',
    '^module_cat_info/([a-zA-Z0-9_\\.]+)/([a-zA-Z0-9_\\.\\+-]+)$',
    '^plot_region/$',
    '^robots.txt$',
    '^search_bar$',
    '^search_bar/([a-zA-Z0-9_\\.\\-]+)/([a-zA-Z0-9_\\.\\-]+)',
    '^search_suggest/.*$',
}

urls = [
    '/about',
    '/amr_comparison',
    '/blast/',
    '/circos/',
    '/circos_main/',
    '/cog_barchart/',
    '/cog_comparison',
    '/cog_phylo_heatmap/False',
    '/cog_phylo_heatmap/True',
    '/cog_venn_subset/L?h=1&h=2',
    '/entry_list_amr',
    '/entry_list_cog',
    '/entry_list_ko',
    '/entry_list_pfam',
    '/entry_list_vf',
    '/extract_amr/',
    '/extract_cog/',
    '/extract_contigs/1',
    '/extract_ko/',
    '/extract_orthogroup/',
    '/extract_pfam/',
    '/extract_vf/',
    '/fam_amr/ybtP',
    '/fam_cog/COG0775',
    '/fam_ko/K01241',
    '/fam_pfam/PF10423',
    '/fam_vf/VFG029227',
    '/FAQ',
    '/genomes',
    '/get_cog/3/L?h=1&h=2&h=3',
    '/gwas_amr/',
    '/gwas_cog/',
    '/gwas_ko/',
    '/gwas_orthogroup/',
    '/gwas_pfam/',
    '/gwas_vf/',
    '/help',
    '/home/',
    '/index_comp/cog',
    '/index_comp/ko',
    '/index_comp/orthogroup',
    '/index_comp/pfam',
    '/index_comp/vf',
    '/invalid/url',
    '/kegg/',
    '/kegg_genomes/',
    '/kegg_genomes_modules/',
    '/KEGG_mapp_ko',
    '/KEGG_mapp_ko/map00230/1',
    '/KEGG_mapp_ko/map01100',
    '/kegg_module/',
    '/KEGG_module_map/M00035',
    '/kegg_module_subcat',
    '/ko_barchart/',
    '/ko_comparison',
    '/ko_venn_subset/Cofactor+and+vitamin+metabolism?h=1&h=2',
    '/locusx',
    '/locusx/CHUV_00025',
    '/locusx/CHUV_00025/True',
    '/module_comparison/',
    '/orthogroup/group_85',
    '/orthogroup_comparison',
    '/pan_genome/amr',
    '/pan_genome/cog',
    '/pan_genome/ko',
    '/pan_genome/orthogroup',
    '/pan_genome/pfam',
    '/pan_genome/vf',
    '/pfam_comparison',
    '/phylogeny',
    '/plot_heatmap/cog',
    '/plot_heatmap/ko',
    '/plot_heatmap/orthogroup',
    '/plot_heatmap/pfam',
    '/plot_heatmap/amr',
    '/plot_heatmap/vf',
    '/venn_amr/',
    '/venn_cog/',
    '/venn_ko/',
    '/venn_orthogroup/',
    '/venn_pfam/',
    '/venn_vf/',
    '/vf_comparison/',
]


def maybe_dump_html(resp, fname_postfix=""):
    if dump_html:
        if not os.path.isdir(html_dumps_dir):
            os.mkdir(html_dumps_dir)
        fname = resp.request["PATH_INFO"].strip("/").replace("/", "-")
        if fname_postfix:
            fname += "_" + fname_postfix
        fname += ".html"

        # Replace csrf token
        match = re.search(b'<input type="hidden" name="csrfmiddlewaretoken" value="(.*?)"', resp.content)
        if match:
            resp.content = resp.content.replace(match.groups()[0], b"XXXXX")

        with open(os.path.join(html_dumps_dir, fname), "wb") as dumpfile:
            dumpfile.write(resp.content)


class TestViewsAreHealthy(SimpleTestCase):

    def test_no_broken_views(self):
        for url in urls:
            resp = self.client.get(url)
            self.assertEqual(200, resp.status_code, f"{url} is broken")

    def test_all_urlpatterns_are_tested(self):
        all_patterns = set(el.pattern._regex for el in urlpatterns)
        tested_patterns = set()
        for url in urls:
            tested_patterns.add(resolve(url.rsplit("?")[0]).route)

        covered_patterns = tested_patterns | untested_patterns
        self.assertFalse(
            all_patterns - covered_patterns,
            "Some patterns are not covered in the tests: please add them to "
            "untested_patterns or urls")

    def test_all_views_render_valid_html(self):
        for url in urls:
            resp = self.client.get(url)
            try:
                self.assertContains(resp, "", html=True)
            except Exception as exc:
                print(f"\n\nInvalid html for {url}")
                print(f"Templates {[el.name for el in resp.templates]}\n\n")
                raise exc
            maybe_dump_html(resp)

    def test_all_comparison_views_render_navigation(self):

        comp_pattern = "/extract_[^s]+/|/[^m].+_comparison|/venn_|"\
                       "/plot_heatmap|/pan_genome|/entry_list"
        for url in urls:
            if re.match(comp_pattern, url):
                resp = self.client.get(url)
                self.assertContains(resp, "<nav>",
                                    msg_prefix=f"{url} does not contain navigation")


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
        self.assertContains(resp, '<a class="link_boxes" href="/index_comp/orthogroup"><span class="link"></span>Orthogroups</a>', html=True)
        self.assertContains(resp, '<a class="link_boxes" href="/index_comp/cog"><span class="link"></span>COGs</a>', html=True)
        self.assertContains(resp, '<a class="link_boxes" href="/index_comp/ko"><span class="link"></span>Kegg Orthologs</a>', html=True)
        self.assertContains(resp, '<a class="link_boxes" href="/index_comp/pfam"><span class="link"></span>Pfam domains</a>', html=True)
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

    def test_blast(self):
        resp = self.client.get("/blast/")
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/blast.html')
        self.assertTitle(resp, "Homology search: Blast")
        self.assertNotIn("envoi", resp.context.keys())
        self.assertNotContains(resp, "Help to interpret the results")
        self.assertNotContains(resp, 'id="phylo_distrib"')
        self.assertNotContains(resp, 'id="blast_details"')

        data = {
            "blast_input": "ATCGCCACGGTGGTGCAGGCGCAGAAAGCGGGCAAAACGCTCAGCGTCG",
            "blast": "blastn_ffn",
            "max_number_of_hits": 10,
            "target": "all"
        }
        resp = self.client.post("/blast/", data=data)
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/blast.html')
        self.assertTitle(resp, "Homology search: Blast")
        self.assertIn("envoi", resp.context.keys())
        self.assertContains(resp, "Help to interpret the results")
        self.assertContains(resp, 'id="phylo_distrib"')
        self.assertContains(resp, 'id="blast_details"')


class ComparisonViewsTestMixin():

    selection_html = """
        <option value="0" {0}>Klebsiella pneumoniae R6724_16313</option>
        <option value="1" {0}>Klebsiella pneumoniae R6726_16314</option>
        <option value="2">Klebsiella pneumoniae R6728_16315</option>
        """

    table_html = '<table class="hover" id="mytable"  style="padding-top: 1em;">'
    venn_html = '<div id="venn_diagram" '
    heatmap_html = '<div id="heatmap" '
    rarefaction_plot_html = 'id="rarefaction_plot"'
    gwas_table_html = '<table id="gwas_table"'

    def assertPageTitle(self, resp, title):
        self.assertContains(
            resp, f'<p class="home-title">{title}</p>', html=True)

    def assertNav(self, resp):
        self.assertContains(resp, "<nav>")

    def assertNoCompTable(self, resp):
        self.assertFalse(resp.context.get("show_comparison_table", False))
        self.assertNotContains(resp, self.table_html)

    def assertCompTable(self, resp):
        self.assertTrue(resp.context.get("show_comparison_table", False))
        self.assertContains(resp, self.table_html)

    def assertSelection(self, resp, selected=False):
        expected = self.selection_html.format(
            'selected=""' if selected else '')
        self.assertContains(resp, expected, html=True)

    def assertNoVennDiagram(self, resp):
        self.assertFalse(resp.context.get("show_results", False))
        self.assertTemplateNotUsed('chlamdb/venn_representation_template.html')
        self.assertNotContains(resp, self.venn_html)

    def assertVennDiagram(self, resp):
        self.assertTrue(resp.context.get("show_results", False))
        self.assertTemplateUsed('chlamdb/venn_representation_template.html')
        self.assertContains(resp, self.venn_html)

    def assertNoHeatmap(self, resp):
        self.assertFalse(resp.context.get("envoi_heatmap", False))
        self.assertNotContains(resp, self.heatmap_html)

    def assertHeatmap(self, resp):
        self.assertTrue(resp.context.get("envoi_heatmap", False))
        self.assertContains(resp, self.heatmap_html)

    def assertNoRarefactionPlot(self, resp):
        self.assertFalse(resp.context.get("envoi", False))
        self.assertNotContains(resp, self.rarefaction_plot_html)

    def assertRarefactionPlot(self, resp):
        self.assertTrue(resp.context.get("envoi", False))
        self.assertContains(resp, self.rarefaction_plot_html)

    def assertNoGwasTable(self, resp):
        self.assertFalse(resp.context.get("show_results", False))
        self.assertNotContains(resp, self.gwas_table_html)

    def assertGwasTable(self, resp):
        self.assertTrue(resp.context.get("show_results", False))
        self.assertContains(resp, self.gwas_table_html)

    @property
    def tab_comp_view(self):
        return f"/{self.view_type}_comparison"

    @property
    def tab_comp_form_data_list(self):
        return [{"targets": ["0", "1"]}]

    @property
    def venn_view(self):
        return f"/venn_{self.view_type}/"

    @property
    def venn_form_data(self):
        return {"targets": ["0", "1"]}

    @property
    def heatmap_view(self):
        return f"/plot_heatmap/{self.view_type}/"

    @property
    def heatmap_form_data(self):
        return {"targets": ["0", "1"]}

    @property
    def gwas_view(self):
        return f"/gwas_{self.view_type}/"

    @property
    def gwas_form_data(self):
        phenotype = SimpleUploadedFile(
            "phenotype.csv", b"1,1\n2,1\n3,0", content_type="text/csv")
        return {"max_number_of_hits": "all",
                "bonferroni_cutoff": 0.99,
                "phenotype_file": phenotype}

    def test_comparison_index_view(self):
        resp = self.client.get(f"/index_comp/{self.view_type}")
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/index_comp.html')
        self.assertPageTitle(resp, self.page_title)
        # Check that all links work
        links = re.findall('location.href="(.*?)";', str(resp.content))
        for link in links:
            resp = self.client.get(link)
            self.assertEqual(200, resp.status_code)
            # Make sure we landed on the desired view
            self.assertEqual(link.lstrip("/").split("/")[0],
                             resp.resolver_match.view_name)
            self.assertPageTitle(resp, self.page_title)

    def test_tabular_comparison_view(self):
        resp = self.client.get(self.tab_comp_view)
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/tabular_comparison.html')
        self.assertPageTitle(resp, self.page_title)
        self.assertSelection(resp)
        self.assertNoCompTable(resp)
        self.assertNav(resp)

        for tab_comp_form_data in self.tab_comp_form_data_list:
            resp = self.client.post(self.tab_comp_view,
                                    data=tab_comp_form_data)
            self.assertEqual(200, resp.status_code)
            self.assertTemplateUsed(resp, 'chlamdb/tabular_comparison.html')
            self.assertPageTitle(resp, self.page_title)
            self.assertSelection(resp, selected=True)
            self.assertCompTable(resp)
            self.assertNav(resp)

        maybe_dump_html(resp, "with_results")

    def test_venn_view(self):
        resp = self.client.get(self.venn_view)
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/venn_generic.html')
        self.assertPageTitle(resp, self.page_title)
        self.assertSelection(resp)
        self.assertNoVennDiagram(resp)
        self.assertNav(resp)

        resp = self.client.post(self.venn_view,
                                data=self.venn_form_data)
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/venn_generic.html')
        self.assertPageTitle(resp, self.page_title)
        self.assertSelection(resp, selected=True)
        self.assertVennDiagram(resp)
        self.assertNav(resp)

        maybe_dump_html(resp, "with_results")

    def test_plot_heatmap_view(self):
        resp = self.client.get(self.heatmap_view)
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/plot_heatmap.html')
        self.assertPageTitle(resp, self.page_title)
        self.assertSelection(resp)
        self.assertNoHeatmap(resp)
        self.assertNav(resp)

        resp = self.client.post(self.heatmap_view,
                                data=self.heatmap_form_data)
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/plot_heatmap.html')
        self.assertPageTitle(resp, self.page_title)
        self.assertSelection(resp, selected=True)
        self.assertHeatmap(resp)
        self.assertNav(resp)

        maybe_dump_html(resp, "with_results")

    def test_pan_genome_view(self):
        resp = self.client.get(f"/pan_genome/{self.view_type}")
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/pan_genome.html')
        self.assertEqual(self.view_type, resp.context["object_type"])
        self.assertPageTitle(resp, self.page_title)
        self.assertNoRarefactionPlot(resp)

        resp = self.client.post(f"/pan_genome/{self.view_type}",
                                data={"targets": ["0", "1"]})
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/pan_genome.html')
        self.assertEqual(self.view_type, resp.context["object_type"])
        self.assertPageTitle(resp, self.page_title)
        self.assertRarefactionPlot(resp)

        maybe_dump_html(resp, "with_results")

    def test_gwas_view(self):
        resp = self.client.get(self.gwas_view)
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/gwas.html')
        self.assertPageTitle(resp, self.page_title)
        self.assertNoGwasTable(resp)
        self.assertNav(resp)
        self.assertNotContains(resp, 'id="error"')

        resp = self.client.post(self.gwas_view,
                                data=self.gwas_form_data)
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/gwas.html')
        self.assertPageTitle(resp, self.page_title)
        # No results because there is no hit.
        self.assertNoGwasTable(resp)
        self.assertNav(resp)
        self.assertContains(resp, 'id="error"')
        self.assertContains(resp, 'No significant association')

        maybe_dump_html(resp, "with_results")


class TestPfamViews(SimpleTestCase, ComparisonViewsTestMixin):

    view_type = "pfam"
    page_title = "Comparisons: PFAM domains"

    pass


class TestKOViews(SimpleTestCase, ComparisonViewsTestMixin):

    view_type = "ko"
    page_title = "Comparisons: Kegg Orthologs (KO)"

    def test_venn_subset_view(self):
        subset_view = f"/{self.view_type}_venn_subset/Cofactor+and+vitamin+metabolism?h=1&h=2"
        resp = self.client.get(subset_view)
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/venn_generic.html')
        self.assertPageTitle(resp, self.page_title)
        self.assertSelection(resp)
        self.assertVennDiagram(resp)
        self.assertNav(resp)

        maybe_dump_html(resp)


class TestCOGViews(SimpleTestCase, ComparisonViewsTestMixin):

    view_type = "cog"
    page_title = "Comparisons: Clusters of Orthologous groups (COGs)"

    def test_venn_subset_view(self):
        subset_view = f"/{self.view_type}_venn_subset/L?h=1&h=2"
        resp = self.client.get(subset_view)
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/venn_generic.html')
        self.assertPageTitle(resp, self.page_title)
        self.assertSelection(resp)
        self.assertVennDiagram(resp)
        self.assertNav(resp)

        maybe_dump_html(resp)


class TestAMRViews(SimpleTestCase, ComparisonViewsTestMixin):

    view_type = "amr"
    page_title = "Comparisons: Antimicrobial Resistance"

    @property
    def tab_comp_form_data_list(self):
        """Tests for class and subclass are not worth much as none
        of the AMRs in the DB have a class or subclass...
        """
        return [{"targets": ["0", "1"], "comp_type": "gene"},
                {"targets": ["0", "1"], "comp_type": "class"},
                {"targets": ["0", "1"], "comp_type": "subclass"}]

    @skip("Heatmap plot fails because the test data does not provide enough hits")
    def test_plot_heatmap_view(self):
        super(TestAMRViews, self).test_plot_heatmap_view()


class TestVFViews(SimpleTestCase, ComparisonViewsTestMixin):

    view_type = "vf"
    page_title = "Comparisons: Virulence Factors"

    @property
    def tab_comp_form_data_list(self):
        return [{"targets": ["0", "1"], "comp_type": "vf_gene_id"},
                {"targets": ["0", "1"], "comp_type": "vfid"},
                {"targets": ["0", "1"], "comp_type": "category"},]


class TestOrthogroupViews(SimpleTestCase, ComparisonViewsTestMixin):

    view_type = "orthogroup"
    page_title = "Comparisons: orthologous groups"

    pass
