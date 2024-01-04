from chlamdb.urls import urlpatterns
from django.test import SimpleTestCase
from django.urls import resolve

broken_views = [
    '/orthogroup_list_cog_barchart/',
    '/orthogroup_list_cog_barchart/True/',
    '/circos_main/',
    '/genomes_intro/',
    '/hydropathy/CHUV_00025'
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


class TestViews(SimpleTestCase):

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
