from django.conf import settings
from django.test import SimpleTestCase
from lib.db_utils import DB


class TestGroupsViews(SimpleTestCase):

    def test_groups_overview(self):
        resp = self.client.get("/groups/")
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/groups_overview.html')

        groups = [
            '<a href=/groups/positive>positive</a>',
            '<a href=/groups/negative>negative</a>',
            '<a href=/groups/all>all</a>']
        self.assertEqual(groups, resp.context["groups"])
        for group in groups:
            self.assertContains(resp, group)

        self.assertContains(resp, 'Add new group')

    def test_add_and_delete_group_view(self):
        resp = self.client.get("/groups/add/")
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/group_add.html')
        self.assertEqual(list(resp.context.get("form").fields.keys()),
                         ['group_name', 'genomes'])

        db = DB.load_db_from_name(settings.BIODB_DB_PATH)
        data = {"group_name": "test_group", "genomes": ["group:negative", "1"]}
        resp = self.client.post("/groups/add/", data=data)
        self.assertEqual(302, resp.status_code)
        self.assertEqual('/groups/test_group', resp.url)
        self.assertTrue(db.get_group("test_group"))

        resp = self.client.post("/groups/test_group/delete/")
        self.assertEqual(302, resp.status_code)
        self.assertEqual('/groups/', resp.url)
        self.assertFalse(db.get_group("test_group"))

    def test_group_details_view(self):
        resp = self.client.get("/groups/positive")
        self.assertEqual(200, resp.status_code)
        self.assertTemplateUsed(resp, 'chlamdb/group_details.html')
        genomes = [
            '<a href="/extract_contigs/1">Klebsiella pneumoniae R6724_16313</a>',
            '<a href="/extract_contigs/2">Klebsiella pneumoniae R6726_16314</a>']
        self.assertEqual(
            [row.accession for i, row in resp.context["genome_table"]["table_data"].iterrows()],
            genomes)
        for genome in genomes:
            self.assertContains(resp, genome, html=True)
