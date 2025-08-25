import json
from contextlib import contextmanager

from django.conf import settings
from django.test import SimpleTestCase
from lib.db_utils import DB


class BaseAutocompleteTestCase(SimpleTestCase):
    def make_request(self, **kwargs):
        url = f"{self.base_url}?forward={json.dumps(kwargs)}"
        return self.client.get(url)


class TestAutocompleteTaxid(BaseAutocompleteTestCase):
    base_url = "/autocomplete_taxid/"

    taxons = [
        {"id": "1", "text": "Klebsiella pneumoniae R6724_16313"},
        {"id": "2", "text": "Klebsiella pneumoniae R6726_16314"},
        {"id": "3", "text": "Klebsiella pneumoniae R6728_16315"},
    ]

    groups = [
        {"id": "group:positive", "text": "positive"},
        {"id": "group:negative", "text": "negative"},
        {"id": "group:all", "text": "all"},
    ]

    plasmids = [
        {"id": "plasmid:1", "text": "Klebsiella pneumoniae R6724_16313 plasmid"},
        {"id": "plasmid:2", "text": "Klebsiella pneumoniae R6726_16314 plasmid"},
        {"id": "plasmid:3", "text": "Klebsiella pneumoniae R6728_16315 plasmid"},
    ]

    @contextmanager
    def add_plasmid_for_taxids(self, taxids):
        """We need to commit this change so that it is picked up
        by the autocomplete view, so we need to cleanup afterwards.
        I guess we could also have isolated the DB by making a backup in setUp
        and restoring it in tearDown, but seemed like overkill for now.
        """
        try:
            plasmid_term_id = self.db.server.adaptor.execute_one(
                "SELECT term_id FROM term WHERE name='plasmid'"
            )[0]
            for taxid in taxids:
                self.db.server.adaptor.execute(
                    f"UPDATE bioentry_qualifier_value SET value=1 "
                    f"WHERE bioentry_id={taxid} AND term_id={plasmid_term_id};"
                )
                self.db.server.commit()
            yield
        finally:
            for taxid in taxids:
                self.db.server.adaptor.execute(
                    f"UPDATE bioentry_qualifier_value SET value=0 "
                    f"WHERE bioentry_id={taxid} AND term_id={plasmid_term_id};"
                )
                self.db.server.commit()

    def assertItemsEqual(self, expected, actual):
        self.assertEqual(
            sorted(expected, key=lambda x: x["id"]),
            sorted(actual, key=lambda x: x["id"]),
        )

    def setUp(self):
        biodb_path = settings.BIODB_DB_PATH
        self.db = DB.load_db_from_name(biodb_path)

    def tearDown(self):
        self.db.server.close()

    def test_handles_include_plasmids(self):
        with self.add_plasmid_for_taxids([1, 3]):
            resp = self.make_request()
            self.assertItemsEqual(self.taxons + self.groups, resp.json()["results"])

            resp = self.make_request(include_plasmids=True)
            self.assertItemsEqual(
                self.taxons + self.groups + self.plasmids[:2], resp.json()["results"]
            )

    def test_handles_exclude(self):
        resp = self.make_request(exclude=["3"])
        self.assertItemsEqual(
            [self.taxons[0], self.taxons[1], self.groups[0]], resp.json()["results"]
        )

        resp = self.make_request(exclude=["group:positive"])
        self.assertItemsEqual([self.taxons[2], self.groups[1]], resp.json()["results"])

        resp = self.make_request(exclude=["3", "group:positive"])
        self.assertItemsEqual([], resp.json()["results"])

    def test_handles_exclude_taxids_in_groups(self):
        # ignored because these are not groups
        resp = self.make_request(exclude_taxids_in_groups=["1", "3"])
        self.assertItemsEqual(self.taxons + self.groups, resp.json()["results"])

        resp = self.make_request(exclude_taxids_in_groups=["group:positive"])
        self.assertItemsEqual([self.taxons[2]] + self.groups, resp.json()["results"])

        resp = self.make_request(
            exclude_taxids_in_groups=["group:positive", "group:negative"]
        )
        self.assertItemsEqual(self.groups, resp.json()["results"])


class TestAutocompleteNMissing(BaseAutocompleteTestCase):
    base_url = "/autocomplete_n_missing/"

    @staticmethod
    def get_expected_response(n):
        return [{"id": i, "text": i} for i in range(n)]

    def test_handles_include_plasmids(self):
        included = ["1", "2", "plasmid:2", "plasmid:3"]
        resp = self.make_request(included=included)
        self.assertEqual(self.get_expected_response(2), resp.json()["results"])

        resp = self.make_request(included=included, include_plasmids=True)
        self.assertEqual(self.get_expected_response(4), resp.json()["results"])

    def test_handles_groups(self):
        resp = self.make_request(included=["group:positive"])
        self.assertEqual(self.get_expected_response(2), resp.json()["results"])

        resp = self.make_request(included=["group:positive", "group:negative"])
        self.assertEqual(self.get_expected_response(3), resp.json()["results"])
