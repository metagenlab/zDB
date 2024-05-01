from django.test import SimpleTestCase

from webapp.views.utils import AccessionFieldHandler


class TestAccessionFieldHandler(SimpleTestCase):

    taxons = [('1', 'Klebsiella pneumoniae R6724_16313'),
              ('2', 'Klebsiella pneumoniae R6726_16314'),
              ('3', 'Klebsiella pneumoniae R6728_16315')]

    plasmids = [('plasmid:1', 'Klebsiella pneumoniae R6724_16313 plasmid'),
                ('plasmid:2', 'Klebsiella pneumoniae R6726_16314 plasmid'),
                ('plasmid:3', 'Klebsiella pneumoniae R6728_16315 plasmid')]

    groups = [('group:all', 'all'),
              ('group:negative', 'negative'),
              ('group:positive', 'positive')]

    def setUp(self):
        self.handler = AccessionFieldHandler()
        # Because we will not commit, we need to make all modification
        # on the handler's database
        self.db = self.handler.db

    def tearDown(self):
        self.db.server.close()

    def add_plasmid_for_taxids(self, taxids):
        plasmid_term_id = self.db.server.adaptor.execute_one(
            "SELECT term_id FROM term WHERE name='plasmid'")[0]
        for taxid in taxids:
            self.db.server.adaptor.execute(
                f"UPDATE bioentry_qualifier_value SET value=1 "
                f"WHERE bioentry_id={taxid} AND term_id={plasmid_term_id};")

    def assertItemsEqual(self, expected, choices):
        self.assertEqual(sorted(expected), sorted(choices))

    def test_get_choices_handles_plasmids(self):
        self.assertItemsEqual(
            self.taxons + self.groups,
            self.handler.get_choices())

        self.add_plasmid_for_taxids([1, 3])
        self.assertItemsEqual(
            self.taxons + self.groups + self.plasmids[::2],
            self.handler.get_choices())

        self.assertItemsEqual(
            self.taxons + self.groups,
            self.handler.get_choices(with_plasmids=False))

    def test_get_choices_handles_groups(self):
        self.assertItemsEqual(
            self.taxons + self.groups,
            self.handler.get_choices())

    def test_get_choices_handles_taxid_exclusion(self):
        to_exclude = [el[0] for el in self.taxons[1:]]
        self.assertItemsEqual(
            self.taxons[:1] + self.groups,
            self.handler.get_choices(to_exclude=to_exclude))

        self.assertItemsEqual(
            self.taxons[:1] + self.groups,
            self.handler.get_choices(to_exclude=to_exclude))

        self.add_plasmid_for_taxids([1, 2, 3])
        self.assertItemsEqual(self.taxons + self.groups + self.plasmids,
                              self.handler.get_choices())

        to_exclude = [self.taxons[-1][0]]
        self.assertItemsEqual(self.taxons[:-1] + self.groups + self.plasmids,
                              self.handler.get_choices(to_exclude=to_exclude))
        self.assertItemsEqual(self.taxons[:-1] + self.groups,
                              self.handler.get_choices(to_exclude=to_exclude,
                                                       with_plasmids=False))

    def test_get_choices_handles_plasmid_exclusion(self):
        self.add_plasmid_for_taxids([1, 2, 3])

        to_exclude = [self.plasmids[-1][0]]
        self.assertItemsEqual(self.taxons + self.groups + self.plasmids[:-1],
                              self.handler.get_choices(to_exclude=to_exclude))

        to_exclude = [el[0] for el in self.plasmids]
        self.assertItemsEqual(self.taxons + self.groups,
                              self.handler.get_choices(to_exclude=to_exclude))

    def test_get_choices_handles_group_exclusion(self):
        to_exclude = [self.groups[1][0]]
        self.assertItemsEqual(self.taxons[:-1] + [self.groups[0], self.groups[2]],
                              self.handler.get_choices(to_exclude=to_exclude))

        to_exclude.append(self.groups[2][0])
        self.assertItemsEqual([self.groups[0]],
                              self.handler.get_choices(to_exclude=to_exclude))

    def test_extract_choices_returns_none_when_include_plasmids_is_false(self):
        self.assertEqual(
            ([1, 3], None),
            self.handler.extract_choices(["1", "3"], False))

        self.assertEqual(
            ([1, 3], []),
            self.handler.extract_choices(["1", "3"], True))

    def test_extract_choices_handles_plasmids(self):
        self.assertEqual(
            ([2], [1, 3]),
            self.handler.extract_choices(["plasmid:1", "2", "plasmid:3"], True))

        self.assertEqual(
            ([2], None),
            self.handler.extract_choices(["plasmid:1", "2", "plasmid:3"], False))

    def test_extract_choices_handles_groups(self):
        self.assertEqual(
            ([3], []),
            self.handler.extract_choices(["group:negative"], True))

        self.assertEqual(
            ([2, 3], [1]),
            self.handler.extract_choices(["plasmid:1", "2", "group:negative"], True))
