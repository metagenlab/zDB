from django.conf import settings
from django.test import SimpleTestCase
from lib.db_utils import DB

from webapp.views.utils import AccessionFieldHandler, EntryIdParser


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
        exclude = [el[0] for el in self.taxons[1:]]
        self.assertItemsEqual(
            [self.taxons[0]],
            self.handler.get_choices(exclude=exclude))

        self.add_plasmid_for_taxids([1, 2, 3])
        self.assertItemsEqual(self.taxons + self.groups + self.plasmids,
                              self.handler.get_choices())

        exclude = [self.taxons[-1][0]]
        self.assertItemsEqual(self.taxons[:-1] + [self.groups[2]] + self.plasmids,
                              self.handler.get_choices(exclude=exclude))

        self.assertItemsEqual(self.taxons[:-1] + [self.groups[2]],
                              self.handler.get_choices(exclude=exclude,
                                                       with_plasmids=False))

    def test_get_choices_handles_plasmid_exclusion(self):
        self.add_plasmid_for_taxids([1, 2, 3])

        exclude = [self.plasmids[-1][0]]
        self.assertItemsEqual(self.taxons + self.groups + self.plasmids[:-1],
                              self.handler.get_choices(exclude=exclude))

        exclude = [el[0] for el in self.plasmids]
        self.assertItemsEqual(self.taxons + self.groups,
                              self.handler.get_choices(exclude=exclude))

    def test_get_choices_handles_group_exclusion(self):
        exclude = [self.groups[1][0]]
        self.assertItemsEqual(self.taxons[:-1] + [self.groups[2]],
                              self.handler.get_choices(exclude=exclude))

        exclude.append(self.groups[2][0])
        self.assertItemsEqual([],
                              self.handler.get_choices(exclude=exclude))

    def test_get_choices_handles_exclude_taxids_in_groups(self):
        exclude = [self.groups[1][0]]
        self.assertItemsEqual(
            self.taxons[:-1] + self.groups,
            self.handler.get_choices(exclude_taxids_in_groups=exclude))

        exclude.append(self.groups[2][0])
        self.assertItemsEqual(
            self.groups,
            self.handler.get_choices(exclude_taxids_in_groups=exclude))

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


class TestEntryIdParser(SimpleTestCase):

    def setUp(self):
        biodb_path = settings.BIODB_DB_PATH
        self.db = DB.load_db_from_name(biodb_path)
        self.entry_id_parser = EntryIdParser(self.db)

    def tearDown(self):
        self.db.server.close()

    def test_handles_orthogroups(self):
        self.assertEqual(
            ('orthogroup', 1),
            self.entry_id_parser.id_to_object_type("group_1"))

    def test_handles_orthogroup_0(self):
        self.assertEqual(
            ('orthogroup', 0),
            self.entry_id_parser.id_to_object_type("group_0"))

    def test_handles_cog(self):
        self.assertEqual(
            ('cog', 1234),
            self.entry_id_parser.id_to_object_type("COG1234"))

    def test_handles_pfam(self):
        self.assertEqual(
            ('pfam', 10423),
            self.entry_id_parser.id_to_object_type("PF10423"))

    def test_handles_ko(self):
        self.assertEqual(
            ('ko', 1241),
            self.entry_id_parser.id_to_object_type("K01241"))

    def test_handles_vf(self):
        self.assertEqual(
            ('vf', "VFG048797"),
            self.entry_id_parser.id_to_object_type("VFG048797"))

    def test_handles_amr(self):
        self.assertEqual(
            ('amr', "ybtP"),
            self.entry_id_parser.id_to_object_type("ybtP"))
