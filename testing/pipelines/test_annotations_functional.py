import os
import sys
from io import StringIO

from testing.pipelines.base import BaseTestCase

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))), "bin"))

from annotations import InputHandler, InvalidInput  # noqa


class TestInputHandler(BaseTestCase):

    def setUp(self):
        self._cwd = os.getcwd()
        os.chdir(os.path.join(self.test_dir, "assets"))

    def tearDown(self):
        os.chdir(self._cwd)

    def assert_names(self, handler, names):
        self.assertEqual(names, [entry.name for entry in handler.csv_entries])

    def assert_files(self, handler, files):
        self.assertEqual(files, [entry.file for entry in handler.csv_entries])

    def assert_groups(self, handler, groups):
        self.assertEqual(groups, [entry.groups for entry in handler.csv_entries])

    def test_handles_list_of_files(self):
        input_file = "test_input.csv"
        handler = InputHandler(input_file)
        self.assert_names(handler, [None, None, None])
        self.assert_files(handler, ['R6724_16313_small.gbk',
                                    'R6726_16314_small.gbk',
                                    'R6728_16315_small.gbk'])

    def test_handles_name_column(self):
        input_file = StringIO("name, file\nfoo,R6724_16313_small.gbk\nbar, R6726_16314_small.gbk")
        handler = InputHandler(input_file)
        self.assertEqual(2, len(handler.csv_entries))
        self.assert_names(handler, ["foo", "bar"])
        self.assert_files(handler, ['R6724_16313_small.gbk', 'R6726_16314_small.gbk'])

    def test_handles_incomplete_name_column(self):
        input_file = StringIO("name,file\n,R6724_16313_small.gbk\nbar,R6726_16314_small.gbk")
        handler = InputHandler(input_file)
        self.assertEqual(2, len(handler.csv_entries))
        self.assert_names(handler, [None, "bar"])
        self.assert_files(handler, ['R6724_16313_small.gbk', 'R6726_16314_small.gbk'])

    def test_handles_group_columns(self):
        input_file = StringIO("name,file,group:first one,group:second_group\n"
                              "foo,R6724_16313_small.gbk,1,\n"
                              "bar,R6726_16314_small.gbk,false,True\n")
        handler = InputHandler(input_file)
        self.assertEqual(2, len(handler.csv_entries))
        self.assert_names(handler, ["foo", "bar"])
        self.assert_files(handler, ['R6724_16313_small.gbk', 'R6726_16314_small.gbk'])
        self.assert_groups(handler, [["first one"], ["second_group"]])

    def test_raises_for_missing_file(self):
        input_file = StringIO("file\nfoo.gbk\nR6724_16313_small.gbk")
        with self.assertRaises(InvalidInput) as exc:
            InputHandler(input_file)
        self.assertEqual('File "foo.gbk" cannot be accessed. Please check that'
                         ' foo.gbk exists and is accessible.',
                         str(exc.exception))

    def test_raises_for_dupplicate_filename(self):
        input_file = StringIO("file\nR6724_16313_small.gbk\nR6724_16313_small.gbk")
        with self.assertRaises(InvalidInput) as exc:
            InputHandler(input_file)
        self.assertEqual('File "R6724_16313_small.gbk" appears twice in the input file.',
                         str(exc.exception))

    def test_raises_for_dupplicate_name(self):
        input_file = StringIO("file,name\nR6724_16313_small.gbk,foo\nR6726_16314_small.gbk,foo")
        with self.assertRaises(InvalidInput) as exc:
            InputHandler(input_file)
        self.assertEqual('Name "foo" appears twice in the input file.',
                         str(exc.exception))

    def test_raises_for_invalid_column_header(self):
        input_file = StringIO("filee,name\nR6724_16313_small.gbk,foo\nR6726_16314_small.gbk,foo")
        with self.assertRaises(InvalidInput) as exc:
            InputHandler(input_file)
        self.assertEqual('Invalid column header "filee" in input file.',
                         str(exc.exception))

    def test_is_header_allowed(self):
        allowed_headers = ["file", "name", "group:foo", "group:foo bar",
                           "group:foo_bar", "group:foo1"]
        for header in allowed_headers:
            self.assertTrue(InputHandler.is_header_allowed(header))

        disallowed_headers = ["filee", "names", "group foo", "group-foo",
                              "foo:bar"]
        for header in disallowed_headers:
            self.assertFalse(InputHandler.is_header_allowed(header))

    def test_is_in_group(self):
        for value in ["true", "True", "Y", "yes", "x", "1"]:
            self.assertTrue(InputHandler.is_in_group({"g1": value}, "g1"))
        for value in ["false", "False", "n", "No", "0", ""]:
            self.assertFalse(InputHandler.is_in_group({"g1": value}, "g1"))
        with self.assertRaises(InvalidInput) as exc:
            InputHandler.is_in_group({"g1": "maybe"}, "g1")
        self.assertEqual(
            'Invalid entry "maybe" in group column. Should be one of '
            '"yes", "true", "x", "1", "no", "false", "0" or empty',
            str(exc.exception))
