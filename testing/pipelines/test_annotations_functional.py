import os
import sys
from io import StringIO
from testing.pipelines.base import BaseTestCase

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))), "bin"))

from annotations import InputHandler, InvalidInput  # noqa


class TestInputHandler(BaseTestCase):

    def assert_names(self, handler, names):
        self.assertEqual(names, [entry.name for entry in handler.csv_entries])

    def assert_files(self, handler, files):
        self.assertEqual(files, [entry.file for entry in handler.csv_entries])

    def test_handles_list_of_files(self):
        input_file = os.path.join(self.test_dir, "assets", "test_input.csv")
        handler = InputHandler(input_file)
        self.assert_names(handler, [None, None, None])
        self.assert_files(handler, ['R6724_16313_small.gbk',
                                    'R6726_16314_small.gbk',
                                    'R6728_16315_small.gbk'])

    def test_handles_name_column(self):
        input_file = StringIO("name,file\nfoo,file1.gbk\nbar,file2.gbk")
        handler = InputHandler(input_file)
        self.assertEqual(2, len(handler.csv_entries))
        self.assert_names(handler, ["foo", "bar"])
        self.assert_files(handler, ['file1.gbk', 'file2.gbk'])

    def test_handles_incomplete_name_column(self):
        input_file = StringIO("name,file\n,file1.gbk\nbar,file2.gbk")
        handler = InputHandler(input_file)
        self.assertEqual(2, len(handler.csv_entries))
        self.assert_names(handler, [None, "bar"])
        self.assert_files(handler, ['file1.gbk', 'file2.gbk'])

    def test_raises_for_dupplicate_filename(self):
        input_file = StringIO("file\nfile1.gbk\nfile1.gbk")
        with self.assertRaises(InvalidInput) as exc:
            InputHandler(input_file)
        self.assertEqual('File "file1.gbk" appears twice in the input file.',
                         str(exc.exception))

    def test_raises_for_dupplicate_name(self):
        input_file = StringIO("file,name\nfile1.gbk,foo\nfile2.gbk,foo")
        with self.assertRaises(InvalidInput) as exc:
            InputHandler(input_file)
        self.assertEqual('Name "foo" appears twice in the input file.',
                         str(exc.exception))

    def test_raises_for_invalid_column_header(self):
        input_file = StringIO("filee,name\nfile1.gbk,foo\nfile2.gbk,foo")
        with self.assertRaises(InvalidInput) as exc:
            InputHandler(input_file)
        self.assertEqual('Invalid column header "filee" in input file.',
                         str(exc.exception))