from unittest import TestCase


class BaseTestCase(TestCase):

    def assertItemsEqual(self, actual, expected):
        self.assertEqual(sorted(actual), sorted(expected))
