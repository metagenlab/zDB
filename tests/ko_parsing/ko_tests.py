
# as the code is run from within the annotation_pipeline container,
# this tests the version of metagenlab installed within the container
from metagenlab_libs import KO_module, db_utils


class KO_test_parser:
    def __init__(self, db):
        self.db = db

    def run(self):
        pass

    def cleanup(self):
        return


class KO_test:
    def __init__(self, expr, kos, expected):
        self.expr = expr
        self.kos = {ko: 1 for ko in kos}
        self.expected = expected

    def cleanup(self):
        return

    def run(self):
        parser = KO_module.ModuleParser(self.expr)
        tree = parser.parse()
        n_missing = tree.get_n_missing(self.kos)
        self.passed = n_missing==self.expected


ko_tests = []
