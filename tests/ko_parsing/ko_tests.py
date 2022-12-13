
# as the code is run from within the annotation_pipeline container,
# this tests the version of metagenlab installed within the container
from metagenlab_libs import KO_module, db_utils


class KO_test_parser:
    """
    Opens a connection to the skeleton db present in the 
    annotation pipeline container and tries
    to parse all the KEGG modules. Sanity check
    to avoid unhandled cases at runtime.
    """
    def __init__(self):
        db_config = {"zdb.db_name": "George", "zdb.db_type": "sqlite", "zdb.db_psswd": ""}
        self.db = db_utils.DB.load_db("/home/metagenlab/zdb_base", db_config)
        self.passed = False

    def run(self):
        definitions = self.db.get_all_modules_definition()
        passed_parsing = True
        for KO_id, definition in definitions:
            try:
                mod = KO_module.ModuleParser(definition)
                mod.parse()
            except Exception as e:
                print(str(e))
                print(f"Failed on KO{KO_id}" + definition)
                passed_parsing = False
        self.passed = passed_parsing

    def cleanup(self):
        return


class KO_test:
    def __init__(self, expr, expected, kos):
        self.expr = expr
        self.kos = {ko: 1 for ko in kos}
        self.expected = expected

    def cleanup(self):
        return

    def run(self):
        parser = KO_module.ModuleParser(self.expr)
        tree = parser.parse()
        n_missing = tree.get_n_missing(self.kos)
        self.passed = (n_missing==self.expected)


ko_tests = [
    # simple and
    KO_test("K00001 K00002", 0, {1, 2}),
    KO_test("K00001 K00002", 1, {2}),
    KO_test("K00001 K00002", 1, {1}),
    KO_test("K00001 K00002", 2, {}),

    # simple or
    KO_test("K00001,K00002", 0, {1, 2}),
    KO_test("K00001,K00002", 0, {2}),
    KO_test("K00001,K00002", 0, {2}),
    KO_test("K00001,K00002", 1, {}),

    # simple complex
    KO_test("K00001+K00002", 0, {1, 2}),
    KO_test("K00001+K00002", 1, {2}),
    KO_test("K00001+K00002", 1, {1}),
    KO_test("K00001+K00002", 2, {}),

    KO_test("K00001-K00002", 0, {1, 2}),
    KO_test("K00001-K00002", 1, {2}),
    KO_test("K00001-K00002", 0, {1}),
    KO_test("K00001-K00002", 1, {}),

    # complex with optional component
    KO_test("K00001+K00002-K00003", 0, {1, 2}),
    KO_test("K00001+K00002-K00003", 0, {1, 2, 3}),
    KO_test("K00001-K00002-K00003+K00004", 0, {1, 2, 3, 4}), # OK
    KO_test("K00001-K00002-K00003+K00004", 0, {1, 3, 4}), # failed
    KO_test("K00001-K00002-K00003+K00004", 0, {1, 2, 4}), # OK
    KO_test("K00001-K00002-K00003+K00004", 1, {1, 2, 3}), # failed
    KO_test("K00001-K00002-K00003+K00004", 1, {2, 3, 4}), # OK

    # complex with two possible subunits
    KO_test("K00001+(K00004,K00005)-K00003", 2, {}),
    KO_test("K00001+(K00004,K00005)-K00003", 1, {1, 3}),
    KO_test("K00001+(K00004,K00005)-K00003", 0, {1, 3, 4}),
    KO_test("K00001+(K00004,K00005)-K00003", 0, {1, 3, 5}),

    # complex with two possible subunits
    KO_test("(K00004,K00005)+K00003", 0, {4, 3}),
    KO_test("(K00004,K00005)+K00003", 0, {5, 3}),
    KO_test("(K00004,K00005)+K00003", 1, {3}),
    KO_test("(K00004,K00005)+K00003", 1, {4, 5}),
    KO_test("(K00004,K00005)+K00003", 2, {}),

    # bizarre expression
    KO_test("K00001 K00002 -K00003", 0, {1, 2}),
    KO_test("K00001 K00002 -K00003", 0, {1, 2, 3}),
    KO_test("K00001 K00002 -K00003", 1, {1, 3}),
    KO_test("K00001 K00002 -K00003", 1, {2, 3}),
    KO_test("K00001 K00002 -K00003", 2, {}),

    # another bizarre expression
    KO_test("K00001 -- -K00003", 1, {}),
    KO_test("K00001 -- -K00003", 0, {1}),
    KO_test("K00001 -- -K00003", 0, {1, 3}),

    # or / and + parentheses combination
    KO_test("(K00001,K00002) K00003", 0, {1, 2, 3}),
    KO_test("(K00001,K00002) K00003", 0, {2, 3}),
    KO_test("(K00001,K00002) K00003", 0, {1, 3}),
    KO_test("(K00001,K00002) K00003", 1, {1, 2}),
    KO_test("(K00001,K00002) K00003", 1, {3}),
    KO_test("(K00001,K00002,K00004) K00003", 2, {}),

    # or / and + parentheses combination
    KO_test("K00003 (K00001,K00002)", 0, {1, 2, 3}),
    KO_test("K00003 (K00001,K00002)", 0, {2, 3}),
    KO_test("K00003 (K00001,K00002)", 0, {1, 3}),
    KO_test("K00003 (K00001,K00002)", 1, {1}),
    KO_test("K00003 (K00001,K00002)", 1, {1, 2}),
    KO_test("K00003 (K00001,K00002)", 2, {}),

    # connects to the db and tries to parse all modules definitions
    KO_test_parser()
]
