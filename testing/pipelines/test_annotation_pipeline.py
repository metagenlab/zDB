import os

from sqlalchemy import MetaData, create_engine
from sqlalchemy.orm import Session
from testing.base import BasePipelineTestCase

base_tables = [
    'biodatabase',
    'biodb_config',
    'bioentry',
    'bioentry_dbxref',
    'bioentry_path',
    'bioentry_qualifier_value',
    'bioentry_reference',
    'bioentry_relationship',
    'biosequence',
    'cog_functions',
    'cog_names',
    'comment',
    'dbxref',
    'dbxref_qualifier_value',
    'filenames',
    'gene_phylogeny',
    'genome_summary',
    'ko_class',
    'ko_def',
    'ko_module_def',
    'ko_pathway_def',
    'ko_to_module',
    'ko_to_pathway',
    'location',
    'location_qualifier_value',
    'og_hits',
    'ontology',
    'orthology_identity',
    'reference',
    'reference_phylogeny',
    'seqfeature',
    'seqfeature_dbxref',
    'seqfeature_path',
    'seqfeature_qualifier_value',
    'seqfeature_relationship',
    'sequence_hash_dictionnary',
    'taxon',
    'taxon_name',
    'term',
    'term_dbxref',
    'term_path',
    'term_relationship',
    'term_relationship_term',
    'term_synonym'
]


class TestAnnotationPipeline(BasePipelineTestCase):

    nf_filename = "annotation_pipeline.nf"

    @property
    def ref_db_dir(self):
        return os.path.join(self.basedir, "zdb_ref")

    @property
    def db_dir(self):
        return os.path.join(self.test_dir, "zdb", "results", "db")

    def load_db(self, execution):
        self.metadata_obj = MetaData()
        self.engine = create_engine(
            "sqlite:////" + os.path.join(self.db_dir, execution.identifier))
        self.session = Session(self.engine)
        self.metadata_obj.reflect(bind=self.engine)

    @property
    def test_dir(self):
        return os.path.dirname(os.path.abspath(__file__))

    def query(self, tablename):
        return self.session.query(self.metadata_obj.tables[tablename])

    def setUp(self):
        super(TestAnnotationPipeline, self).setUp()
        filedir = os.path.dirname(os.path.abspath(__file__))
        self.nf_params["input"] = os.path.join(
            filedir, "assets", "test_input.csv")

    def assert_created_files(self, proc, files):
        created_files = os.listdir(proc.path)
        for file in files:
            self.assertIn(file, created_files)

    def assert_db_base_table_row_counts(self):
        self.assertEqual(3, self.query("taxon").count())
        self.assertEqual(4, self.query("bioentry").count())
        self.assertEqual(637, self.query("seqfeature").count())
        self.assertEqual(304, self.query("og_hits").count())
        self.assertEqual(304, self.query("sequence_hash_dictionnary").count())
        self.assertEqual(32, self.query("term").count())

    def test_base_pipeline(self):
        execution = self.execute_pipeline()
        self.assert_success(execution)

        self.load_db(execution)

        # Let's check that tables were correctly created and filled
        self.assertItemsEqual(base_tables, self.metadata_obj.tables.keys())
        self.assert_db_base_table_row_counts()

    def test_cog_hits(self):
        self.nf_params["cog"] = "true"
        execution = self.execute_pipeline()
        self.assert_success(execution)
        self.load_db(execution)

        # Let's check that tables were correctly created and filled
        self.assertItemsEqual(base_tables + ['cog_hits'],
                              self.metadata_obj.tables.keys())
        self.assert_db_base_table_row_counts()
        self.assertEqual(189, self.query("cog_hits").count())

    def test_ko_hits(self):
        self.nf_params["ko"] = "true"
        execution = self.execute_pipeline()
        self.assert_success(execution)
        self.load_db(execution)

        # Let's check that tables were correctly created and filled
        self.assertItemsEqual(base_tables + ['module_completeness', 'ko_hits'],
                              self.metadata_obj.tables.keys())
        self.assert_db_base_table_row_counts()
        self.assertEqual(3, self.query("module_completeness").count())
        self.assertEqual(137, self.query("ko_hits").count())

    def test_pfam_hits(self):
        self.nf_params["pfam"] = "true"
        execution = self.execute_pipeline()
        self.assert_success(execution)
        self.load_db(execution)

        # Let's check that tables were correctly created and filled
        self.assertItemsEqual(base_tables + ['pfam_hits', 'pfam_table'],
                              self.metadata_obj.tables.keys())
        self.assert_db_base_table_row_counts()
        self.assertEqual(324, self.query("pfam_hits").count())
        self.assertEqual(237, self.query("pfam_table").count())

    def test_swissprot_hits(self):
        self.nf_params["blast_swissprot"] = "true"
        execution = self.execute_pipeline()
        self.assert_success(execution)
        self.load_db(execution)

        # Let's check that tables were correctly created and filled
        self.assertItemsEqual(
            base_tables + ['swissprot_defs', 'swissprot_hits'],
            self.metadata_obj.tables.keys())
        self.assert_db_base_table_row_counts()
        self.assertEqual(19400, self.query("swissprot_defs").count())
        self.assertEqual(29400, self.query("swissprot_hits").count())

    def test_full_pipeline(self):
        self.nf_params["pfam"] = "true"
        self.nf_params["ko"] = "true"
        self.nf_params["blast_swissprot"] = "true"
        self.nf_params["cog"] = "true"
        execution = self.execute_pipeline()
        self.assert_success(execution)
        self.load_db(execution)

        # Let's check that tables were correctly created and filled
        added_tables = [
            'ko_hits',
            'pfam_table',
            'pfam_hits',
            'cog_hits',
            'swissprot_hits',
            'swissprot_defs',
            'module_completeness'
        ]
        self.assertItemsEqual(base_tables + added_tables,
                              self.metadata_obj.tables.keys())
        self.assert_db_base_table_row_counts()
        self.assertEqual(189, self.query("cog_hits").count())
        self.assertEqual(3, self.query("module_completeness").count())
        self.assertEqual(137, self.query("ko_hits").count())
        self.assertEqual(324, self.query("pfam_hits").count())
        self.assertEqual(237, self.query("pfam_table").count())
        self.assertEqual(19400, self.query("swissprot_defs").count())
        self.assertEqual(29400, self.query("swissprot_hits").count())