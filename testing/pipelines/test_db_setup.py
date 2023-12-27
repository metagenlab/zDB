import glob
import os
import tempfile

from testing.base import BasePipelineTestCase


class TestDBSetupPipeline(BasePipelineTestCase):

    nf_filename = "db_setup.nf"
    _ref_db_dir = None

    @property
    def ref_db_dir(self):
        if self._ref_db_dir is None:
            self._ref_db_dir = tempfile.TemporaryDirectory()
        return self._ref_db_dir.name

    def assert_created_files(self, proc, files):
        created_files = os.listdir(proc.path)
        for file in files:
            self.assertIn(file, created_files)

    def test_base_pipeline(self):
        execution = self.execute_pipeline()
        self.assert_success(execution)
        self.assertEqual(execution.process_executions, [])

    def test_creating_pfam_db(self):
        self.nf_params["pfam"] = "true"
        execution = self.execute_pipeline()
        self.assert_success(execution)

        self.assertEqual(
            [proc.name for proc in execution.process_executions],
            ["setup_pfam_db:download_pfam_db", "setup_pfam_db:prepare_hmm"])

        download_process = execution.process_executions[0]
        self.assert_created_files(download_process,
                                  ["Pfam-A.hmm", "Pfam-A.hmm.dat"])

        setup_process = execution.process_executions[1]
        expected_files = ["Pfam-A.hmm", "Pfam-A.hmm.dat", "Pfam-A.hmm.h3f",
                          "Pfam-A.hmm.h3i", "Pfam-A.hmm.h3m", "Pfam-A.hmm.h3p"]
        self.assert_created_files(setup_process, expected_files)
        # Files are copied to db directory
        self.assertItemsEqual(
            expected_files,
            os.listdir(os.path.join(self.ref_db_dir, "pfam")))

    def test_creating_cog_db(self):
        self.nf_params["cog"] = "true"
        execution = self.execute_pipeline()
        self.assert_success(execution)

        self.assertEqual(
            [proc.name for proc in execution.process_executions],
            ["setup_cogg_db:download_cog_cdd", "setup_cogg_db:setup_cog_cdd"])

        download_process = execution.process_executions[0]
        self.assertIn("Cog.pn", os.listdir(download_process.path))
        smp_files = glob.glob(os.path.join(download_process.path, "COG*.smp"))
        self.assertTrue(len(smp_files) > 1000)

        expected_files = ["cdd_to_cog", "cog_db.aux", "cog_db.freq",
                          "cog_db.loo", "cog_db.phr", "cog_db.pin",
                          "cog_db.psi", "cog_db.psq", "cog_db.rps",
                          "cog_db.psd"]
        # Files are moved to db directory
        self.assertItemsEqual(
            expected_files,
            os.listdir(os.path.join(self.ref_db_dir, "cog")))

    def test_creating_ko_db(self):
        self.nf_params["ko"] = "true"
        execution = self.execute_pipeline()
        self.assert_success(execution)

        self.assertEqual([proc.name for proc in execution.process_executions],
                         ["setup_ko_db:download_ko_profiles"])

        # Files are moved to db directory
        db_dir = os.path.join(self.ref_db_dir, "kegg")
        self.assertTrue(os.path.isfile(os.path.join(db_dir, "ko_list")))
        self.assertTrue(os.path.isfile(
            os.path.join(db_dir, "profiles", "prokaryote.hal")))
        hmm_files = glob.glob(os.path.join(db_dir, "profiles", "*.hmm"))
        self.assertTrue(len(hmm_files) > 20000)

    def test_creating_swissprot_db(self):
        self.nf_params["blast_swissprot"] = "true"
        execution = self.execute_pipeline()
        self.assert_success(execution)

        self.assertEqual([proc.name for proc in execution.process_executions],
                         ["setup_swissprot_db:download_swissprot",
                          "setup_swissprot_db:prepare_swissprot"])

        download_process = execution.process_executions[0]
        self.assert_created_files(download_process, ["swissprot.fasta"])

        # Files are moved to db directory
        expected_files = ['swissprot.fasta.phr',
                          'swissprot.fasta',
                          'swissprot.fasta.pin',
                          'swissprot.fasta.psq']

        self.assertItemsEqual(
            expected_files,
            os.listdir(os.path.join(self.ref_db_dir, "uniprot", "swissprot")))
