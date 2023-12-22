import glob
import os
import tempfile
from unittest import skip

import nextflow
from testing.base import BaseTestCase


class TestDBSetupPipeline(BaseTestCase):

    def setUp(self):
        super(TestDBSetupPipeline, self).setUp()
        basedir = os.path.dirname(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))))
        self.nf_file_path = os.path.join(basedir, "db_setup.nf")
        self.db_dir = tempfile.TemporaryDirectory()
        self.singularity_dir = tempfile.TemporaryDirectory()
        # self.nf_dir = tempfile.TemporaryDirectory()

        self.nf_params = {
            "docker": "false",
            "conda": "false",
            "singularity": "true",
            "base_db": self.db_dir.name,
            "singularity_dir": self.singularity_dir.name,
        }
        self.nf_args = {}
        # self.nf_args = {"output_path": self.nf_dir.name}

    def execute_pipeline(self):
        return nextflow.run(self.nf_file_path,
                            params=self.nf_params,
                            **self.nf_args)

    def assert_success(self, execution):
        self.assertEqual(execution.status, "OK",
                         "Pipeline execution failed")
        self.assertEqual(execution.return_code, "0",
                         "Pipeline execution failed")
        for proc in execution.process_executions:
            self.assertEqual(proc.return_code, "0",
                             "Process execution failed")
            self.assertEqual(proc.status, "COMPLETED",
                             "Process execution failed")

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

        self.assertEqual([proc.name for proc in execution.process_executions],
                         ["download_pfam_db", "prepare_hmm"])

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
            os.listdir(os.path.join(self.db_dir.name, "pfam")))

    def test_creating_cog_db(self):
        self.nf_params["cog"] = "true"
        execution = self.execute_pipeline()
        self.assert_success(execution)

        self.assertEqual([proc.name for proc in execution.process_executions],
                         ["download_cog_cdd", "setup_cog_cdd"])

        download_process = execution.process_executions[0]
        self.assertIn("Cog.pn", os.listdir(download_process.path))
        smp_files = glob.glob(os.path.join(download_process.path, "COG*.smp"))
        self.assertTrue(len(smp_files) > 1000)

        setup_process = execution.process_executions[1]
        expected_files = ["cdd_to_cog", "cog_db.aux", "cog_db.freq",
                          "cog_db.loo", "cog_db.phr", "cog_db.pin",
                          "cog_db.psi", "cog_db.psq", "cog_db.rps",
                          "cog_db.psd"]
        # Files are moved to db directory
        self.assertItemsEqual(
            expected_files,
            os.listdir(os.path.join(self.db_dir.name, "cog")))

    def test_creating_ko_db(self):
        self.nf_params["ko"] = "true"
        execution = self.execute_pipeline()
        self.assert_success(execution)

        self.assertEqual([proc.name for proc in execution.process_executions],
                         ["download_ko_profiles"])

        # Files are moved to db directory
        db_dir = os.path.join(self.db_dir.name, "kegg")
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
                         ["download_swissprot", "prepare_swissprot"])

        download_process = execution.process_executions[0]
        self.assert_created_files(download_process, ["swissprot.fasta"])

        setup_process = execution.process_executions[1]
        # Files are moved to db directory
        expected_files = ['swissprot.fasta.phr',
                          'swissprot.fasta',
                          'swissprot.fasta.pin',
                          'swissprot.fasta.psq']

        self.assertItemsEqual(
            expected_files,
            os.listdir(os.path.join(self.db_dir.name, "uniprot", "swissprot")))
