import os
import tempfile
from unittest import TestCase

import nextflow


class BaseTestCase(TestCase):
    def assertItemsEqual(self, actual, expected):
        self.assertEqual(sorted(actual), sorted(expected))

    @property
    def test_dir(self):
        return os.path.dirname(os.path.abspath(__file__))


class BasePipelineTestCase(BaseTestCase):
    nf_filename = None

    def setUp(self):
        super(BasePipelineTestCase, self).setUp()
        self.execution_dir = self.test_dir
        self.basedir = os.path.dirname(os.path.dirname(self.execution_dir))
        self.nf_file_path = os.path.join(self.basedir, self.nf_filename)
        self.singularity_dir = tempfile.TemporaryDirectory()

        self._ref_db_dir = os.environ.get("ZDB_TEST_REF_DIR")

        self.nf_params = {
            "docker": "false",
            "conda": "false",
            "singularity": "true",
            "base_db": self.ref_db_dir,
            "singularity_dir": self.singularity_dir.name,
        }
        self.nf_args = {"run_path": "./testing/pipelines"}

    def execute_pipeline(self):
        return nextflow.run(self.nf_file_path, params=self.nf_params, **self.nf_args)

    def assert_success(self, execution):
        self.assertEqual(execution.status, "OK", "Pipeline execution failed")
        self.assertEqual(execution.return_code, "0", "Pipeline execution failed")
        for proc in execution.process_executions:
            self.assertIn(proc.return_code, ["0", ""], "Process execution failed")
            self.assertIn(proc.status, ["COMPLETED", "-"], "Process execution failed")
