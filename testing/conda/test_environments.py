import glob
import os
import subprocess
from tempfile import TemporaryDirectory

from testing.pipelines.base import BaseTestCase


class TestCondaEnvironments(BaseTestCase):

    env_commands = {
        "amrfinderplus.yaml": "amrfinder --version",
        'annotation.yaml': 'python -c "import Bio"',
        "blast.yaml": "blastp -h",
        "checkm.yaml": "checkm -h",
        "fasttree.yaml": "fasttree -expert",
        "kofamscan.yaml": "exec_annotation -h",
        "mafft.yaml": "mafft --version",
        "main.yaml": "nextflow -v && singularity version",
        "orthofinder.yaml": "orthofinder -h",
        "pfam_scan.yaml": "pfam_scan.pl -h",
        'testing.yaml': "nextflow -v && singularity version",
        'webapp.yaml': 'python -c "import django"',
        'zdb.yaml': "zdb",
        }

    @property
    def conda_dir(self):
        basedir = os.path.dirname(os.path.dirname(self.test_dir))
        return os.path.join(basedir, "conda")

    def create_env(self, filename, envdir):
        filepath = os.path.join(self.conda_dir, filename)
        subprocess.check_call(f"mamba env create -f {filepath} -p {envdir}", shell=True)

    def run_with_env(self, command, envdir):
        subprocess.check_call(f"conda run -p {envdir} {command}", shell=True)

    def test_all_environments_are_tested(self):
        yaml_files = glob.glob("*.yaml", root_dir=self.conda_dir)
        self.assertItemsEqual(self.env_commands.keys(), yaml_files)

    def test_all_environments_work(self):
        for filename, command in self.env_commands.items():
            with TemporaryDirectory() as tempdir:
                self.create_env(filename, tempdir)
                self.run_with_env(command, tempdir)
