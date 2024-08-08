import logging
import os
import sys

logger = logging.getLogger("webapp testing")

base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
testing_dir = os.path.join(base_dir, "testing")
results_dir = os.path.join(testing_dir, "pipelines", "zdb", "results")

comp_run_path = os.path.join(results_dir, ".completed_runs", "_webapp_testing")
if not os.path.exists(comp_run_path):
    logger.error(f"Could not find completed run {comp_run_path}")
    logger.error("Please run TestAnnotationPipeline.test_full_pipeline first.")
    sys.exit(1)

try:
    with open(comp_run_path, "r") as fhandler:
        run_name = fhandler.read().strip()
except OSError as exc:
    logger.error(f"Could not find completed run {comp_run_path}")
    logger.error("Please run TestAnnotationPipeline.test_full_pipeline first.")
    raise exc

os.environ["RUN_NAME"] = run_name

from settings.settings import *  # noqa

BIODB_DB_PATH = results_dir + "/db/" + run_name
SEARCH_INDEX = results_dir + "/search_index/" + run_name
BLAST_DB_PATH = results_dir + "/blast_DB/" + run_name
ASSET_ROOT = os.path.join(testing_dir, "webapp", "served_assets")
STATICFILES_DIRS = (ASSET_ROOT,)
