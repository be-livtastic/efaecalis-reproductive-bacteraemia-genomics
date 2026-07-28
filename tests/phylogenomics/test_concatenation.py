from pathlib import Path
import importlib.util
import sys

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "06_phylogenomics" / "06_concatenate_alignments.py"
sys.path.insert(0, str(SCRIPT.parent))
spec = importlib.util.spec_from_file_location("concat", SCRIPT)
concat = importlib.util.module_from_spec(spec)
spec.loader.exec_module(concat)


def test_site_counts():
    records = {"a": "ACGT", "b": "ATGT", "c": "ATGT", "d": "ACGT"}
    assert concat.site_counts(records) == (1, 1)
