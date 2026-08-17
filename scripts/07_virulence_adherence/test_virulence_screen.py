#!/usr/bin/env python3
"""Small static unit tests for BLAST coverage and status classification."""

import importlib.util
import unittest
from pathlib import Path

import pandas as pd

MODULE_PATH = Path(__file__).with_name("virulence_screen.py")
SPEC = importlib.util.spec_from_file_location("virulence_screen", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class VirulenceHelpersTest(unittest.TestCase):
    def setUp(self):
        self.mapping = pd.DataFrame([{
            "Query_ID": "VFG1|REF1|ace|ace", "Project_target": "ace",
            "Database_gene_symbol": "ace", "Reference_accession": "REF1",
            "Database_annotation": "collagen adhesin", "Mapping_status": "approved_primary",
        }])

    def raw(self, identity, length, start=1, contig="contig1"):
        return pd.DataFrame([{
            "Genome": "GCA_029011255.1", "qseqid": "VFG1|REF1|ace|ace", "sseqid": contig,
            "pident": identity, "length": length, "qlen": 1000, "qstart": start,
            "qend": start + length - 1, "sstart": 100, "send": 100 + length - 1,
            "evalue": 1e-50, "bitscore": 500,
        }])

    def status(self, raw):
        _, statuses = MODULE.classify_candidates(raw, ["GCA_029011255.1"], self.mapping)
        return statuses.loc[statuses.Target_gene == "ace", "Detection_status"].iloc[0]

    def test_interval_union(self):
        self.assertEqual(MODULE.union_length([(1, 100), (90, 150), (200, 220)]), 171)

    def test_distant_hsps_on_one_contig_are_separate_loci(self):
        raw = pd.concat([self.raw(95, 400, start=1), self.raw(95, 400, start=401)], ignore_index=True)
        raw.loc[1, ["sstart", "send"]] = [20000, 20399]
        raw["Strand"] = "+"
        groups = MODULE.split_coordinate_loci(raw, 1000)
        self.assertEqual(len(groups), 2)

    def test_accepted_threshold(self):
        self.assertEqual(self.status(self.raw(80, 800)), "accepted_present")

    def test_partial_threshold(self):
        self.assertEqual(self.status(self.raw(90, 600)), "flagged_partial")

    def test_review_threshold(self):
        self.assertEqual(self.status(self.raw(75, 850)), "review_required")

    def test_multiple_distinct_hits_are_ambiguous(self):
        raw = pd.concat([self.raw(95, 900, contig="c1"), self.raw(95, 900, contig="c2")], ignore_index=True)
        self.assertEqual(self.status(raw), "ambiguous_multiple_hit")

    def test_multiple_references_at_same_locus_are_not_ambiguous(self):
        mapping = pd.concat([self.mapping, self.mapping.assign(
            Query_ID="VFG2|REF2|ace|ace", Reference_accession="REF2")], ignore_index=True)
        first = self.raw(95, 900)
        second = self.raw(96, 900)
        second["qseqid"] = "VFG2|REF2|ace|ace"
        raw = pd.concat([first, second], ignore_index=True)
        candidates, statuses = MODULE.classify_candidates(raw, ["GCA_029011255.1"], mapping)
        status = statuses.loc[statuses.Target_gene == "ace", "Detection_status"].iloc[0]
        self.assertEqual(status, "accepted_present")
        self.assertEqual(candidates.loc[candidates.Target_gene == "ace", "Distinct_target_locus_ID"].nunique(), 1)


if __name__ == "__main__":
    unittest.main()
