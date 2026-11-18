"""Unit tests for primers.py.  Run from the project folder:  python -m unittest discover tests -v"""
import os
import random
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import primers as P


class TestChecks(unittest.TestCase):
    def test_revcomp(self):
        self.assertEqual(P.revcomp("ATGC"), "GCAT")
        self.assertEqual(P.revcomp(P.revcomp("GACGTCAGC")), "GACGTCAGC")

    def test_gc_percent(self):
        self.assertEqual(P.gc_percent("GGCC"), 100.0)
        self.assertEqual(P.gc_percent("ATGC"), 50.0)

    def test_self_dimer(self):
        self.assertEqual(P.longest_complementary("GAATTC", "GAATTC"), 6)   # palindrome pairs with itself
        self.assertEqual(P.longest_complementary("AAAAAAAA", "AAAAAAAA"), 0)
        self.assertEqual(P.longest_complementary("AAAA", "TTTT"), 4)       # perfect cross-dimer

    def test_hairpin(self):
        self.assertEqual(P.hairpin_stem("GGGGAAAACCCC"), 4)
        self.assertEqual(P.hairpin_stem("ACACACACACAC"), 0)
        self.assertEqual(P.hairpin_stem("GGGCAAAAGCCC"), 4)   # GGGC ... GCCC
        self.assertEqual(P.hairpin_stem("GGGAAAAACCC"), 3)

    def test_three_prime_dimer(self):
        self.assertTrue(P.three_prime_self_match("GCCGAATTCGGC"))
        self.assertFalse(P.three_prime_self_match("ATATATATGAAC"))

    def test_gc_clamp(self):
        self.assertTrue(P.gc_clamp_ok("ATATATATATATATATATGC"))
        self.assertFalse(P.gc_clamp_ok("ATATATATATATATATGGGCC"))   # 5 G/C in the last 5
        self.assertFalse(P.gc_clamp_ok("ATATATATATATATATATGA"))    # ends in A

    def test_runs(self):
        self.assertTrue(P.has_run("GAAAAT"))
        self.assertFalse(P.has_run("GAAATC"))

    def test_idt_test_sequence_is_flagged_as_a_poor_primer(self):
        # GACGTCAGCTAGCTAGCTGATCG is a fine thermodynamics test case, but it is highly self-complementary
        chk = P.check_primer("GACGTCAGCTAGCTAGCTGATCG")
        self.assertGreater(chk["self_dimer"], P.MAX_SELF_DIMER)
        self.assertGreater(chk["hairpin"], P.MAX_HAIRPIN)

    def test_problem_messages(self):
        good = "AGAGTTTGATCCTGGCTCAG"   # the 27F 16S primer passes every check
        self.assertEqual(P.problems(P.check_primer(good)), [])
        self.assertIn("run of 4 identical bases", P.problems(P.check_primer("GCGAAAATCGATCGGCTAGC")))


class TestParsing(unittest.TestCase):
    def test_multi_record_fasta_stays_separate(self):
        # regression test: the first version merged all records into one long sequence
        records, errors = P.parse_sequences(">a\nACGT\nACGT\n>b\nGGGGCCCC\n>c\nTTTTAAAA\n")
        self.assertEqual(records, [("a", "ACGTACGT"), ("b", "GGGGCCCC"), ("c", "TTTTAAAA")])
        self.assertEqual(errors, [])

    def test_plain_list_and_cleanup(self):
        records, _ = P.parse_sequences("acgtacgt\n GGGG CCCC 12 \r\n")
        self.assertEqual([s for _, s in records], ["ACGTACGT", "GGGGCCCC"])

    def test_errors_are_reported_not_raised(self):
        records, errors = P.parse_sequences(">ok\nACGTACGT\n>rna\nACGUACGU\n>empty\n>tiny\nAC\n")
        self.assertEqual([n for n, _ in records], ["ok"])
        self.assertEqual(len(errors), 3)
        self.assertTrue(any("U" in e for e in errors))

    def test_empty(self):
        self.assertEqual(P.parse_sequences("  \n\n"), ([], []))

    def test_text_before_first_header_is_kept(self):
        records, _ = P.parse_sequences("ACGTACGT\n>b\nGGGGCCCC\n")
        self.assertEqual([s for _, s in records], ["ACGTACGT", "GGGGCCCC"])


class TestScanAndPairs(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        random.seed(5)
        cls.template = "".join(random.choice("ACGT") for _ in range(1500))
        cls.cands = P.scan_template(cls.template, 20, 40, 60, 55, 65)

    def test_candidates_are_real_and_pass_every_check(self):
        self.assertGreater(len(self.cands), 10)
        for c in self.cands:
            window = self.template[c["start"]:c["end"]]
            expected = window if c["strand"] == "fwd" else P.revcomp(window)
            self.assertEqual(c["seq"], expected)
            self.assertEqual(P.problems(P.check_primer(c["seq"])), [])
            self.assertTrue(55 <= c["tm"] <= 65)

    def test_pairs_respect_every_constraint(self):
        pairs = P.find_pairs(self.cands, 150, 600, 2.0)
        self.assertGreater(len(pairs), 0)
        last = (-1, -1)
        for p in pairs:
            f, r = p["fwd"], p["rev"]
            self.assertEqual(f["strand"], "fwd")
            self.assertEqual(r["strand"], "rev")
            self.assertLessEqual(f["end"], r["start"])                  # no overlap, forward first
            self.assertEqual(p["product"], r["end"] - f["start"])
            self.assertTrue(150 <= p["product"] <= 600)
            self.assertLessEqual(abs(f["tm"] - r["tm"]), 2.0)
            self.assertLessEqual(p["cross"], P.MAX_SELF_DIMER)
            self.assertGreaterEqual((p["dtm"], p["product"]), last)     # sorted best-first
            last = (p["dtm"], p["product"])

    def test_tighter_limit_gives_fewer_pairs(self):
        loose = P.find_pairs(self.cands, 100, 800, 3.0, limit=10**6)
        tight = P.find_pairs(self.cands, 100, 800, 0.5, limit=10**6)
        self.assertLessEqual(len(tight), len(loose))
        self.assertTrue(all(p["dtm"] <= 0.5 for p in tight))


if __name__ == "__main__":
    unittest.main()
