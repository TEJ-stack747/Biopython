"""Unit tests for thermo.py.  Run from the project folder:  python -m unittest discover tests -v"""
import os
import sys
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import thermo

# (sequence, Tm with Na 50 mM / no Mg, Tm with Na 50 / Mg 3 mM / dNTP 0.8 mM) at 250 nM oligo.
# Reference values come from Biopython's nearest-neighbor engine run with the same
# SantaLucia 1998 table (see tests/validate_against_biopython.py).
REFERENCE = [
    ("GACGTCAGCTAGCTAGCTGATCG", 58.86, 65.87),
    ("ATATATATATAT", 7.53, 13.52),
    ("GCGCGCGCGCGC", 61.90, 69.23),
    ("AGCTTAGCATGC", 36.92, 43.68),
    ("TGTAAAACGACGGCCAGT", 53.57, 60.55),
    ("CCAGGGTTTTCCCAGTCACGAC", 59.47, 66.97),
    ("AGAGTTTGATCCTGGCTCAG", 53.31, 60.62),
    ("GGATTAGATACCCTGGTAGTC", 50.31, 57.61),
]
S = "GACGTCAGCTAGCTAGCTGATCG"


class TestTm(unittest.TestCase):
    def test_matches_reference_values(self):
        for seq, tm0, tm3 in REFERENCE:
            self.assertAlmostEqual(thermo.calc_tm(seq, 50, 250)["tm"], tm0, delta=0.01, msg=seq)
            self.assertAlmostEqual(thermo.calc_tm(seq, 50, 250, 3, 0.8)["tm"], tm3, delta=0.01, msg=seq)

    def test_matches_idt_oligoanalyzer(self):
        # Readings taken from IDT OligoAnalyzer for the same sequence (250 nM oligo).
        self.assertAlmostEqual(thermo.calc_tm(S, 50, 250)["tm"], 58.86, delta=0.05)
        self.assertAlmostEqual(thermo.calc_tm(S, 50, 250, 3, 0.8)["tm"], 65.8, delta=0.15)

    def test_dh_ds_hand_calculation(self):
        # Hand-checkable case: 3-mer ACG at 1 M Na (salt term = 0).
        # dimers AC (-8.4, -22.4) + CG (-10.6, -27.2); ends A (2.3, 4.1) and G (0.1, -2.8)
        dh, ds = thermo._nn_sum("ACG")
        self.assertAlmostEqual(dh, -8.4 - 10.6 + 2.3 + 0.1, places=6)
        self.assertAlmostEqual(ds, -22.4 - 27.2 + 4.1 - 2.8, places=6)

    def test_salt_effects(self):
        self.assertLess(thermo.calc_tm(S, 20)["tm"], thermo.calc_tm(S, 50)["tm"])
        self.assertLess(thermo.calc_tm(S, 50)["tm"], thermo.calc_tm(S, 300)["tm"])
        self.assertGreater(thermo.calc_tm(S, 50, 250, 3)["tm"], thermo.calc_tm(S, 50, 250, 0)["tm"])
        # dNTPs bind Mg2+, so adding them lowers Tm
        self.assertLess(thermo.calc_tm(S, 50, 250, 3, 0.8)["tm"], thermo.calc_tm(S, 50, 250, 3, 0)["tm"])

    def test_negligible_mg_is_ignored(self):
        # sqrt(0.0005)/0.5 = 0.045 < 0.22: monovalent ions dominate, Mg2+ changes nothing
        self.assertEqual(thermo.calc_tm(S, 500, 250, 0.5)["tm"], thermo.calc_tm(S, 500, 250, 0)["tm"])

    def test_concentration_and_composition(self):
        self.assertGreater(thermo.calc_tm(S, 50, 1000)["tm"], thermo.calc_tm(S, 50, 50)["tm"])
        self.assertGreater(thermo.calc_tm("GCGCGCGCGCGC")["tm"], thermo.calc_tm("ATATATATATAT")["tm"])
        self.assertGreater(thermo.calc_tm("GACGTCAGCTAGCTAGCTGATCGACG")["tm"], thermo.calc_tm(S)["tm"])

    def test_bad_input(self):
        for bad in ("", "A", "ACGU", "ACGN", "AC GT"):
            with self.assertRaises(ValueError, msg=repr(bad)):
                thermo.calc_tm(bad)
        self.assertEqual(thermo.calc_tm("gacgtcagctagctagctgatcg")["tm"], thermo.calc_tm(S)["tm"])


class TestCurve(unittest.TestCase):
    def test_theta_is_half_at_tm(self):
        for seq, _, _ in REFERENCE:
            for mg, dn in ((0, 0), (3, 0.8)):
                r = thermo.simulate(seq, 50, mg, 250, dn)
                self.assertAlmostEqual(thermo.theta_at(r, r["tm"]), 0.5, delta=0.002, msg=seq)

    def test_theta_shape(self):
        r = thermo.simulate(S)
        th = r["theta"]
        self.assertGreater(th[0], 0.99)
        self.assertLess(th[-1], 0.01)
        self.assertTrue(all(0.0 <= x <= 1.0 for x in th))
        self.assertTrue(all(b <= a + 1e-12 for a, b in zip(th, th[1:])))   # never rises with temperature

    def test_peak_sits_just_above_tm(self):
        for seq, _, _ in REFERENCE:
            r = thermo.simulate(seq)
            offset = thermo.peak_temperature(r) - r["tm"]
            self.assertTrue(0.0 <= offset <= 1.3, "%s offset %.2f" % (seq, offset))

    def test_window_follows_tm(self):
        # low-Tm and high-Tm primers must not be cut off by a fixed window
        for seq in ("ATATATATATAT", "GCGCGCGCGCGCGCGCGCGCGCGC"):
            r = thermo.simulate(seq)
            self.assertGreater(r["theta"][0], 0.99, seq)
            self.assertLess(r["theta"][-1], 0.01, seq)
            self.assertGreater(max(r["dtheta"]), 0.0)


if __name__ == "__main__":
    unittest.main()
