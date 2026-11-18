"""Tests that drive the real web page with Streamlit's AppTest (needs streamlit installed).
Run from the project folder:  python -m unittest discover tests -v"""
import os
import sys
import unittest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, ROOT)
from streamlit.testing.v1 import AppTest

APP = os.path.join(ROOT, "app.py")
PRIMER = "GACGTCAGCTAGCTAGCTGATCG"


def fresh():
    at = AppTest.from_file(APP, default_timeout=60)
    at.run()
    return at


def text_area(at, key):
    return next(w for w in at.text_area if w.key == key)


def number(at, label):
    return next(w for w in at.number_input if w.label.startswith(label))


class TestApp(unittest.TestCase):
    def test_loads_without_error(self):
        at = fresh()
        self.assertEqual(len(at.exception), 0)
        self.assertEqual(at.title[0].value, "Sequence-Based Primer Melting Curve Simulator")
        self.assertEqual(len(at.tabs), 4)

    def test_single_primer_gives_known_tm(self):
        at = fresh()
        text_area(at, "primer_text").set_value(PRIMER).run()
        self.assertEqual(len(at.exception), 0)
        row = at.dataframe[0].value.iloc[0]
        self.assertAlmostEqual(row["Tm (°C)"], 58.86, places=2)
        theta_metric = next(m for m in at.metric if m.label == "θ at Tm")
        self.assertEqual(theta_metric.value, "0.500")

    def test_magnesium_changes_tm_like_idt(self):
        at = fresh()
        text_area(at, "primer_text").set_value(PRIMER).run()
        number(at, "Mg").set_value(3.0).run()
        number(at, "dNTPs").set_value(0.8).run()
        self.assertEqual(len(at.exception), 0)
        self.assertAlmostEqual(at.dataframe[0].value.iloc[0]["Tm (°C)"], 65.87, places=2)

    def test_multi_record_fasta_gives_separate_rows(self):
        at = fresh()
        text_area(at, "primer_text").set_value(">p1\n%s\n>p2\nAGAGTTTGATCCTGGCTCAG\n>p3\nTGTAAAACGACGGCCAGT" % PRIMER).run()
        self.assertEqual(len(at.exception), 0)
        self.assertEqual(len(at.dataframe[0].value), 3)

    def test_bad_input_shows_message_not_crash(self):
        at = fresh()
        text_area(at, "primer_text").set_value("ACGUACGUACGU").run()
        self.assertEqual(len(at.exception), 0)
        self.assertTrue(any("only A, C, G, T" in w.value for w in at.warning))

    def test_template_scan_and_pairs_work_without_nested_buttons(self):
        at = fresh()
        next(b for b in at.button if b.key == "tpl_example").click().run()
        self.assertEqual(len(at.exception), 0)
        names = [w.label for w in at.slider]
        self.assertIn("Primer length", names)
        # candidates table, then the primer-pair table with its own selector
        self.assertGreaterEqual(len(at.dataframe), 2)
        pair_table = at.dataframe[-1].value
        self.assertIn("ΔTm (°C)", pair_table.columns)
        self.assertGreater(len(pair_table), 0)
        self.assertTrue((pair_table["ΔTm (°C)"] <= 2.0).all())
        pick = next(s for s in at.selectbox if s.label.startswith("Show melting curves"))
        pick.select(pair_table["Rank"].iloc[-1]).run()
        self.assertEqual(len(at.exception), 0)


if __name__ == "__main__":
    unittest.main()
