"""Optional check (needs:  pip install biopython).  Compares thermo.py with Biopython's own
nearest-neighbor engine on 300 random primers, using the SAME SantaLucia 1998 table.
Run from the project folder:  python tests/validate_against_biopython.py"""
import math
import os
import random
import statistics
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import thermo
from Bio.SeqUtils import MeltingTemp as mt

# Biopython's DNA_NN4 is the 2004 refinement. Build the 1998 unified table inside Biopython's own engine.
nn98 = dict(mt.DNA_NN4)
nn98["AA/TT"] = (-7.9, -22.2)
nn98["init"] = (0.0, 0.0)
nn98["init_A/T"] = (2.3, 4.1)
nn98["init_G/C"] = (0.1, -2.8)
nn98["init_oneG/C"] = (0.0, 0.0)
nn98["init_allA/T"] = (0.0, 0.0)
nn98["init_5T/A"] = (0.0, 0.0)

def bio(table, seq, na, mg, dn, nM, corr):
    return mt.Tm_NN(seq, nn_table=table, Na=na, Mg=mg, dNTPs=dn, dnac1=nM/2, dnac2=nM/2, saltcorr=corr)

random.seed(7)
seqs = ["".join(random.choice("ACGT") for _ in range(random.randint(18, 30))) for _ in range(300)]

def report(label, table, na, mg, dn, nM, corr):
    d = [thermo.calc_tm(s, na, nM, mg, dn)["tm"] - bio(table, s, na, mg, dn, nM, corr) for s in seqs]
    print("%-52s mean %+.4f  max|diff| %.4f C" % (label, statistics.mean(d), max(abs(x) for x in d)))

print("A) same 1998 parameters inside Biopython's engine (300 random 18-30mers)")
for na in (20, 50, 100, 300):
    report("Na=%d, no Mg" % na, nn98, na, 0, 0, 250, 5)
for nM in (50, 1000):
    report("Na=50, no Mg, oligo %d nM" % nM, nn98, 50, 0, 0, nM, 5)
for na, mg, dn in ((50, 1.5, 0), (50, 3, 0.8), (50, 3, 0), (20, 3, 0.8), (100, 3, 0.8), (50, 5, 2.0)):
    ratio = math.sqrt(mg/1000)/(na/1000)
    report("Na=%d Mg=%.1f dNTP=%.1f (R=%.2f)" % (na, mg, dn, ratio), nn98, na, mg, dn, 250, 7)
print()
print("B) sensitivity: our 1998 engine vs Biopython's own 2004 table (parameter-set choice only)")
for na, mg, dn, corr in ((50, 0, 0, 5), (50, 3, 0.8, 7)):
    report("Na=%d Mg=%.1f dNTP=%.1f vs DNA_NN4 (2004)" % (na, mg, dn), mt.DNA_NN4, na, mg, dn, 250, corr)
