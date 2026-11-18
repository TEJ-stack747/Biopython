"""
thermo.py -- Nearest-neighbor primer melting curve engine

Model
  * Nearest-neighbor (NN) thermodynamics, SantaLucia (1998) unified parameters,
    PNAS 95(4):1460-1465.
  * Tm = dH / (dS + R ln(C_T/4))   (two-state model, non-self-complementary strands)
  * Na+ correction: SantaLucia (1998), added to the entropy.
  * Mg2+ correction: Owczarzy et al. (2008), Biochemistry 47:5336-5353, applied
    to 1/Tm.  Used when Mg2+ is significant (sqrt([Mg])/[Na] >= 0.22);
    otherwise Mg2+ is treated as negligible and only the Na+ correction is used.
  * dNTPs bind Mg2+ (K = 3e4 per M), so only the free Mg2+ is used in the correction.
  * Duplex fraction theta(T) from the quadratic that follows from the
    equilibrium constant K(T) = exp(-(dH - T dS) / (R T)).
"""

import math

# (dH kcal/mol, dS cal/mol/K) for each nearest-neighbor dimer
NN_PARAMS = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2),
    "AT": (-7.2, -20.4),
    "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7),
    "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0),
    "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2),
    "GC": (-9.8, -24.4),
    "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
INIT_GC = (0.1, -2.8)   # terminal G/C initiation, added once per G/C end
INIT_AT = (2.3, 4.1)    # terminal A/T initiation, added once per A/T end
R = 1.987               # gas constant, cal / (mol K)


def check_sequence(seq):
    """Return the sequence in upper case, or raise ValueError if it is not usable DNA."""
    seq = str(seq).upper()
    if len(seq) < 2:
        raise ValueError("Sequence must be at least 2 bases long.")
    bad = sorted(set(seq) - set("ACGT"))
    if bad:
        raise ValueError("Only A, C, G, T are allowed (found: %s)." % ", ".join(bad))
    return seq


def _nn_sum(seq):
    """Add up dH and dS over every neighboring pair, plus the two terminal terms."""
    dH, dS = 0.0, 0.0
    for i in range(len(seq) - 1):
        h, s = NN_PARAMS[seq[i:i + 2]]
        dH += h
        dS += s
    for end in (seq[0], seq[-1]):
        h, s = INIT_GC if end in "GC" else INIT_AT
        dH += h
        dS += s
    return dH, dS


def _salt_correct(dS, length, na_mM):
    """SantaLucia (1998) sodium correction, applied to the entropy."""
    return dS + 0.368 * (length - 1) * math.log(na_mM / 1000.0)


def _tm_kelvin(dH_cal, dS, C_T):
    return dH_cal / (dS + R * math.log(C_T / 4))


def _owczarzy_shift(seq, na_mM, mg_mM):
    """Owczarzy (2008) magnesium term: 1/Tm(Mg) - 1/Tm(1 M Na).  Units: 1/K."""
    mon = na_mM / 1000.0
    mg = mg_mM / 1000.0
    f_gc = (seq.count("G") + seq.count("C")) / len(seq)
    a, b, c, d = 3.92e-5, -9.11e-6, 6.26e-5, 1.42e-5
    e, f, g = -4.82e-4, 5.25e-4, 8.31e-5
    ln_g = math.log(mg)
    if mon > 0:
        ln_m = math.log(mon)
        ratio = math.sqrt(mg) / mon
        if ratio < 6.0:  # Na+ still competes with Mg2+: use the mixed-salt constants
            a = 3.92e-5 * (0.843 - 0.352 * math.sqrt(mon) * ln_m)
            d = 1.42e-5 * (1.279 - 4.03e-3 * ln_m - 8.03e-3 * ln_m ** 2)
            g = 8.31e-5 * (0.486 - 0.258 * ln_m + 5.25e-3 * ln_m ** 3)
    return (a + b * ln_g + f_gc * (c + d * ln_g)
            + (e + f * ln_g + g * ln_g ** 2) / (2.0 * (len(seq) - 1)))


def _free_mg(mg_mM, dntp_mM):
    """Free Mg2+ (mM) left after dNTPs bind some of it (Mg.dNTP association constant 3e4 per M)."""
    if dntp_mM <= 0 or mg_mM <= 0:
        return mg_mM
    ka = 3.0e4
    mg, dn = mg_mM / 1000.0, dntp_mM / 1000.0
    b = ka * dn - ka * mg + 1.0
    free = (-b + math.sqrt(b * b + 4.0 * ka * mg)) / (2.0 * ka)
    return free * 1000.0


def _mg_matters(na_mM, mg_mM):
    """Owczarzy: below sqrt([Mg])/[Na] = 0.22 monovalent ions dominate and Mg2+ is ignored."""
    if mg_mM <= 0:
        return False
    if na_mM <= 0:
        return True
    return math.sqrt(mg_mM / 1000.0) / (na_mM / 1000.0) >= 0.22


def _effective(seq, na_mM, mg_mM, oligo_nM, dntp_mM=0.0):
    """Return (dH kcal/mol, effective dS cal/mol/K, C_T mol/L, Tm in kelvin).

    The salt corrections change Tm.  We keep dH fixed and put the whole shift into
    an 'effective' dS, so the same dH / dS pair drives both Tm and the curve.
    """
    seq = check_sequence(seq)
    dH, dS = _nn_sum(seq)
    C_T = oligo_nM * 1e-9
    dH_cal = dH * 1000.0
    mg_mM = _free_mg(mg_mM, dntp_mM)
    if _mg_matters(na_mM, mg_mM):
        tm_1M = _tm_kelvin(dH_cal, dS, C_T)                       # Tm in 1 M Na+
        tm_k = 1.0 / (1.0 / tm_1M + _owczarzy_shift(seq, na_mM, mg_mM))
        dS_eff = dH_cal / tm_k - R * math.log(C_T / 4)
    else:
        dS_eff = _salt_correct(dS, len(seq), na_mM)
        tm_k = _tm_kelvin(dH_cal, dS_eff, C_T)
    return dH, dS_eff, C_T, tm_k


def calc_tm(seq, na_mM=50.0, oligo_nM=250.0, mg_mM=0.0, dntp_mM=0.0):
    dH, dS, _, tm_k = _effective(seq, na_mM, mg_mM, oligo_nM, dntp_mM)
    return {"dh": dH, "ds": dS, "tm": tm_k - 273.15}


def _theta(K, C_T):
    """Duplex fraction theta for K*C_T*(1-theta)^2 = 2*theta, the root that lies between 0 and 1.

    With x = K*C_T this is x*theta^2 - (2x+2)*theta + x = 0.  Its two roots multiply to 1, so
    the small root is  theta = 2x / (2x + 2 + sqrt(8x + 4)).  Written this way (sums only) it
    stays accurate for tiny and huge x; the textbook (-b - sqrt(b^2-4ac))/2a loses digits.
    """
    x = K * C_T
    if x <= 0:
        return 0.0
    return 2.0 * x / (2.0 * x + 2.0 + math.sqrt(8.0 * x + 4.0))


def melting_curve(seq, na_mM=50.0, oligo_nM=250.0, t_start=None, t_end=None,
                  step=0.2, span=30.0, mg_mM=0.0, dntp_mM=0.0):
    dH, dS, C_T, tm_k = _effective(seq, na_mM, mg_mM, oligo_nM, dntp_mM)
    dH_cal = dH * 1000.0

    # Center the scan window on the analytic Tm so low-Tm (short, AT-rich) and
    # high-Tm (long, GC-rich) primers never get their transition cut off.
    tm_c = tm_k - 273.15
    if t_start is None:
        t_start = tm_c - span
    if t_end is None:
        t_end = tm_c + span

    temps, thetas = [], []
    n = int(round((t_end - t_start) / step))
    for i in range(n + 1):
        t = t_start + i * step
        T_k = t + 273.15
        K = math.exp(-(dH_cal - T_k * dS) / (R * T_k))
        thetas.append(_theta(K, C_T))
        temps.append(t)

    # -d(theta)/dT by central differences: the "melting peak" a qPCR machine shows
    dtheta = [0.0] * len(thetas)
    for i in range(1, len(thetas) - 1):
        dtheta[i] = -(thetas[i + 1] - thetas[i - 1]) / (2 * step)

    return temps, thetas, dtheta


def peak_temperature(result):
    """Temperature at which -d(theta)/dT is largest (the melting-peak position)."""
    i = max(range(len(result["dtheta"])), key=result["dtheta"].__getitem__)
    return result["temps"][i]


def theta_at(result, temp):
    """Duplex fraction at `temp` (deg C), by straight-line interpolation along the simulated curve."""
    t, th = result["temps"], result["theta"]
    step = t[1] - t[0]
    i = max(0, min(int((temp - t[0]) / step), len(t) - 2))
    return th[i] + (temp - t[i]) / step * (th[i + 1] - th[i])


def simulate(seq, na_mM=50.0, mg_mM=0.0, oligo_nM=250.0, dntp_mM=0.0):
    """Shared interface -- this is the function app.py calls."""
    tm_data = calc_tm(seq, na_mM, oligo_nM, mg_mM, dntp_mM)
    temps, theta, dtheta = melting_curve(seq, na_mM, oligo_nM, mg_mM=mg_mM, dntp_mM=dntp_mM)
    return {
        "tm": round(tm_data["tm"], 2),
        "dh": round(tm_data["dh"], 2),
        "ds": round(tm_data["ds"], 2),
        "temps": temps,
        "theta": theta,
        "dtheta": dtheta,
    }


if __name__ == "__main__":
    for mg in (0.0, 3.0):
        print("Na+ = 50 mM, Mg2+ = %.0f mM, oligo = 250 nM" % mg)
        for s in ["GACGTCAGCTAGCTAGCTGATCG", "ATATATATATAT", "GCGCGCGCGCGC", "AGCTTAGCATGC"]:
            r = simulate(s, mg_mM=mg)
            pk = peak_temperature(r)
            print("  %-25s len=%2d  Tm=%6.2fC  dH=%6.1f  dS=%7.1f  peak=%6.1fC  diff=%.2f"
                  % (s, len(s), r["tm"], r["dh"], r["ds"], pk, abs(pk - r["tm"])))
