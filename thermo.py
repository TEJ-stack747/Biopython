"""
thermo.py -- Nearest-neighbor primer melting curve engine
SantaLucia J. (1998), PNAS 95(4):1460-1465 -- unified NN parameters.
"""

import math

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
INIT_GC = (0.1, -2.8)
INIT_AT = (2.3, 4.1)
R = 1.987  # cal / (mol . K)


def _nn_sum(seq):
    seq = seq.upper()
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
    na_M = na_mM / 1000.0
    return dS + 0.368 * (length - 1) * math.log(na_M)


def _tm_kelvin(dH_cal, dS, C_T):
    return dH_cal / (dS + R * math.log(C_T / 4))


def calc_tm(seq, na_mM=50.0, oligo_nM=250.0):
    dH, dS = _nn_sum(seq)
    dS = _salt_correct(dS, len(seq), na_mM)
    C_T = oligo_nM * 1e-9
    tm_k = _tm_kelvin(dH * 1000, dS, C_T)
    return {"dh": dH, "ds": dS, "tm": tm_k - 273.15}


def _theta(K, C_T):
    x = K * C_T
    if x <= 0:
        return 0.0
    A, B, C = x, -(2 * x + 2), x
    disc = max(B * B - 4 * A * C, 0.0)
    theta = (-B - math.sqrt(disc)) / (2 * A)
    return min(max(theta, 0.0), 1.0)


def melting_curve(seq, na_mM=50.0, oligo_nM=250.0, t_start=None, t_end=None, step=0.2, span=30.0):
    dH, dS = _nn_sum(seq)
    dS = _salt_correct(dS, len(seq), na_mM)
    C_T = oligo_nM * 1e-9
    dH_cal = dH * 1000

    tm_c = _tm_kelvin(dH_cal, dS, C_T) - 273.15
    if t_start is None:
        t_start = tm_c - span
    if t_end is None:
        t_end = tm_c + span

    temps, thetas = [], []
    t = t_start
    while t <= t_end:
        T_k = t + 273.15
        K = math.exp(-(dH_cal - T_k * dS) / (R * T_k))
        thetas.append(_theta(K, C_T))
        temps.append(t)
        t += step

    dtheta = [0.0] * len(thetas)
    for i in range(1, len(thetas) - 1):
        dtheta[i] = -(thetas[i + 1] - thetas[i - 1]) / (2 * step)

    return temps, thetas, dtheta


def simulate(seq, na_mM=50.0, mg_mM=0.0, oligo_nM=250.0):
    """Shared interface -- this is the function app.py calls.
    mg_mM accepted for signature compatibility; Mg2+ correction is a stretch goal, not applied yet."""
    tm_data = calc_tm(seq, na_mM, oligo_nM)
    temps, theta, dtheta = melting_curve(seq, na_mM, oligo_nM)
    return {
        "tm": round(tm_data["tm"], 2),
        "dh": round(tm_data["dh"], 2),
        "ds": round(tm_data["ds"], 2),
        "temps": temps,
        "theta": theta,
        "dtheta": dtheta,
    }


if __name__ == "__main__":
    for seq in ["GACGTCAGCTAGCTAGCTGATCG", "ATATATATATAT", "GCGCGCGCGCGC", "AGCTTAGCATGC"]:
        r = simulate(seq)
        peak_i = r["dtheta"].index(max(r["dtheta"]))
        peak_t = r["temps"][peak_i]
        print(f"{seq:25s} len={len(seq):2d}  Tm={r['tm']:6.2f}C  dH={r['dh']:6.1f}  dS={r['ds']:7.1f}  peak={peak_t:6.1f}C  diff={abs(peak_t-r['tm']):.2f}")
