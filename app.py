"""
app.py -- Streamlit web app: Sequence-Based Primer Melting Curve Simulator

Streamlit re-runs this whole file from top to bottom every time the user touches a
widget, so the app is written as one top-to-bottom script.  The maths lives in
thermo.py, the primer checks live in primers.py; this file only draws the page.
"""

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from primers import MAX_TEMPLATE, check_primer, find_pairs, parse_sequences, problems, scan_template
from thermo import peak_temperature, simulate, theta_at

EXAMPLE_PRIMERS = (
    ">27F\nAGAGTTTGATCCTGGCTCAG\n"
    ">1492R\nGGTTACCTTGTTACGACTT\n"
    ">M13F(-21)\nTGTAAAACGACGGCCAGT\n"
    ">M13R\nCAGGAAACAGCTATGAC"
)   # widely used primers; the checks pass some and flag others
# Synthetic example (random bases, not a real gene) so the template tools can be tried quickly
EXAMPLE_TEMPLATE = (
    "TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCCCCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAAGTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCCCTGAAGCATTGCTTTGTGAAGAGGGACTTCAGCCAATAGACCTGCATACCGGCTCATTCTTCATGTGCAACCTAGGGAGAATGTGTACATACGCTCTTACTGCGGTCGCGTCTAATAATATACATTTGCTTCGTTGACTAGCAACCCAGGGCTATAGCTATTCCCCCCGCGGCCCACCCAGTATTCCTAACGGAGCATAAATCCCACCCGAACTAAGTTTGTCGAACCTTGGTCCAAGATCGGGACTCGGTCTCCAGGTAAGACGGGCTCATTCATAAACGTTACTAAGGGGTATAATCTTCTATTTGTGGGTGGGAACACTTAGTAGACTTGCAATCCAATTACAGCAGTCTTGTGCGC"
)
COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]

st.set_page_config(page_title="Primer Melting Curve Simulator", layout="wide")
st.title("Sequence-Based Primer Melting Curve Simulator")
st.markdown(
    "Nearest-neighbor thermodynamics gives each primer's duplex fraction θ(T). "
    "The melting peak −dθ/dT is the curve a qPCR instrument shows."
)

# ------------------------------------------------------------------ sidebar
with st.sidebar:
    st.header("Reaction conditions")
    na_mm = st.number_input("Na⁺ (mM)", min_value=1.0, max_value=1000.0, value=50.0, step=5.0)
    mg_mm = st.number_input("Mg²⁺ (mM)", min_value=0.0, max_value=50.0, value=0.0, step=0.5,
                            help="PCR buffers are usually 1.5-3 mM. 0 means no Mg²⁺ correction.")
    dntp_mm = st.number_input("dNTPs (mM)", min_value=0.0, max_value=10.0, value=0.0, step=0.1,
                              help="dNTPs bind Mg²⁺, leaving less free Mg²⁺. Only matters when Mg²⁺ > 0.")
    oligo_nm = st.number_input("Oligo concentration (nM)", min_value=1.0, max_value=10000.0, value=250.0, step=50.0)
    st.caption("IDT OligoAnalyzer defaults for comparison: Na⁺ 50, Mg²⁺ 3, dNTPs 0.8, oligo 250 nM.")


def run(seq):
    """One call to the engine with the sidebar conditions."""
    return simulate(seq, na_mM=na_mm, mg_mM=mg_mm, oligo_nM=oligo_nm, dntp_mM=dntp_mm)


def read_input(key, label, placeholder, example):
    """Text box + FASTA upload (+ example button). Returns the usable (name, sequence) records."""
    st.button("Use example", key=key + "_example",
              on_click=lambda: st.session_state.update({key + "_text": example}))
    text = st.text_area(label, key=key + "_text", height=110, placeholder=placeholder)
    uploaded = st.file_uploader("or upload a FASTA file", type=["fa", "fasta", "fna", "txt"], key=key + "_file")
    if uploaded is not None:
        text = uploaded.getvalue().decode("utf-8", errors="replace")
    records, errors = parse_sequences(text)
    for message in errors:
        st.warning(message)
    return records


def curve_figure(items, title="Melting curve simulation"):
    """items: list of (label, simulate() result). Solid = θ(T), dashed = −dθ/dT."""
    fig = go.Figure()
    peak_max = max(max(res["dtheta"]) for _, res in items)
    for i, (label, res) in enumerate(items):
        color = COLORS[i % len(COLORS)]
        fig.add_trace(go.Scatter(x=res["temps"], y=res["theta"], name=label + "  θ(T)",
                                 mode="lines", line=dict(color=color, width=2), yaxis="y1"))
        fig.add_trace(go.Scatter(x=res["temps"], y=res["dtheta"], name=label + "  −dθ/dT",
                                 mode="lines", line=dict(color=color, width=2, dash="dash"), yaxis="y2"))
    fig.update_layout(
        title=title, hovermode="x unified", height=480,
        xaxis=dict(title="Temperature (°C)"),
        yaxis=dict(title="Duplex fraction θ(T)", side="left", range=[0, 1.02]),
        yaxis2=dict(title="Melting peak −dθ/dT (per °C)", side="right", overlaying="y1", showgrid=False,
                    range=[0, peak_max * 1.08]),
        legend=dict(orientation="h", y=-0.2),
    )
    return fig


@st.cache_data(show_spinner=False)
def cached_scan(template, length, gc_lo, gc_hi, tm_lo, tm_hi, na, mg, oligo, dntp):
    return scan_template(template, length, gc_lo, gc_hi, tm_lo, tm_hi, na, mg, oligo, dntp)


def table(df, **options):
    st.dataframe(df, hide_index=True, **options)


def csv_button(df, filename):
    st.download_button("Download table (CSV)", df.to_csv(index=False), file_name=filename, mime="text/csv")


tab_single, tab_scan, tab_pairs, tab_about = st.tabs(
    ["Analyze primers", "Scan a template", "Primer pairs", "How it works"])

# ------------------------------------------------------- tab 1: analyze primers
with tab_single:
    records = read_input("primer", "Primer sequence(s), 5'→3' DNA - one per line, or FASTA",
                         "GACGTCAGCTAGCTAGCTGATCG", EXAMPLE_PRIMERS)
    if not records:
        st.info("Enter one or more primer sequences above (or press Use example).")
    else:
        rows, results = [], {}
        for name, seq in records:
            res = run(seq)
            chk = check_primer(seq)
            results[name] = res
            rows.append({
                "Name": name, "Sequence": seq, "Length": len(seq),
                "GC%": round(chk["gc"], 1), "Tm (°C)": res["tm"],
                "Peak (°C)": round(peak_temperature(res), 1),
                "Notes": "; ".join(problems(chk)) or "passes all checks",
                "ΔH° (kcal/mol)": res["dh"], "ΔS° (cal/mol·K)": res["ds"],
                "GC clamp": "ok" if chk["clamp"] else "weak",
                "Self-dimer (bp)": chk["self_dimer"], "Hairpin stem (bp)": chk["hairpin"],
            })
        df = pd.DataFrame(rows)
        names = [r["Name"] for r in rows]
        shown = st.multiselect("Primers to plot (up to 6)", names, default=names[:3], max_selections=6)
        if shown:
            st.plotly_chart(curve_figure([(n, results[n]) for n in shown]))
        table(df, column_config={"Notes": st.column_config.TextColumn("Notes", width="large")})
        csv_button(df, "primer_results.csv")

        if len(shown) == 1:
            row = next(r for r in rows if r["Name"] == shown[0])
            res = results[shown[0]]
            c1, c2, c3, c4, c5 = st.columns(5)
            c1.metric("Tm (analytic)", "%.2f °C" % row["Tm (°C)"])
            c2.metric("θ at Tm", "%.3f" % theta_at(res, row["Tm (°C)"]))
            c3.metric("Peak of −dθ/dT", "%.1f °C" % row["Peak (°C)"])
            c4.metric("Peak − Tm", "%+.1f °C" % (row["Peak (°C)"] - row["Tm (°C)"]))
            c5.metric("GC content", "%.1f %%" % row["GC%"])
            st.caption("Self-check: θ must equal 0.500 at Tm (half the strands are paired). The −dθ/dT peak, "
                       "which melt-curve instruments usually report as Tm, sits slightly above Tm (about 0.5-1 °C) "
                       "because the curve is not symmetric. Compare Tm with "
                       "[IDT OligoAnalyzer](https://www.idtdna.com/calc/analyzer/) using the same Na⁺, Mg²⁺, dNTP and oligo values.")

# ------------------------------------------------------- tab 2: scan a template
candidates, template, template_name = None, None, None
with tab_scan:
    records = read_input("tpl", "Template DNA, 5'→3' (FASTA allowed, up to %d bases)" % MAX_TEMPLATE,
                         "Paste a longer sequence here...", EXAMPLE_TEMPLATE)
    if not records:
        st.info("Paste a template sequence above to find candidate primers.")
    else:
        template_name, template = records[0]
        if len(records) > 1:
            template_name = st.selectbox("The file has several records - scan which one?", [n for n, _ in records])
            template = dict(records)[template_name]
        if len(template) > MAX_TEMPLATE:
            st.error("Template is %d bases; the limit is %d to keep the app responsive." % (len(template), MAX_TEMPLATE))
            template = None
    if template:
        c1, c2, c3 = st.columns(3)
        length = c1.slider("Primer length", 18, 30, 20)
        gc_lo, gc_hi = c2.slider("GC% range", 0, 100, (40, 60))
        tm_lo, tm_hi = c3.slider("Tm range (°C)", 30, 90, (55, 65))
        if len(template) < length:
            st.error("Template is shorter than the primer length.")
        else:
            candidates = cached_scan(template, length, gc_lo, gc_hi, tm_lo, tm_hi,
                                    na_mm, mg_mm, oligo_nm, dntp_mm)
            st.write("**%d** candidate primers pass every check (GC range, no 4-base runs, GC clamp, "
                     "no self-dimer or hairpin, 3′ end, Tm range) in %s (%d bases)."
                     % (len(candidates), template_name, len(template)))
            if candidates:
                cdf = pd.DataFrame([{
                    "Strand": "forward" if c["strand"] == "fwd" else "reverse",
                    "Start": c["start"] + 1, "End": c["end"], "Primer (5'→3')": c["seq"],
                    "GC%": round(c["chk"]["gc"], 1), "Tm (°C)": round(c["tm"], 2),
                    "Self-dimer (bp)": c["chk"]["self_dimer"], "Hairpin stem (bp)": c["chk"]["hairpin"],
                } for c in candidates]).sort_values(["Start", "Strand"]).reset_index(drop=True)
                table(cdf)
                csv_button(cdf, "template_primers.csv")

                labels = ["%s %d-%d  (%.1f °C)" % (r["Strand"], r["Start"], r["End"], r["Tm (°C)"])
                          for _, r in cdf.iterrows()]
                picked = st.multiselect("Plot these primers (up to 6)", labels, default=labels[:2], max_selections=6)
                if picked:
                    seqs = dict(zip(labels, cdf["Primer (5'→3')"]))
                    st.plotly_chart(curve_figure([(lab, run(seqs[lab])) for lab in picked], "Selected candidate primers"))
            else:
                st.info("No window passes all checks. Try widening the GC or Tm range, or change the primer length.")

# ------------------------------------------------------- tab 3: primer pairs
with tab_pairs:
    if not candidates:
        st.info("Scan a template in the previous tab first; pairs are built from those candidates.")
    else:
        c1, c2, c3 = st.columns(3)
        prod_lo, prod_hi = c1.slider("Product size (bp)", 50, 2000, (100, 500), step=10)
        max_dtm = c2.slider("Largest Tm difference (°C)", 0.0, 5.0, 2.0, step=0.5)
        c3.caption("A pair is a forward and a reverse primer with a product in range, similar Tm, "
                   "and no stretch of 5+ complementary bases between the two primers.")
        pairs = find_pairs(candidates, prod_lo, prod_hi, max_dtm)
        if not pairs:
            st.info("No pair found. Widen the product size range or allow a larger Tm difference.")
        else:
            pdf = pd.DataFrame([{
                "Rank": i + 1,
                "Forward (5'→3')": p["fwd"]["seq"], "Fwd Tm (°C)": round(p["fwd"]["tm"], 2),
                "Reverse (5'→3')": p["rev"]["seq"], "Rev Tm (°C)": round(p["rev"]["tm"], 2),
                "ΔTm (°C)": round(p["dtm"], 2), "Product (bp)": p["product"],
                "Fwd start": p["fwd"]["start"] + 1, "Rev end": p["rev"]["end"],
                "Cross-dimer (bp)": p["cross"],
            } for i, p in enumerate(pairs)])
            table(pdf)
            csv_button(pdf, "primer_pairs.csv")
            rank = st.selectbox("Show melting curves for pair", list(pdf["Rank"]))
            chosen = pairs[rank - 1]
            st.plotly_chart(curve_figure([
                ("Forward", run(chosen["fwd"]["seq"])), ("Reverse", run(chosen["rev"]["seq"]))],
                "Pair %d: the two melting peaks should nearly overlap" % rank))

# ------------------------------------------------------- tab 4: explanation
with tab_about:
    st.markdown("""
**Model.** Each pair of neighboring bases contributes an enthalpy ΔH° and entropy ΔS°
(SantaLucia 1998 nearest-neighbor table), plus a term for each end of the duplex.

**Melting temperature.** Tm = ΔH° / (ΔS° + R ln(Cₜ/4)), with R = 1.987 cal/mol·K and Cₜ the oligo concentration.

**Salt.** Na⁺ shifts ΔS° (SantaLucia 1998). When Mg²⁺ is present and significant, the Owczarzy (2008)
correction is used instead, with dNTPs reducing the free Mg²⁺.

**Curve.** At temperature T the equilibrium constant is K = exp(−(ΔH° − TΔS°)/RT).
The duplex fraction θ solves K·Cₜ·(1−θ)² = 2θ. At Tm, θ = 0.5 exactly. The melting peak is −dθ/dT; its maximum sits about 0.5-1 °C above Tm.

**Primer checks.** GC 40-60 %, no run of 4 identical bases, GC clamp (3′ end is G or C, at most 3 G/C in the last 5 bases),
self-dimer and hairpin stems of at most 4 bases, and the last 4 bases at the 3′ end must not pair with the primer itself.
These are screening rules based on the longest complementary stretch, not free-energy calculations.

**Limits.** Two-state model for perfectly matched duplexes; no mismatches, no DMSO or formamide, no dNTP-free Mg²⁺
outside the Owczarzy range, and the primer is assumed to pair with a different strand (not with itself).
""")
