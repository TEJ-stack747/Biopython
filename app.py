import streamlit as st
import plotly.graph_objects as go
from thermo import simulate
import csv
from io import StringIO

st.set_page_config(page_title="Primer Melting Curve Simulator", layout="wide")
st.title("Sequence-Based Primer Melting Curve Simulator")

st.markdown("""
**Nearest-neighbor thermodynamics** — simulates the duplex-fraction melting curve θ(T) of a PCR primer
and plots the melting peak (−dθ/dT), the same curve shape a qPCR instrument displays.
""")

with st.sidebar:
    st.header("⚙️ Parameters")
    na_mm = st.number_input("Na⁺ concentration (mM)", min_value=1, max_value=1000, value=50, step=5)
    oligo_nm = st.number_input("Oligo concentration (nM)", min_value=1, max_value=10000, value=250, step=50)
    st.markdown("---")
    st.markdown("**Stretch goal:** Mg²⁺ correction (not yet implemented)")

col1, col2 = st.columns(2)

with col1:
    st.subheader("📝 Input Sequence")
    input_type = st.radio("Input method:", ["Paste sequence", "Upload FASTA"])
    
    if input_type == "Paste sequence":
        seq_input = st.text_area("Enter primer sequence (5'→3', DNA only):", height=100, placeholder="GACGTCAGCTAGCTAGCTGATCG")
    else:
        uploaded = st.file_uploader("Upload FASTA file", type=["fa", "fasta"])
        if uploaded:
            seq_input = uploaded.read().decode("utf-8")
            lines = seq_input.split("\n")
            seq_input = "".join([line.strip() for line in lines if line and not line.startswith(">")])
        else:
            seq_input = ""

with col2:
    st.subheader("🔍 Primer Scanning (Optional)")
    template_input = st.text_area("Or enter template DNA to scan for primers:", height=100, placeholder="Longer template sequence...")
    
    if template_input.strip():
        primer_len = st.slider("Primer length", min_value=18, max_value=30, value=20)
        gc_min = st.slider("GC% min", min_value=0, max_value=100, value=40)
        gc_max = st.slider("GC% max", min_value=0, max_value=100, value=60)
        
        if st.button("Scan template for primers"):
            template = template_input.upper().replace(" ", "").replace("\n", "")
            if len(template) < primer_len:
                st.error(f"Template too short for primer length {primer_len}")
            else:
                primers_found = []
                for i in range(len(template) - primer_len + 1):
                    primer = template[i:i + primer_len]
                    if all(base in "ATGC" for base in primer):
                        gc_count = primer.count("G") + primer.count("C")
                        gc_pct = (gc_count / primer_len) * 100
                        if gc_min <= gc_pct <= gc_max:
                            skip = False
                            for base in "ATGC":
                                if base * 4 in primer:
                                    skip = True
                            if not skip:
                                primers_found.append((i, primer, gc_pct))
                
                st.write(f"Found {len(primers_found)} primers matching filter")
                if primers_found:
                    st.dataframe({
                        "Position": [p[0] for p in primers_found],
                        "Primer": [p[1] for p in primers_found],
                        "GC%": [f"{p[2]:.1f}" for p in primers_found],
                    }, use_container_width=True)
                    
                    if st.button("Analyze first 5 primers"):
                        seq_input = "\n".join([p[1] for p in primers_found[:5]])

if seq_input.strip():
    sequences = []
    labels = []
    
    lines = seq_input.strip().split("\n")
    for line in lines:
        if line.startswith(">"):
            continue
        line = line.upper().replace(" ", "").replace("\n", "")
        if line and all(base in "ATGC" for base in line):
            sequences.append(line)
            if len(sequences) == 1:
                labels.append(line[:15] + ("..." if len(line) > 15 else ""))
            else:
                labels.append(f"Seq {len(sequences)}")
    
    if not sequences:
        st.error("❌ Invalid sequence: contains non-ATGC characters")
    else:
        st.subheader("📊 Results")
        
        results_data = []
        fig = go.Figure()
        colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
        
        for idx, seq in enumerate(sequences):
            result = simulate(seq, na_mM=na_mm, oligo_nM=oligo_nm)
            
            results_data.append({
                "Sequence": seq,
                "Length": len(seq),
                "Tm (°C)": f"{result['tm']:.2f}",
                "ΔH° (kcal/mol)": f"{result['dh']:.1f}",
                "ΔS° (cal/mol·K)": f"{result['ds']:.1f}",
            })
            
            color = colors[idx % len(colors)]
            
            fig.add_trace(go.Scatter(
                x=result["temps"], y=result["theta"],
                name=f"{labels[idx]} θ(T)", mode="lines",
                line=dict(color=color, width=2), yaxis="y1",
            ))
            fig.add_trace(go.Scatter(
                x=result["temps"], y=result["dtheta"],
                name=f"{labels[idx]} −dθ/dT (peak)", mode="lines",
                line=dict(color=color, width=2, dash="dash"), yaxis="y2",
            ))
        
        fig.update_layout(
            title="Melting Curve Simulation",
            xaxis=dict(title="Temperature (°C)"),
            yaxis=dict(title="Duplex Fraction θ(T)", side="left"),
            yaxis2=dict(title="Melting Peak −dθ/dT", side="right", overlaying="y1"),
            hovermode="x unified", height=500, width=None,
        )
        
        st.plotly_chart(fig, use_container_width=True)
        st.dataframe(results_data, use_container_width=True)
        
        col_csv, col_png = st.columns(2)
        with col_csv:
            csv_buf = StringIO()
            writer = csv.DictWriter(csv_buf, fieldnames=results_data[0].keys())
            writer.writeheader()
            writer.writerows(results_data)
            st.download_button(
                label="📥 Download results (CSV)", data=csv_buf.getvalue(),
                file_name="primer_results.csv", mime="text/csv",
            )
        with col_png:
            st.info("💡 Use browser's right-click → Save image to download the plot")
        
        if len(sequences) == 1:
            st.subheader("📈 Detailed Analysis")
            r = simulate(sequences[0], na_mM=na_mm, oligo_nM=oligo_nm)
            peak_i = r["dtheta"].index(max(r["dtheta"]))
            peak_t = r["temps"][peak_i]
            
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("Tm (analytic)", f"{r['tm']:.2f}°C")
            col2.metric("Peak (−dθ/dT)", f"{peak_t:.2f}°C")
            col3.metric("Difference", f"{abs(peak_t - r['tm']):.2f}°C")
            col4.metric("Sequence length", len(sequences[0]))
            
            st.markdown("**Validation:** Peak of −dθ/dT should land within ~1°C of analytic Tm. "
                       "Cross-check with [IDT OligoAnalyzer](https://www.idtdna.com/calc/analyzer/) for external validation.")
else:
    st.info("👆 Enter a primer sequence above to get started")
