"""tools/R5_bio_rows.py — righe di BIO_REFERENCES per R5 dalle tabelle finali (results/R5). Uso: .venv/bin/python tools/R5_bio_rows.py <primo id>"""
import sys, pandas as pd, numpy as np
start = int(sys.argv[1]); T = "results/R5"
P = pd.read_csv(f"{T}/R5_predictions_edges.csv"); K = pd.read_csv(f"{T}/R5_kernel_descriptive.csv").set_index("key")
H = pd.read_csv(f"{T}/R5_PH2_inner_summary.csv"); N = pd.read_csv(f"{T}/R5_predictions_nc.csv")
DS = {"A1": "Visium_HD_Mouse_Small_Intestine", "A23": "Visium_HD_Human_Colon_Cancer", "A4": "Visium_HD_Human_Lymph_Node_FFPE", "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}
ARC = {"A1": "A1", "A23": "A2, A3", "A4": "A4", "A5": "A5", "A6": "A6"}
f = lambda x, d=3: f"{x:.{d}f}"
rows = []
def add(arc, q, v, iv, unit, src, note): rows.append((arc, q, v, iv, unit, src, note))
for k in ["A1", "A23", "A4", "A5", "A6"]:
    p = P[P.key == k].set_index("prediction"); kd = K.loc[k]; nt = int(p.loc["P-1", "n_tiles"])
    src = f"tools/R5_analysis.py@R5 (results/R5/R5_predictions_edges.csv, R5_kernel_descriptive.csv) + {DS[k]}, {nt} finestre di bordo tessuto/vetro"
    add(ARC[k], "larghezza della parte gaussiana del nucleo di diffusione ai bordi del tessuto (K2 = (1 - w) gaussiana(σ) + w esponenziale(ℓ)), σ", f(p.loc["P-1", "estimate"], 2),
        f"[{f(p.loc['P-1','lo'],2)}, {f(p.loc['P-1','hi'],2)}]; P-1 {p.loc['P-1','verdict']} (atteso [1, 4])", "µm", src, "mediana delle finestre, bootstrap per blocchi; dipende dal bordo H&E (maschera alternativa +0.6...+2.1 µm)")
    add(ARC[k], "massa e lunghezza della coda esponenziale del nucleo di diffusione (w, ℓ)", f"{f(p.loc['P-2','estimate'])} / {f(kd['ell'],1)}",
        f"w [{f(p.loc['P-2','lo'])}, {f(p.loc['P-2','hi'])}]; P-2 {p.loc['P-2','verdict']}", "adimensionale / µm", src, "ℓ mediana delle finestre")
    add(ARC[k], "rapporto σ_riga / σ_colonna (K3)", f(p.loc["P-3", "estimate"]), f"[{f(p.loc['P-3','lo'])}, {f(p.loc['P-3','hi'])}]; P-3 {p.loc['P-3','verdict']}", "adimensionale", src, "isotropia rispetto alla griglia")
    add(ARC[k], "scarto di registrazione immagine H&E → griglia (righe, colonne; modulo del vettore delle mediane)", f"{f(p.loc['P-4','dr_median'],2)}, {f(p.loc['P-4','dc_median'],2)}; |δ| {f(p.loc['P-4','estimate'],2)}",
        f"|δ| [{f(p.loc['P-4','lo'],2)}, {f(p.loc['P-4','hi'],2)}]; P-4 {p.loc['P-4','verdict']}; scarto interquartile fra finestre {f(kd['delta_spread_iqr_dr'],1)} / {f(kd['delta_spread_iqr_dc'],1)}", "µm", src,
        "relativo alla registrazione di Space Ranger (tissue_positions.parquet)")
    add(ARC[k], "spostamento medio dei trascritti con K2; frazione oltre 2, 4, 8, 16, 32 µm", f(kd["mean_displacement_um"], 1),
        " / ".join(f(kd[c]) for c in ("frac_beyond_2um", "frac_beyond_4um", "frac_beyond_8um", "frac_beyond_16um", "frac_beyond_32um")), "µm; frazione", src, "D-2")
    add(ARC[k], "fondo sul vetro lontano (UMI per bin a 2 µm)", f(kd["b_far"]), f"regola: {kd['b_far_rule']}; mediana del fondo b nelle finestre {f(kd['b_tiles_median'])}", "UMI per bin", src, "D-3")
for a in ["A1", "A2", "A3", "A4", "A5", "A6"]:
    h = H[H.tag.str.startswith(a + "_")]; ok = h[h.sigma_in_star.notna()]
    v = f"{f(ok.sigma_in_star.min(),2)}-{f(ok.sigma_in_star.max(),2)}" if len(ok) else "non stimabile"
    iv = f"soluzione in {len(ok)}/5 ROI; r prima della diffusione {f(ok.r_at_star.min(),2)}-{f(ok.r_at_star.max(),2)}" if len(ok) else "nessuna soluzione in 5/5 ROI"
    add(a, "sfocatura gaussiana dentro il tessuto dal contrasto nucleo/corona (σ_in, post hoc PH-2) e rapporto nucleo/corona prima della diffusione a σ_in", v, iv, "µm; adimensionale",
        "tools/R5_ph2_inner.py@R5 (results/R5/R5_PH2_inner_summary.csv) + 5 ROI di R2, nuclei di R4", "post hoc; controprova sintetica PH-2s (σ e r ritrovati entro 0.35 µm e 2 %); da confermare")
n5 = N[(N.archetype == "A5") & (N.prediction == "P-5") & N.variant.str.startswith("as written")].iloc[0]
add("A5", "rapporto nucleo/corona prima della diffusione con il nucleo dei bordi (P-5, come scritta)", f(n5.estimate, 2), f"[{f(n5.lo,2)}, {f(n5.hi,2)}]; CP-1 confermata (1.026)", "adimensionale",
    "tools/R5_nc.py, R5_analysis.py@R5 (results/R5/R5_predictions_nc.csv) + Visium_HD_6p5mm_Mouse_Brain", "fuori dal cervello non calcolabile come scritto (nucleo dei bordi troppo largo per il contrasto interno)")
for i, r in enumerate(rows):
    print(f"| B-{start + i} | " + " | ".join(str(x) for x in r[:6]) + " | R5 | " + r[6] + " |")
