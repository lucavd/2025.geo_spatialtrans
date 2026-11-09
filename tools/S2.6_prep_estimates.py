"""tools/S2.6_prep_estimates.py - preparazione di S2.6 (scritta alla chiusura di S2.5, su delega di Luca, DEC-071). NON e' un risultato di step:
stime da parametri e misure gia' registrati, per impostare la prova di prestazione. Legge: la tabella dei parametri del documento di progetto
dello Step 2 (sezione 9), results/S2.5/S2.5_verdicts_main.csv (P-8), results/S2.4/S2.4_R3_fits.csv, tissue_positions.parquet a 8 um (sola lettura).
Uso: .venv/bin/python tools/S2.6_prep_estimates.py -> results/S2.6/S2.6_prep_estimates.csv, S2.6_prep_overhead.csv"""
import os, glob, math, csv
import pandas as pd
DESIGN = "docs/superpowers/specs/2026-11-03-step2-visium-grid-design.md"
T = ["A1", "A2", "A3", "A4", "A5", "A6"]
SEC = {"A1": "A1", "A2": "A23", "A3": "A23", "A4": "A4", "A5": "A5", "A6": "A6"}
DS = {"A1": "Visium_HD_Mouse_Small_Intestine", "A23": "Visium_HD_Human_Colon_Cancer", "A4": "Visium_HD_Human_Lymph_Node_FFPE",
      "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}
lines = open(DESIGN, encoding="utf-8").read().split("\n")
def row(prefix):
    l = [x for x in lines if x.startswith(prefix)]; assert len(l) == 1, (prefix, len(l))
    return [float(x.strip().replace("\u2020", "").replace(" ", "")) for x in l[0].split("|")[2:8]]
w, ell = row("| coda esponenziale: massa w |"), row("| coda esponenziale: lunghezza \u2113 (\u00b5m) |")
w2, ell2 = row("| seconda coda: massa w\u2082"), row("| seconda coda: lunghezza \u2113\u2082 (\u00b5m) |")
beyond = lambda d, l: (1 + d / l) * math.exp(-d / l)   # raggio Gamma(2, l): massa oltre la distanza d
V = pd.read_csv("results/S2.5/S2.5_verdicts_main.csv"); P8 = V[V.id == "P-8"]
F = pd.read_csv("results/S2.4/S2.4_R3_fits.csv").set_index("key")
area = {}
for k, ds in DS.items():
    f = glob.glob(f"/mnt/micron/geo_spatialtrans/data_real/{ds}/**/square_008um/spatial/tissue_positions.parquet", recursive=True)
    t = pd.read_parquet(f[0], columns=["in_tissue"]); area[k] = (int(t.in_tissue.sum()) * 64e-6, len(t) * 64e-6)
os.makedirs("results/S2.6", exist_ok=True); out = []
for i, a in enumerate(T):
    s = float(P8[(P8.archetype == a) & (P8.metric == "elapsed_sim_s")].sim_est.iloc[0]); mb = float(P8[(P8.archetype == a) & (P8.metric == "peak_mb")].sim_est.iloc[0])
    tm, am = area[SEC[a]]
    out.append(dict(archetype=a, section=SEC[a], tissue_mm2_section=round(tm, 2), array_mm2=round(am, 2), w=w[i], ell_um=ell[i], w2=w2[i], ell2_um=ell2[i],
                    frac_all_beyond_100um_tail1=w[i] * beyond(100, ell[i]), frac_all_beyond_280um_tail1=w[i] * beyond(280, ell[i]),
                    frac_all_beyond_280um_tail2=w2[i] * beyond(280, ell2[i]), frac_all_beyond_500um_tail2=w2[i] * beyond(500, ell2[i]),
                    frac_all_beyond_1000um_tail2=w2[i] * beyond(1000, ell2[i]), share_of_tail2_beyond_280um=beyond(280, ell2[i]),
                    tail2_mass_in_capture_S24=float(F.loc[SEC[a], "tail_mass_in_capture"]),
                    s_per_mm2_area_S25_main=s, peak_mb_per_mm2_area_S25_main=mb, est_process_hours_section_no_halo=tm * s / 3600))
pd.DataFrame(out).to_csv("results/S2.6/S2.6_prep_estimates.csv", index=False)
ov = [dict(tile_um=L, halo_um=h, simulated_over_useful_area=((L + 2 * h) / L) ** 2) for L in (500, 1000, 2000) for h in (0, 100, 280)]
pd.DataFrame(ov).to_csv("results/S2.6/S2.6_prep_overhead.csv", index=False)
print(pd.DataFrame(out).round(4).to_string()); print(pd.DataFrame(ov).round(2).to_string())
