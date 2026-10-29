"""tools/S1.6_real_refs.py - S1.6 (BL-163): riferimenti per tipo della taglia del nucleo sullo stesso insieme di nuclei del modello di taglia.
Nel linfonodo (A4) e nel cuore (A6) il modello di taglia dei preset e i riferimenti del tessuto (P-3, P-8a: results/S1.5/S1.5_real_refs_roi.csv) sono
stimati sulle celle interne con un nucleo di forma StarDist della maschera CO di R3b ([DEC-006]); i riferimenti per tipo (P-4, P-8b:
results/S1.5b/S1.5b_real_refs_roi_type.csv; P-14: results/S1.5d/S1.5d_real_q90.csv) su tutte le celle interne della maschera CO (esplorazione E1 di
S1.6). Qui le righe per tipo di A4 e A6 si ricalcolano sull'insieme del modello: celle interne, src == "sd" nella mappa della maschera CO, area > 0;
stesse definizioni di tools/S1.5b_real_refs.py e tools/S1.5d_size_refs.py. Gli altri quattro tessuti restano identici (stesso insieme gia' oggi).
Calcolato prima di ogni simulazione di S1.6. Uscite: results/S1.6/S1.6_real_refs_roi_type.csv (formato di S1.5b), results/S1.6/S1.6_real_q90_roi.csv,
results/S1.6/S1.6_real_q90.csv (formato di S1.5d)."""
import importlib.util, os, numpy as np, pandas as pd
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
cal = mod("cal", "tools/S1.5b_calibration.py"); R3B = "/mnt/micron/geo_spatialtrans/R3b"; OUT = "results/S1.6"; os.makedirs(OUT, exist_ok=True)
MAIN = {"A1": ["epiteliale", "stromale e muscolo liscio", "immunitario"], "A2": ["epiteliale", "fibroblasti"],
        "A3": ["epiteliale", "fibroblasti", "endoteliale", "immunitario"], "A4": ["linfociti B", "linfociti T"],
        "A5": ["neuroni", "oligodendrociti", "astrociti"], "A6": ["cardiomiociti", "fibroblasti", "immunitari"]}
C = pd.read_csv("results/S1.5b/S1.5b_type_calibration.csv")
RT_old = pd.read_csv("results/S1.5b/S1.5b_real_refs_roi_type.csv"); Q_old = pd.read_csv("results/S1.5d/S1.5d_real_q90_roi.csv")
def sd_flag(A, d):
    parts = []
    for r in sorted(d.roi_id.unique()):
        g = pd.read_parquet(f"{R3B}/R3/real/{A}_{r}_consensus_gen.parquet"); mp = pd.read_parquet(f"{R3B}/masks/{A}_{r}_consensus_map.parquet")
        g["sd"] = g.label.map(dict(zip(mp.label, mp.src))) == "sd"; g["roi_id"] = r; parts.append(g[["roi_id", "idx", "x", "y", "sd"]])
    e = d.merge(pd.concat(parts), on=["roi_id", "idx"], how="left", suffixes=("", "_g"), validate="one_to_one")
    assert e.sd.notna().all() and np.allclose(e.x, e.x_g) and np.allclose(e.y, e.y_g), A
    return e
rows, qrows = [], []
for A in ["A4", "A6"]:
    _, d = cal.load(A); d = sd_flag(A, d)
    ev = d[d.interior & d.sd.astype(bool) & (d.a_n > 0)].copy()
    g = C[C.archetype == A].set_index("cell_type").size_group
    ev["size_group"] = ev.type.map(g); ev["nt"] = ev.a_n / ev.a_t
    for r, e in ev.groupby("roi_id"):
        units = [("type", t, x) for t, x in e.groupby("type")] + [("size_group", int(s), x) for s, x in e.groupby("size_group")]
        for k, v, x in units:
            rows.append(dict(archetype=A, roi_id=r, level=k, unit=str(v), n=len(x), median_area_nuc=x.a_n.median(), median_nt=x.nt.median(),
                             median_a_t=x.a_t.median(), frac_nt095=float((x.nt >= 0.95).mean())))
    for t in MAIN[A]:
        for r, gr in ev[ev.type == t].groupby("roi_id"):
            qrows.append(dict(archetype=A, cell_type=t, roi_id=r, n=len(gr), q90_an=np.quantile(gr.a_n, 0.9), q10_an=np.quantile(gr.a_n, 0.1)))
RT = pd.concat([RT_old[~RT_old.archetype.isin(["A4", "A6"])], pd.DataFrame(rows)], ignore_index=True)
Qr = pd.concat([Q_old[~Q_old.archetype.isin(["A4", "A6"])], pd.DataFrame(qrows)], ignore_index=True)
RT.to_csv(f"{OUT}/S1.6_real_refs_roi_type.csv", index=False); Qr.to_csv(f"{OUT}/S1.6_real_q90_roi.csv", index=False)
Qs = Qr.groupby(["archetype", "cell_type"]).agg(q90_min=("q90_an", "min"), q90_max=("q90_an", "max"), n_roi=("q90_an", "size")).reset_index()
Qs.to_csv(f"{OUT}/S1.6_real_q90.csv", index=False)
T = RT[RT.level == "type"]; T = T[T.apply(lambda z: z.unit in MAIN[z.archetype], axis=1)]
s = T.groupby(["archetype", "unit"]).agg(area_min=("median_area_nuc", "min"), area_max=("median_area_nuc", "max"), nt_min=("median_nt", "min"), nt_max=("median_nt", "max"), n_roi=("roi_id", "size")).reset_index()
o = RT_old[RT_old.level == "type"].groupby(["archetype", "unit"]).agg(area_min_old=("median_area_nuc", "min"), area_max_old=("median_area_nuc", "max"), nt_min_old=("median_nt", "min"), nt_max_old=("median_nt", "max")).reset_index()
s = s.merge(o, on=["archetype", "unit"]).merge(Qs.rename(columns={"cell_type": "unit"}), on=["archetype", "unit"])
pd.set_option("display.width", 250); print(s[s.archetype.isin(["A4", "A6"])].round(4).to_string(index=False))
a = RT_old[~RT_old.archetype.isin(["A4", "A6"])].reset_index(drop=True); b = pd.read_csv(f"{OUT}/S1.6_real_refs_roi_type.csv"); b = b[~b.archetype.isin(["A4", "A6"])].reset_index(drop=True)
a2 = pd.read_csv("results/S1.5b/S1.5b_real_refs_roi_type.csv"); a2 = a2[~a2.archetype.isin(["A4", "A6"])].reset_index(drop=True)
num = [k for k in a2.columns if a2[k].dtype.kind in "fi"]; txt = [k for k in a2.columns if k not in num]
same = len(a2) == len(b) and a2[txt].astype(str).equals(b[txt].astype(str)) and all(np.allclose(a2[k], b[k], rtol=0, atol=1e-12, equal_nan=True) for k in num)
print("righe fuori da A4 e A6 uguali a S1.5b entro 1e-12 (giro di lettura e scrittura del CSV):", same, len(a2), len(b))
q2 = pd.read_csv("results/S1.5d/S1.5d_real_q90.csv"); q3 = Qs[~Qs.archetype.isin(["A4", "A6"])].reset_index(drop=True); q2 = q2[~q2.archetype.isin(["A4", "A6"])].reset_index(drop=True)
print("90 percentile fuori da A4 e A6 identico a S1.5d:", bool(np.allclose(q2.q90_min, q3.q90_min, rtol=0, atol=1e-12) and np.allclose(q2.q90_max, q3.q90_max, rtol=0, atol=1e-12)), len(q2))
