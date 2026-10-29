"""tools/S1.6_explore.py - S1.6: esplorazione prima della pre-registrazione (dichiarata in docs/steps/S1.6.md). Nessuna simulazione nuova.
E0: tools/S1.6_lib.py riproduce i verdetti di S1.5e (results/S1.5e/S1.5e_verdicts.csv) riga per riga.
E1 (BL-163): 90 percentile e mediana dell'area del nucleo per tipo principale, per ROI, su due insiemi di nuclei veri: celle interne (riferimento di
    S1.5d, tools/S1.5d_size_refs.py; e' anche l'insieme di P-8a e P-8b, tools/S1.5b_real_refs.py, dove `matched` e' sempre vero) e l'insieme su cui e'
    stimato il modello di taglia del preset: uguale nei quattro tessuti con Cellpose, nel linfonodo e nel cuore le celle interne con un nucleo di forma
    StarDist della maschera CO di R3b ([DEC-006], tools/S1.5d_export_r3b_cells.R: src == "sd" nella mappa della maschera), area > 0.
E3: errore standard tipico della differenza appaiata fra due bracci con la stessa legge: meta' dei semi (1-5) contro l'altra meta' (6-10) dello
    stesso braccio di S1.5e, riportato a 10 semi per lato (diviso per sqrt(2)), contro un quarto della larghezza dell'intervallo vero.
E1c: verdetti dei bracci U e Q di S1.5e con i riferimenti per tipo di tools/S1.6_real_refs.py (quali cambiano, prima di ogni simulazione nuova).
E2 (BL-173) e' in tools/S1.6_explore_parts.py, E4 (BL-161) in tools/S1.6_explore_bound.R, E5 (BL-169) in tools/S1.6_explore_refit.py.
Uscite: results/S1.6/S1.6_explore_q90_roi.csv, S1.6_explore_q90.csv, S1.6_explore_se.csv; stampa i riepiloghi."""
import importlib.util, os, sys, numpy as np, pandas as pd
def mod(n, p):
    s = importlib.util.spec_from_file_location(n, p); m = importlib.util.module_from_spec(s); s.loader.exec_module(m); return m
L = mod("s16lib", "tools/S1.6_lib.py")
OUT = "results/S1.6"; os.makedirs(OUT, exist_ok=True)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
# ---- E0 ----
R = L.load_refs(); B = L.load_bench("results/S1.5e")
old = pd.read_csv("results/S1.5e/S1.5e_verdicts.csv")
new = pd.concat([L.verdicts(L.roi_values(B, R, a), a) for a in old.arm.unique()], ignore_index=True)
k = ["id", "archetype", "unit", "arm"]
m = old.merge(new, on=k, suffixes=("_o", "_n"), how="outer", indicator=True)
both = m[m._merge == "both"]
dv = (both.verdict_o != both.verdict_n).sum(); de = np.nanmax(np.abs(both.est_o - both.est_n)); ds = np.nanmax(np.abs(both.se_o - both.se_n))
print(f"E0: righe S1.5e {len(old)}, lib {len(new)}, comuni {len(both)}, solo una parte {int((m._merge != 'both').sum())}, verdetti diversi {dv}, max |est| {de:.3g}, max |se| {ds:.3g}")
# ---- E1 ----
cal = mod("cal", "tools/S1.5b_calibration.py"); R3B = "/mnt/micron/geo_spatialtrans/R3b"
rows, fm, SETS = [], [], {}
for A in cal.ARCH:
    _, d = cal.load(A); e = d[d.interior.astype(bool)]
    if A in ("A4", "A6"):
        sd = []
        for r in sorted(d.roi_id.unique()):
            g = pd.read_parquet(f"{R3B}/R3/real/{A}_{r}_consensus_gen.parquet"); mp = pd.read_parquet(f"{R3B}/masks/{A}_{r}_consensus_map.parquet")
            g["sd"] = g.label.map(dict(zip(mp.label, mp.src))) == "sd"; g["roi_id"] = r; sd.append(g[["roi_id", "idx", "x", "y", "sd"]])
        sd = pd.concat(sd); e = e.merge(sd, on=["roi_id", "idx"], how="left", suffixes=("", "_g"), validate="one_to_one")
        assert e.sd.notna().all() and np.allclose(e.x, e.x_g) and np.allclose(e.y, e.y_g), A
        ms = e[e.sd.astype(bool) & (e.a_n > 0)]
    else:
        ms = e[e.matched.astype(bool) & (e.a_n > 0)]
    sets = {"interior": e, "model_set": ms}; SETS[(A, "interior")] = e; SETS[(A, "model_set")] = ms
    fm.append(dict(archetype=A, n_interior=len(e), n_model_set=len(sets["model_set"]), frac=len(sets["model_set"]) / len(e)))
    for nm, x in sets.items():
        for t in L.MAIN[A]:
            for r, gr in x[x.type == t].groupby("roi_id"):
                rows.append(dict(archetype=A, set=nm, cell_type=t, roi_id=r, n=len(gr), q90_an=np.quantile(gr.a_n, 0.9), median_an=np.median(gr.a_n)))
Q = pd.DataFrame(rows); Q.to_csv(f"{OUT}/S1.6_explore_q90_roi.csv", index=False)
S = Q.groupby(["archetype", "cell_type", "set"]).agg(q90_min=("q90_an", "min"), q90_max=("q90_an", "max"), n_roi=("q90_an", "size")).reset_index()
S.to_csv(f"{OUT}/S1.6_explore_q90.csv", index=False)
print("E1 frazione dell'insieme del modello fra le celle interne:"); print(pd.DataFrame(fm).round(4).to_string(index=False))
ref = pd.read_csv("results/S1.5d/S1.5d_real_q90.csv")
chk = S[S.set == "interior"].merge(ref, on=["archetype", "cell_type"], suffixes=("", "_s15d"))
print("E1 controllo: celle interne = S1.5d_real_q90.csv, max scarto", float(np.max(np.abs(chk.q90_min - chk.q90_min_s15d)) + np.max(np.abs(chk.q90_max - chk.q90_max_s15d))))
S2 = Q.groupby(["archetype", "cell_type", "set"]).agg(med_min=("median_an", "min"), med_max=("median_an", "max")).reset_index(); S = S.merge(S2)
S.to_csv(f"{OUT}/S1.6_explore_q90.csv", index=False)
w = S.pivot_table(index=["archetype", "cell_type"], columns="set", values=["q90_min", "q90_max", "med_min", "med_max"]).round(3)
print(w.to_string())

# ---- E1b: riferimenti per ROI del tessuto (P-3, P-8a: results/S1.5/S1.5_real_refs_roi.csv) e per tipo (P-4, P-8b: S1.5b_real_refs_roi_type.csv)
RR = pd.read_csv("results/S1.5/S1.5_real_refs_roi.csv"); RT = pd.read_csv("results/S1.5b/S1.5b_real_refs_roi_type.csv"); RT = RT[RT.level == "type"]
tis = []
for (A, nm), x in SETS.items():
    for r, gr in x.groupby("roi_id"):
        tis.append(dict(archetype=A, set=nm, roi_id=r, med_an=np.median(gr.a_n), med_nc=np.median(gr.a_n / gr.a_t), n=len(gr)))
TI = pd.DataFrame(tis).merge(RR[["archetype", "roi_id", "median_area_nuc", "median_nc"]], on=["archetype", "roi_id"])
TI["d_an"] = TI.med_an - TI.median_area_nuc; TI["d_nc"] = TI.med_nc - TI.median_nc
print("E1b tessuto: scarto massimo dal riferimento di S1.5 per insieme"); print(TI.groupby(["archetype", "set"])[["d_an", "d_nc"]].agg(lambda v: np.max(np.abs(v))).round(4).to_string())
QT = Q.merge(RT[["archetype", "roi_id", "unit", "median_area_nuc"]].rename(columns={"unit": "cell_type"}), on=["archetype", "roi_id", "cell_type"])
QT["d"] = QT.median_an - QT.median_area_nuc
print("E1b tipo: scarto massimo dal riferimento di S1.5b per insieme"); print(QT.groupby(["archetype", "set"]).d.agg(lambda v: np.max(np.abs(v))).round(4).to_string())
TI.to_csv(f"{OUT}/S1.6_explore_tissue_refs.csv", index=False)

# ---- E1c: verdetti dei bracci U e Q di S1.5e con i riferimenti per tipo nuovi (tools/S1.6_real_refs.py) ----
R2 = L.load_refs("results/S1.6/S1.6_real_q90.csv", "results/S1.6/S1.6_real_refs_roi_type.csv")
cmp = []
for arm in ["U", "Q"]:
    a = L.verdicts(L.roi_values(B, R, arm), arm); b = L.verdicts(L.roi_values(B, R2, arm), arm)
    m = a.merge(b, on=["id", "archetype", "unit", "arm"], suffixes=("_old", "_new")); cmp.append(m)
CM = pd.concat(cmp); CM.to_csv(f"{OUT}/S1.6_explore_refs_verdicts.csv", index=False)
ch = CM[(CM.verdict_old != CM.verdict_new) | (CM.lo_old != CM.lo_new) | (CM.hi_old != CM.hi_new)]
print("E1c: righe con riferimento cambiato:", len(ch), "; verdetti cambiati:", int((CM.verdict_old != CM.verdict_new).sum()))
print(ch[["arm", "id", "archetype", "unit", "est_old", "se_old", "lo_old", "hi_old", "lo_new", "hi_new", "verdict_old", "verdict_new"]].round(4).to_string(index=False))
# ---- E3 ----
se_rows = []
for arm in ["U", "Q"]:
    ra = L.roi_values(B, R, arm, seeds=range(1, 6)); rb = L.roi_values(B, R, arm, seeds=range(6, 11))
    sk = rb.level == "seed"; rb.loc[sk, "key"] = rb.loc[sk, "key"].map(lambda s: f"s{int(s.lstrip('s')) - 5:02d}")   # P-1: seme s+5 appaiato a s
    p = L.paired(ra, rb); p["arm"] = arm; p["se10"] = p.se / np.sqrt(2); p["ratio_2se_margin"] = 2 * p.se10 / p.margin
    se_rows.append(p)
P = pd.concat(se_rows, ignore_index=True); P.to_csv(f"{OUT}/S1.6_explore_se.csv", index=False)
print("E3: 2 ES tipici (10 semi) / margine (un quarto della larghezza), per previsione, braccio Q:")
print(P[P.arm == "Q"].groupby("id").ratio_2se_margin.describe()[["count", "50%", "max"]].round(3).to_string())
print("E3: righe con 2 ES > margine:", int((P.ratio_2se_margin > 1).sum()), "su", len(P))
print(P[(P.arm == "Q") & (P.ratio_2se_margin > 1)][["arm", "id", "archetype", "unit", "diff", "se10", "margin", "ratio_2se_margin"]].round(4).to_string(index=False))
