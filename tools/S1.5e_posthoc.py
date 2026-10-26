"""tools/S1.5e_posthoc.py - S1.5e, PH-1 (dichiarato dopo il banco, docs/steps/S1.5e.md): scomposizione della perdita della calibrazione di tutte le ROI
per bersaglio, nel cuore e nel linfonodo, fra l'ottimo di XA e il miglior punto con pochi vuoti (phi <= 0.05), con la perdita di ciascuna delle 2 repliche.
Ingressi: /mnt/micron/geo_spatialtrans/S1.5e/calpos_met/all_units.csv, results/S1.5e/S1.5e_het_real.csv, S1.5e_het_loss.csv, griglie.
Uscita: results/S1.5e/S1.5e_posthoc_loss.csv. Uso: .venv/bin/python -B tools/S1.5e_posthoc.py"""
import importlib.util, pandas as pd
s = importlib.util.spec_from_file_location("m", "tools/S1.5e_metrics.py"); m = importlib.util.module_from_spec(s); s.loader.exec_module(m)
RES = "results/S1.5e"
L = pd.read_csv(f"{RES}/S1.5e_het_loss.csv"); G = pd.concat([pd.read_csv(f"{RES}/S1.5e_het_grid.csv"), pd.read_csv(f"{RES}/S1.5e_het_grid_ext.csv")]).drop_duplicates("combo")
L = L.merge(G[["combo", "void_frac", "void_ell_um", "void_aspect"]], on="combo", how="left")
U = pd.read_csv("/mnt/micron/geo_spatialtrans/S1.5e/calpos_met/all_units.csv"); real = pd.read_csv(f"{RES}/S1.5e_het_real.csv"); rows = []
for A in ["A4", "A6"]:
    x = L[(L.archetype == A) & (L.fold == "pooled") & (L.family == "XA")].sort_values("loss")
    picks = {"ottimo XA": x.iloc[0].combo, "migliore con phi <= 0.05": x[x.void_frac <= 0.05].iloc[0].combo}
    T = [t for t in m.CAL_TARGETS_E if not (A == "A4" and t == "elong_wmed")]; ra = real[real.archetype == A]; var = ra[T].var(ddof=1)
    for lab, cb in picks.items():
        u = U[(U.archetype == A) & (U.combo == cb)]
        z = u.groupby("roi_id")[T].mean().reset_index().merge(ra[["roi_id"] + T], on="roi_id", suffixes=("", "_r"))
        o = dict(archetype=A, which=lab, combo=cb, loss=float(x[x.combo == cb].loss.iloc[0]))
        for rep in [1, 2]:
            zz = u[u.rep == rep].merge(ra[["roi_id"] + T], on="roi_id", suffixes=("", "_r")); o[f"loss_rep{rep}"] = float(sum(((zz[t] - zz[t + "_r"]) ** 2 / var[t]).sum() for t in T))
        for t in T: o[f"contrib_{t}"] = float(((z[t] - z[t + "_r"]) ** 2 / var[t]).sum())
        rows.append(o)
D = pd.DataFrame(rows); D.to_csv(f"{RES}/S1.5e_posthoc_loss.csv", index=False); pd.set_option("display.width", 250); print(D.round(2).T.to_string())
