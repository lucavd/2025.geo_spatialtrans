#!/usr/bin/env python3
"""tools/R3b_report_parts.py — estrae dal file dello step (docs/steps/R3b.md) le sezioni incluse nel report (il testo esiste una
volta sola) e costruisce la tabella «Da dove vengono i numeri del riassunto» rileggendo ogni valore dai CSV di results/R3b.
Uso (root del repo): .venv/bin/python -B tools/R3b_report_parts.py"""
import pathlib, re
import numpy as np, pandas as pd
R = pathlib.Path("results/R3b"); out = pathlib.Path("reports")
t = pathlib.Path("docs/steps/R3b.md").read_text(encoding="utf-8").split("\n---\n", 1)[1]
def sec(start, stop_re):
    i = t.index(start); m = re.search(stop_re, t[i + len(start):]); j = i + len(start) + (m.start() if m else len(t))
    return t[i:j].strip() + "\n"
(out / "_R3b_inbreve.md").write_text(sec("\n# In breve\n", r"\n# "), encoding="utf-8")
(out / "_R3b_rinviati.md").write_text(sec("\n# Dubbi rinviati\n", r"\n# ").replace("# Dubbi rinviati\n", "", 1), encoding="utf-8")
(out / "_R3b_verifiche.md").write_text(sec("\n# Verifiche di Luca\n", r"\n# ").replace("# Verifiche di Luca\n", "", 1), encoding="utf-8")
rev = sec("\n## Revisione avversariale (cancello G-5)\n", r"\n#{1,2} ").replace("## Revisione avversariale (cancello G-5)\n", "", 1)
post = sec("\n## Dopo il risultato e la revisione avversariale", r"\n#{1,2} ")
(out / "_R3b_review.md").write_text(rev + "\n" + post.replace("## Dopo il risultato", "## Note dopo il risultato", 1), encoding="utf-8")
# ---- fonti dei numeri del riassunto
s = pd.read_csv(R / "R3b_p1p2_summary.csv").set_index(["archetype", "mask"]); p = pd.read_csv(R / "R3b_predictions.csv")
pv = lambda i, a, m: float(p[(p.id == i) & (p.archetype == a) & (p["mask"] == m)].est.iloc[0])
sp = pd.read_csv(R / "R3b_posthoc_srprec.csv").groupby(["archetype", "src"])[["n_counted", "n_matched_to_point"]].sum()
ag = pd.read_csv(R / "R3b_posthoc_agreement.csv").set_index(["archetype", "mask"])
pt = pd.read_csv(R / "R3b_posthoc_partners.csv"); pt4 = pt[(pt.archetype == "A4") & pt.contact]
pa = pd.read_csv(R / "R3b_posthoc_area.csv").set_index("archetype"); rf = pd.read_csv(R / "R3b_refs.csv")
rv = lambda a, m, c: float(rf[(rf.archetype == a) & (rf.metric == m)][c].iloc[0])
w = pd.read_csv(R / "R3b_posthoc_sdarea_windows.csv"); sda = w[w.archetype == "A4"].groupby("seed").s1.mean().mean()
mp = pd.concat([pd.read_parquet(f) for f in sorted(pathlib.Path("/mnt/micron/geo_spatialtrans/R3b/masks").glob("A4_*_map.parquet"))])
ms = mp[mp.src == "sd"]; ratio = float((ms.n_px_sr / ms.n_px_sd).median())
terr = rf[rf.metric.isin(["area", "area_q25", "area_q75", "area_roi_min", "area_roi_max", "eq_r", "cv_loc", "ecc_T"])]
vd = pd.read_csv(R / "R3b_posthoc_valid.csv").groupby("mask")[["n_labels", "n_out_raw", "n_out_after_cp2", "n_gen_out"]].sum()
rows = [
 ("nuclei della maschera scelta nelle finestre di linfonodo / punti di Luca", "R3b_p1p2_summary.csv A4 consensus n_obj, n_manual", f"{s.loc[('A4','consensus')].n_obj} / {s.loc[('A4','consensus')].n_manual}"),
 ("StarDist da solo nel linfonodo", "R3b_p1p2_summary.csv A4 stardist_he n_obj; p1_ratio_mean", f"{s.loc[('A4','stardist_he')].n_obj}; {s.loc[('A4','stardist_he')].p1_ratio_mean:.3f}"),
 ("cuore: maschera scelta, StarDist, Luca", "R3b_p1p2_summary.csv A6 n_obj, n_manual", f"{s.loc[('A6','consensus')].n_obj}, {s.loc[('A6','stardist_he')].n_obj}, {s.loc[('A6','consensus')].n_manual}"),
 ("poligoni Space Ranger rimasti nelle finestre di linfonodo, di cui con un punto di Luca", "R3b_posthoc_srprec.csv A4 sr", f"{sp.loc[('A4','sr')].n_counted}, {sp.loc[('A4','sr')].n_matched_to_point}"),
 ("idem nel cuore", "R3b_posthoc_srprec.csv A6 sr", f"{sp.loc[('A6','sr')].n_counted}, {sp.loc[('A6','sr')].n_matched_to_point}"),
 ("contatti simulati nel linfonodo: punti di Luca, maschera scelta, StarDist", "R3b_predictions.csv CP-1, P-3 A4", f"{pv('CP-1','A4','manuale+CO'):.3f}, {pv('P-3','A4','CO'):.3f}, {pv('P-3','A4','SD'):.3f}"),
 ("contatti simulati nel cuore (maschera scelta)", "R3b_predictions.csv P-3 A6 CO", f"{pv('P-3','A6','CO'):.4f}"),
 ("valore manuale di Luca e margine (linfonodo, cuore)", "results/R2c/R2c_posthoc_ci.csv, due stadi", "; ".join(f"{r.est:.3f} [{r.lo:.2f}, {r.hi:.3f}]" for r in pd.read_csv("results/R2c/R2c_posthoc_ci.csv").query("ci == 'due stadi' and archetype in ['A4','A6']").itertuples())),
 ("frazione di nuclei con forma Space Ranger (celle interne)", "R3b_posthoc_area.csv frac_sr_ev", f"{pa.loc['A4'].frac_sr_ev:.3f}, {pa.loc['A6'].frac_sr_ev:.3f}"),
 ("poligono Space Ranger / StarDist sullo stesso nucleo, linfonodo (mediana)", "masks/A4_*_consensus_map.parquet, src = sd", f"{ratio:.3f}"),
 ("variabilita' della taglia nel linfonodo: tutti, solo forma StarDist, oggi", "R3b_posthoc_area.csv sigma_all, sigma_sd; R3b_refs.csv sigma old_recomputed", f"{pa.loc['A4'].sigma_all:.3f}, {pa.loc['A4'].sigma_sd:.3f}, {rv('A4','sigma','old_recomputed'):.3f}"),
 ("contatti simulati nel linfonodo con la taglia dai soli nuclei di forma StarDist", "R3b_posthoc_sdarea_windows.csv A4", f"{sda:.3f}"),
 ("nuclei focali del linfonodo: contatti nella maschera, di Luca, separati per Luca", "R3b_posthoc_agreement.csv A4 consensus", f"{ag.loc[('A4','consensus')].mask_touch}, {ag.loc[('A4','consensus')].manual_touch}, {ag.loc[('A4','consensus')].mask_only} su {ag.loc[('A4','consensus')].n}"),
 ("contatti che passano da un poligono Space Ranger", "R3b_posthoc_partners.csv A4, contact = True, any_sr_partner", f"{pt4.any_sr_partner.mean():.3f}"),
 ("cambiamento massimo delle misure dei territori", "R3b_refs.csv rel_change, metriche dei territori", f"{terr.rel_change.abs().max():.4f}"),
 ("rapporto nucleo/cellula mediano nel linfonodo: oggi, nuovo", "R3b_refs.csv A4 nc", f"{rv('A4','nc','old_recomputed'):.3f}, {rv('A4','nc','consensus'):.3f}"),
 ("tetto del nucleo nel cuore: oggi, nuovo (tutti / solo forma StarDist)", "R3b_refs.csv A6 phi; R3b_posthoc_area.csv phi_sd", f"{rv('A6','phi','old_recomputed'):.3f}, {rv('A6','phi','consensus'):.3f} / {pa.loc['A6'].phi_sd:.3f}"),
 ("nuclei con centro fuori dal tessuto valido: maschera scelta, Space Ranger, StarDist; dopo il filtro; generatori", "R3b_posthoc_valid.csv", f"{vd.loc['consensus'].n_out_raw}, {vd.loc['spaceranger'].n_out_raw}, {vd.loc['stardist_he'].n_out_raw}; {int(vd.n_out_after_cp2.sum())}; {int(vd.n_gen_out.sum())}"),
 ("nuclei nella maschera scelta (10 ROI)", "R3b_posthoc_valid.csv consensus n_labels", f"{vd.loc['consensus'].n_labels}"),
]
f = "## Da dove vengono i numeri del riassunto\n\n| frase del riassunto | fonte (file, riga) | valore letto |\n|---|---|---|\n" + "".join(f"| {a} | {b} | {c} |\n" for a, b, c in rows)
(out / "_R3b_fonti.md").write_text(f, encoding="utf-8"); print(f)
