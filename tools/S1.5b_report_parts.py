#!/usr/bin/env python3
"""tools/S1.5b_report_parts.py - estrae dal file dello step (docs/steps/S1.5b.md) le sezioni incluse nel report e costruisce la tabella
«Da dove vengono i numeri del riassunto» rileggendo ogni valore dai file dei risultati. Uso (root del repo): .venv/bin/python -B tools/S1.5b_report_parts.py"""
import pathlib, re, numpy as np, pandas as pd
R = pathlib.Path("results/S1.5b"); out = pathlib.Path("reports")
t = pathlib.Path("docs/steps/S1.5b.md").read_text(encoding="utf-8").split("\n---\n", 1)[1]
def sec(start, stop_re):
    i = t.index(start); m = re.search(stop_re, t[i + len(start):]); j = i + len(start) + (m.start() if m else len(t))
    return t[i:j].strip() + "\n"
(out / "_S1.5b_inbreve.md").write_text(sec("\n# In breve\n", r"\n# "), encoding="utf-8")
(out / "_S1.5b_rinviati.md").write_text(sec("\n# Dubbi rinviati\n", r"\n# ").replace("# Dubbi rinviati\n", "", 1), encoding="utf-8")
(out / "_S1.5b_verifiche.md").write_text(sec("\n# Verifiche di Luca\n", r"\n# ").replace("# Verifiche di Luca\n", "", 1), encoding="utf-8")
(out / "_S1.5b_dopo.md").write_text(sec("\n# Dopo il risultato\n", r"\n# In breve\n").replace("# Dopo il risultato\n", "", 1).replace("](../decisions/", "](../docs/decisions/"), encoding="utf-8")
V = pd.read_csv(R / "S1.5b_verdicts.csv"); W = pd.read_csv(R / "_v1_before_review/S1.5b_verdicts.csv"); AR = pd.read_csv(R / "S1.5b_arm_rule.csv").set_index("arm")
def v(i, a, u, arm="recentred"): return V[(V.id == i) & (V.archetype == a) & (V.unit == u) & (V.arm == arm)].iloc[0]
cal = pd.read_csv(R / "S1.5b_type_calibration.csv"); a5 = pd.read_csv(R / "S1.5b_a5_composition.csv"); ph = pd.read_csv(R / "S1.5b_posthoc_a5.csv")
sm = pd.read_csv(R / "S1.5b_posthoc_a5_small.csv"); nn = pd.read_csv(R / "S1.5b_desc_nn.csv")
nnm = nn.groupby(["archetype", "cell_type", "source"]).nn_median.mean()
r6 = pd.read_csv("results/R6/R6_refs_roi_type_ph1.csv"); r6 = r6[r6.archetype == "A5"]
fr = (r6[r6.cell_type == "neuroni"].set_index("roi_id").n / r6.groupby("roi_id").n.sum())
m = V.merge(W, on=["id", "archetype", "unit", "arm"], suffixes=("", "_v1")); nchg = int((m.verdict != m.verdict_v1).sum())
rc = pathlib.Path("/mnt/micron/geo_spatialtrans/S1.5b/review_code/rc_territory_by_type.csv")
mye = "file del revisore non trovato"
if rc.exists():
    d = pd.read_csv(rc); k = d[(d.astype(str).apply(lambda c: c.str.contains("mieloidi")).any(axis=1))]
    mye = k.to_string(index=False)[:300].replace("\n", " ; ")
nb = cal.set_index(["archetype", "cell_type"])
rows = [
 ("taglia del nucleo per tipo in 14 tipi su 17 (senza tipi in 9)", "S1.5b_arm_rule.csv, p2_confirmed", f"recentred {AR.loc['recentred','p2_confirmed']}, group {AR.loc['group','p2_confirmed']}, single {AR.loc['single','p2_confirmed']} su {AR.loc['recentred','p2_n']}"),
 ("ha deciso lo stroma (pareggio, poi P-1)", "S1.5b_arm_rule.csv; S1.5b_verdicts.csv P-1 A3 group", f"decided_by {AR.loc['recentred','decided_by']}; P-1 recentred {AR.loc['recentred','p1_confirmed']}, group {AR.loc['group','p1_confirmed']}; A3 group {v('P-1','A3','tessuto','group').est:.3f} ± 2·{v('P-1','A3','tessuto','group').se:.3f} contro {v('P-1','A3','tessuto','group').hi:.3f}"),
 ("cellule quasi senza citoplasma da 3,6 a 2,2 su 100; massimo reale 2,9", "S1.5b_verdicts.csv, P-5 A2", f"single {v('P-5','A2','tessuto','single').est:.4f}, recentred {v('P-5','A2','tessuto').est:.4f}, estremo {v('P-5','A2','tessuto').hi}"),
 ("astrociti 34 µm² contro 48,5-54", "S1.5b_verdicts.csv, P-2 A5 astrociti", f"{v('P-2','A5','astrociti').est:.2f}; {v('P-2','A5','astrociti').lo:.2f}-{v('P-2','A5','astrociti').hi:.2f}"),
 ("neuroni: media geometrica sotto la mediana", "step, PH-3 (calcolo sull'insieme di calibrazione)", "mediana 111.6, media geometrica 99.4, simulato " + f"{v('P-2','A5','neuroni').est:.1f}"),
 ("B 16,9 contro T 15,9 µm² (6 %)", "S1.5b_verdicts.csv, P-2 A4", f"B {v('P-2','A4','linfociti B').est:.2f}, T {v('P-2','A4','linfociti T').est:.2f}, rapporto {v('P-2','A4','linfociti B').est / v('P-2','A4','linfociti T').est:.3f}"),
 ("ultimo tipo seminato con territori più piccoli del 10 % (mieloidi del linfonodo)", "review_code/rc_territory_by_type.csv (revisore, seme 1, v1); test T12", mye),
 ("cambiati 6 esiti su 165", "S1.5b_verdicts.csv contro _v1_before_review/S1.5b_verdicts.csv", f"{nchg} su {len(m)}"),
 ("neuroni 3,5 volte la glia; aree vere 2,3-7", "S1.5b_verdicts.csv, P-3", f"{v('P-3','A5','neuroni/glia').est:.2f}; {v('P-3','A5','neuroni/glia').lo:.2f}-{v('P-3','A5','neuroni/glia').hi:.2f}"),
 ("oggetti di 4-8 µm², soprattutto in un'area con sostanza bianca; dispersione della glia giù di quasi un quinto", "S1.5b_posthoc_a5_small.csv", f"frazione < 8 µm² in r1 {sm.frac_small_r1.iloc[0]:.3f}; sigma {sm.sigma_glia.iloc[0]:.3f} -> {sm.sigma_glia.iloc[1]:.3f} ({1 - sm.sigma_glia.iloc[1] / sm.sigma_glia.iloc[0]:.0%})"),
 ("rapporto nucleo/territorio sbaglia in 11 tipi su 17", "S1.5b_verdicts.csv, P-7 recentred", f"smentite {int(((V.id=='P-7')&(V.arm=='recentred')&(V.verdict=='smentita')).sum())} su {int(((V.id=='P-7')&(V.arm=='recentred')).sum())}"),
 ("immunitarie del cuore 8,8 µm dal vicino contro 14,5 simulati", "S1.5b_desc_nn.csv (media dei ROI)", f"reale {nnm[('A6','immunitari','reale')]:.2f}, simulato {nnm[('A6','immunitari','recentred')]:.2f}"),
 ("quota di neuroni dal 65 al 70 %", "S1.5b_type_calibration.csv, A5 neuroni", f"{nb.loc[('A5','neuroni'),'fraction_r6']:.3f} -> {nb.loc[('A5','neuroni'),'fraction_tissue']:.3f}"),
 ("indovina 9 casi su 10 sui nuclei trovati", "S1.5b_a5_composition.csv, Luca", f"{a5[a5.rater=='Luca'].accuracy_matched.iloc[0]:.3f}"),
 ("66 % con la correzione; 61 % con le due valutazioni", "S1.5b_posthoc_a5.csv; (1-q)·R6 + q·quota corretta dei persi, q dei due valutatori insieme", f"{ph[(ph.rater=='Luca')&(ph.cell_type=='neuroni')].composition_if_corrected.iloc[0]:.3f}; " + (lambda q: f"{(1 - q) * nb.loc[('A5','neuroni'),'fraction_r6'] + q * ph[(ph.rater=='both')&(ph.cell_type=='neuroni')].corrected_share.iloc[0]:.3f} (q dei due valutatori = {q:.3f})")(a5.groupby('rater').n_missed.first().sum() / a5.groupby('rater').n_points.first().sum())),
 ("variabilità fra le aree vere dal 53 al 76 % di neuroni", "results/R6/R6_refs_roi_type_ph1.csv, A5", f"{fr.min():.3f}-{fr.max():.3f}"),
]
tab = "| frase del riassunto | fonte | valore letto |\n|---|---|---|\n" + "\n".join(f"| {a} | {b} | {c} |" for a, b, c in rows) + "\n"
(out / "_S1.5b_fonti.md").write_text(tab, encoding="utf-8"); print(tab)
