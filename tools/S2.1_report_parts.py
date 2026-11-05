#!/usr/bin/env python3
"""tools/S2.1_report_parts.py - estrae dal file dello step (docs/steps/S2.1.md) le sezioni incluse nel report e costruisce la tabella
«Da dove vengono i numeri del riassunto»: ogni numero di «In breve» e' riletto dai file dei risultati e dai log e confrontato con il valore
scritto (si ferma se uno non coincide). Uso (root del repo): .venv/bin/python -B tools/S2.1_report_parts.py
Uscite: reports/_S2.1_*.md, results/S2.1/S2.1_verification_log.csv"""
import pathlib, re, pandas as pd
R = pathlib.Path("results/S2.1"); out = pathlib.Path("reports"); M = pathlib.Path("/mnt/micron/geo_spatialtrans/S2.1")
full = pathlib.Path("docs/steps/S2.1.md").read_text(encoding="utf-8"); t = "\n" + full.split("\n---\n", 1)[1]
def sec(start, stop_re="\n# "):
    i = t.index(start); m = re.search(stop_re, t[i + len(start):]); j = i + len(start) + (m.start() if m else len(t)); return t[i:j].strip() + "\n"
fix = lambda s: s.replace("](/decisions/", "](../docs/decisions/").replace("](/steps/", "](../docs/steps/")
body = lambda s: s.split("\n", 1)[1]
(out / "_S2.1_inbreve.md").write_text(fix(sec("\n# In breve\n")), encoding="utf-8")
for name, head in (("rinviati", "\n# Dubbi rinviati\n"), ("verifiche", "\n# Verifiche di Luca\n"), ("posthoc", "\n# Analisi post hoc"),
                   ("deviazioni", "\n# Deviazioni dichiarate\n"), ("review", "\n# Revisione avversariale\n"), ("perimetro", "\n# Perimetro\n"),
                   ("esiti", "\n# Esiti\n"), ("descr", "\n# Riferimenti descrittivi: risultati")):
    (out / f"_S2.1_{name}.md").write_text(fix(body(sec(head))).replace("\n## ", "\n### "), encoding="utf-8")
pre = sec("\n# Attese pre-registrate\n"); (out / "_S2.1_prereg.md").write_text(fix(body(pre)).replace("\n## ", "\n### "), encoding="utf-8")
inb = sec("\n# In breve\n")
cp = pd.read_csv(R / "S2.1_CP1_roi.csv"); gr = pd.read_csv(R / "S2.1_CP1_grid.csv"); ph1 = pd.read_csv(R / "S2.1_PH1_G2c_cells_only.csv")
ph2 = pd.read_csv(R / "S2.1_PH2_density_area.csv"); ds = pd.read_csv(R / "S2.1_desc_tissue.csv"); tt = pd.read_csv(R / "S2.1_D4_time.csv")
g2 = pd.read_csv(R / "S2.1_G2.csv"); rep = (M / "repro.log").read_text()
ok = cp[cp.grid == "none"]; wr = cp[cp.grid != "none"]
a2 = ph2.groupby("archetype").mean(numeric_only=True)
d2 = lambda a, s: ds[(ds.quantity == "d2c") & (ds.archetype == a) & (ds.source == s)]["mean"].iloc[0]
glass = tt.n_bins_glass_margin / tt.n_bins
rows = [
 ("25 aree di cinque tessuti", "S2.1_CP1_roi.csv, griglia giusta: aree", len(ok), 25, len(ok) == 25),
 ("esattamente i quadratini che Space Ranger gli assegna (accordo)", "S2.1_CP1_roi.csv, jaccard minimo, griglia giusta", ok.jaccard.min(), 1.0, ok.jaccard.min() == 1.0),
 ("nessun quadratino cade a meta'", "S2.1_CP1_roi.csv, frac_mid massimo, griglia giusta", ok.frac_mid.max(), 0.0, ok.frac_mid.max() == 0),
 ("griglie sbagliate al piu' al 70 %", "S2.1_CP1_roi.csv, jaccard massimo, griglie sbagliate", round(wr.jaccard.max(), 4), 0.70, round(wr.jaccard.max(), 2) == 0.70),
 ("passo 2 um entro un millesimo", "S2.1_CP1_grid.csv, pitch_um", f"{gr.pitch_um.min():.5f}-{gr.pitch_um.max():.5f}", "|p/2 - 1| < 0.001", bool((abs(gr.pitch_um / 2 - 1) < 0.001).all())),
 ("blocchi di 8 um come Space Ranger", "S2.1_CP1_grid.csv, n_bm_block_mismatch", int(gr.n_bm_block_mismatch.fillna(0).sum()), 0, gr.n_bm_block_mismatch.fillna(0).sum() == 0),
 ("6 aree su 30 con regioni senza cellule", "S2.1_PH1_G2c_cells_only.csv, n_regions_no_cells > 0", int((ph1.n_regions_no_cells > 0).sum()), 6, (ph1.n_regions_no_cells > 0).sum() == 6),
 ("macchia di 560 um2 nel cervello", "S2.1_PH1_G2c_cells_only.csv, A5_r1 area_no_cells_um2", round(ph1.loc[ph1.tag == "A5_r1", "area_no_cells_um2"].iloc[0], 1), 560, round(ph1.loc[ph1.tag == "A5_r1", "area_no_cells_um2"].iloc[0]) == 560),
 ("632 um2 su 26 mm2", "S2.1_PH1_G2c_cells_only.csv, somme", f"{ph1.area_no_cells_um2.sum():.1f} / {ph1.area_tissue_um2.sum() / 1e6:.2f}", "632 / 26", round(ph1.area_no_cells_um2.sum()) == 632 and round(ph1.area_tissue_um2.sum() / 1e6) == 26),
 ("rifatto sulle regioni con cellule torna in tutte e 30", "S2.1_PH1_G2c_cells_only.csv, pass", int(ph1["pass"].sum()), 30, ph1["pass"].sum() == 30),
 ("controllo come scritto: 24 aree su 30", "S2.1_G2.csv, pass", int(g2["pass"].sum()), 24, g2["pass"].sum() == 24),
 ("file identici in un'altra cartella", "repro.log", "PASS G-3" in rep, True, "PASS G-3" in rep),
 ("stroma tre quarti, cervello 1.2 volte (quadratini con centro nel nucleo)", "S2.1_PH2_density_area.csv, ratio_d1 medio A3 e A5", f"{a2.loc['A3', 'ratio_d1']:.3f} / {a2.loc['A5', 'ratio_d1']:.3f}", "0.75 / 1.2", round(a2.loc['A3', 'ratio_d1'], 2) in (0.75, 0.76) and round(a2.loc['A5', 'ratio_d1'], 1) == 1.2),
 ("taglia dei nuclei entro il 4 %", "S2.1_PH2_density_area.csv, ratio_area medio per tessuto", f"{a2.ratio_area.min():.3f}-{a2.ratio_area.max():.3f}", "0.96-1.04", a2.ratio_area.min() > 0.96 and a2.ratio_area.max() < 1.04),
 ("nuclei per blocco 1.3 nel tumore simulato come nel vero", "S2.1_desc_tissue.csv, d2c A2 simulato / Space Ranger", f"{d2('A2', 'simulated'):.3f} / {d2('A2', 'Space Ranger 4'):.3f}", "1.3 / 1.3", round(d2('A2', 'simulated'), 1) == 1.3 and round(d2('A2', 'Space Ranger 4'), 1) == 1.3),
 ("3.0 contro 3.3 nel linfonodo", "S2.1_desc_tissue.csv, d2c A4 simulato / Space Ranger", f"{d2('A4', 'simulated'):.3f} / {d2('A4', 'Space Ranger 4'):.3f}", "3.0 / 3.3", round(d2('A4', 'simulated'), 1) == 3.0 and round(d2('A4', 'Space Ranger 4'), 1) == 3.3),
 ("circa un minuto, al piu' 1.5 GB", "S2.1_D4_time.csv, elapsed_s mediana e vm_hwm_mb massimo", f"{tt.elapsed_s.median():.1f} s / {tt.vm_hwm_mb.max():.0f} MB", "~60 s / <= 1.5 GB", 40 < tt.elapsed_s.median() < 80 and tt.vm_hwm_mb.max() <= 1536),
 ("fra il 59 e il 73 % di vetro", "S2.1_D4_time.csv, n_bins_glass_margin / n_bins", f"{glass.min():.3f}-{glass.max():.3f}", "0.59-0.73", round(glass.min(), 2) == 0.59 and round(glass.max(), 2) == 0.73),
 ("46-81 secondi, mediana 58", "S2.1_D4_time.csv, elapsed_s", f"{tt.elapsed_s.min():.1f}-{tt.elapsed_s.max():.1f}, {tt.elapsed_s.median():.1f}", "46-81, 58", round(tt.elapsed_s.min()) == 46 and round(tt.elapsed_s.max()) == 81 and round(tt.elapsed_s.median()) == 58),
 ("carico medio fra 14 e 17 processi", "S2.1_D4_time.csv, load1", f"{tt.load1.min():.1f}-{tt.load1.max():.1f}", "14-17", round(tt.load1.min()) == 14 and round(tt.load1.max()) == 17),
]
V = pd.DataFrame(rows, columns=["frase", "fonte", "valore letto", "valore nel testo", "coincide"])
V.to_csv(R / "S2.1_verification_log.csv", index=False)
md = "| frase del riassunto | fonte (tabella o file, riga) | valore letto | valore nel testo | coincide |\n|---|---|---|---|---|\n" + "".join(
    f"| {r.frase} | {r.fonte} | {r._3} | {r._4} | {'sì' if r.coincide else 'NO'} |\n" for r in V.itertuples())
(out / "_S2.1_fonti.md").write_text(md, encoding="utf-8")
assert V.coincide.all(), V[~V.coincide]
print(f"{len(V)} numeri del riassunto verificati, tutti coincidono")
