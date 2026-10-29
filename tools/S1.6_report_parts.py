#!/usr/bin/env python3
"""tools/S1.6_report_parts.py - estrae dal file dello step (docs/steps/S1.6.md) le sezioni incluse nel report e costruisce la tabella
«Da dove vengono i numeri del riassunto»: ogni numero di «In breve» e' riletto dai file dei risultati e dai log e confrontato con il testo
(si ferma se uno non coincide). Uso (root del repo): .venv/bin/python -B tools/S1.6_report_parts.py
Uscite: reports/_S1.6_*.md, results/S1.6/S1.6_verification_log.csv"""
import pathlib, re, os, datetime, pandas as pd
R = pathlib.Path("results/S1.6"); out = pathlib.Path("reports"); M = pathlib.Path("/mnt/micron/geo_spatialtrans/S1.6")
full = pathlib.Path("docs/steps/S1.6.md").read_text(encoding="utf-8"); t = "\n" + full.split("\n---\n", 1)[1]
def sec(start, stop_re):
    i = t.index(start); m = re.search(stop_re, t[i + len(start):]); j = i + len(start) + (m.start() if m else len(t)); return t[i:j].strip() + "\n"
fix = lambda s: s.replace("](/decisions/", "](../docs/decisions/").replace("](/steps/", "](../docs/steps/")
(out / "_S1.6_inbreve.md").write_text(fix(sec("\n# In breve\n", r"\n# ")), encoding="utf-8")
(out / "_S1.6_rinviati.md").write_text(fix(sec("\n# Dubbi rinviati\n", r"\n# ").replace("# Dubbi rinviati\n", "", 1)), encoding="utf-8")
(out / "_S1.6_verifiche.md").write_text(sec("\n# Verifiche di Luca\n", r"\n# ").replace("# Verifiche di Luca\n", "", 1), encoding="utf-8")
(out / "_S1.6_dopo.md").write_text(fix(sec("\n# Dopo il risultato\n", r"\n# In breve\n").replace("# Dopo il risultato\n", "", 1)).replace("\n## ", "\n### "), encoding="utf-8")
(out / "_S1.6_prereg.md").write_text(fix(sec("\n# Attese pre-registrate\n", r"\n# Dopo il risultato\n").replace("# Attese pre-registrate\n", "", 1)).replace("\n## ", "\n### "), encoding="utf-8")
(out / "_S1.6_explore.md").write_text(fix(sec("\n# Esplorazione dichiarata (prima della pre-registrazione)\n", r"\n# ").split("\n", 2)[2]), encoding="utf-8")
(out / "_S1.6_received.md").write_text(fix(sec("\n# Dubbi ricevuti\n", r"\n# ").split("\n", 2)[2]), encoding="utf-8")
(out / "_S1.6_perimetro.md").write_text(fix(sec("\n# Perimetro\n", r"\n# ").split("\n", 2)[2]), encoding="utf-8")

G3 = pd.read_csv(R / "S1.6_G3_check.csv").set_index("check").value
P1 = pd.read_csv(R / "S1.6_P1.csv"); CP1 = pd.read_csv(R / "S1.6_CP1.csv"); CP2 = pd.read_csv(R / "S1.6_CP2.csv").set_index("archetype")
V = pd.read_csv(R / "S1.6_verdicts.csv"); E = pd.read_csv(R / "S1.6_explore_refs_verdicts.csv")
K2 = pd.read_csv(R / "S1.6_explore_parts_S1.4e_K2.csv").iloc[0]
g2 = (M / "tests2/G2.log").read_text(); rev = (M / "review_code/rapporto.md").read_text()
gen = re.search(r"\[gen\] (\d+) unita', 0 errori, ([\d.]+) min", (M / "logs/gen.log").read_text())
mt = lambda p: datetime.datetime.fromtimestamp(os.path.getmtime(p))
bench_min = (mt(M / "logs/md5after.log") - mt(M / "logs/md5before.log")).total_seconds() / 60
def span(d):
    f = [M / d / x for x in os.listdir(M / d) if (M / d / x).is_file()]; ts = sorted(mt(x) for x in f); return (ts[-1] - ts[0]).total_seconds() / 60
workers = (M / "WORKERS").read_text().strip()
a5 = CP1[(CP1.id == "P-15") & (CP1.archetype == "A5")].iloc[0]
others = CP2.drop("A5").ratio
u6 = lambda a, u: V[(V.id == "P-14") & (V.archetype == a) & (V.unit == u) & (V.arm == "U6")].verdict.iloc[0]
ev = lambda a, u: E[(E.id == "P-14") & (E.archetype == a) & (E.unit == u) & (E.arm == "Q")].iloc[0]
eB, eT = ev("A4", "linfociti B"), ev("A4", "linfociti T")
rows = [
 ("300 simulazioni su 300 identiche bit per bit", "S1.6_G3_check.csv, U6_equals_S1.5e_U_md5", str(G3["U6_equals_S1.5e_U_md5"]), "300"),
 ("2 avvisi su 300 simulazioni prima, nessuno ora", "S1.6_G3_check.csv, desc_C7_warn_U_S1.5e e desc_C7_warn_U6", f"{G3['desc_C7_warn_U_S1.5e']} {G3['desc_C7_warn_U6']}", "2 0"),
 ("179 misure: 169 confermate, 10 dove la precisione non basta", "S1.6_P1.csv, esiti", f"{len(P1)} {(P1.verdict == 'confermata').sum()} {(P1.verdict == 'indeterminata').sum()} {(P1.verdict == 'smentita').sum()}", "179 169 10 0"),
 ("dall'8 al 10 % di proposte in più in cinque tessuti, 11 % in meno nel cervello", "S1.6_CP2.csv, ratio", f"{round(100 * (others.min() - 1))}-{round(100 * (others.max() - 1))} {round(100 * (1 - CP2.loc['A5'].ratio))}", "8-10 11"),
 ("cuore, predefinita: fibroblasti e immunitari indeterminati", "S1.6_verdicts.csv, P-14 U6 cuore", f"{u6('A6', 'fibroblasti')} {u6('A6', 'immunitari')}", "indeterminata indeterminata"),
 ("zone dei tipi: linfociti B migliorano, T peggiorano", "S1.6_explore_refs_verdicts.csv, P-14 Q linfonodo (vecchio → nuovo)", f"B {eB.verdict_old}→{eB.verdict_new}, T {eT.verdict_old}→{eT.verdict_new}", "B smentita→confermata, T confermata→smentita"),
 ("nuclei davvero in due parti 180 su quasi due milioni e mezzo", "S1.6_explore_parts_S1.4e_K2.csv, n_multipart, n_nuclei", f"{K2.n_multipart} {K2.n_nuclei}", "180 2478884"),
 ("9 alterazioni su 10 del revisore passavano inosservate; ora i test le vedono tutte", "review_code/rapporto.md (punto 4); tests2/G2.log",
  ("1 " if "ne rileva 1" in rev else "? ") + re.search(r"mutanti di S1.6 rilevati (\d+/\d+)", g2).group(1), "1 16/16"),
 ("in tutto in 34 misure su 179", "S1.6_CP1.csv, smentite", f"{(CP1.verdict == 'smentita').sum()} {len(CP1)}", "34 179"),
 ("cervello: differenza circa 390 nuclei per mm², soglia 340, 4 aree confrontate", "S1.6_CP1.csv, P-15 cervello (diff, margin, n_pairs)", f"{-a5['diff']:.0f} {a5.margin:.0f} {a5.n_pairs}", "390 340 4"),
 ("600 simulazioni, con raster e misure, 23 minuti con 20 processori", "logs/md5before.log → md5after.log (ora di modifica); WORKERS; gen.log", f"{bench_min:.0f} {workers} gen {gen.group(2)} min", "23 20 gen 11.9 min"),
 ("test, mutanti e ambiente: due giri di circa un'ora", "tests/ e tests2/ (ora di modifica del primo e dell'ultimo file)", f"{span('tests'):.0f} {span('tests2'):.0f}", "52 55"),
]
bad = [r for r in rows if r[2] != r[3]]
pd.DataFrame(rows, columns=["frase", "fonte", "valore letto", "testo"]).to_csv(R / "S1.6_verification_log.csv", index=False)
tab = "| frase del riassunto | fonte (tabella o file, riga) | valore letto |\n|---|---|---|\n" + "".join(f"| {a} | {b} | {c} |\n" for a, b, c, _ in rows)
(out / "_S1.6_fonti.md").write_text(tab, encoding="utf-8")
if bad:
    for r in bad: print("DIVERSO:", r)
    raise SystemExit(1)
print(f"{len(rows)} frasi verificate dai file, tutte coincidono")
