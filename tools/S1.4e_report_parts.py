#!/usr/bin/env python3
"""tools/S1.4e_report_parts.py - estrae dal file dello step (docs/steps/S1.4e.md) le sezioni incluse nel report e costruisce la tabella
«Da dove vengono i numeri del riassunto»: ogni numero di «In breve» e' riletto dai file dei risultati e confrontato con il testo (si ferma
se uno non coincide). Uso (root del repo): .venv/bin/python -B tools/S1.4e_report_parts.py
Uscite: reports/_S1.4e_*.md, results/S1.4e/S1.4e_verification_log.csv"""
import pathlib, re, numpy as np, pandas as pd
R = pathlib.Path("results/S1.4e"); out = pathlib.Path("reports"); LOGS = pathlib.Path("/mnt/micron/geo_spatialtrans/S1.4e/logs")
full = pathlib.Path("docs/steps/S1.4e.md").read_text(encoding="utf-8"); t = "\n" + full.split("\n---\n", 1)[1]
def sec(start, stop_re):
    i = t.index(start); m = re.search(stop_re, t[i + len(start):]); j = i + len(start) + (m.start() if m else len(t)); return t[i:j].strip() + "\n"
fix = lambda s: s.replace("](/decisions/", "](../docs/decisions/").replace("](/steps/", "](../docs/steps/")
(out / "_S1.4e_inbreve.md").write_text(fix(sec("\n# In breve\n", r"\n# ")), encoding="utf-8")
(out / "_S1.4e_rinviati.md").write_text(sec("\n# Dubbi rinviati\n", r"\n# ").replace("# Dubbi rinviati\n", "", 1), encoding="utf-8")
(out / "_S1.4e_verifiche.md").write_text(sec("\n# Verifiche di Luca\n", r"\n# ").replace("# Verifiche di Luca\n", "", 1), encoding="utf-8")
(out / "_S1.4e_dopo.md").write_text(fix(sec("\n# Dopo il risultato\n", r"\n# In breve\n").replace("# Dopo il risultato\n", "", 1)).replace("\n## ", "\n### "), encoding="utf-8")
(out / "_S1.4e_prereg.md").write_text(fix(sec("\n# Attese pre-registrate\n", r"\n# Dopo il risultato\n").replace("# Attese pre-registrate\n", "", 1)).replace("\n## ", "\n### "), encoding="utf-8")
(out / "_S1.4e_explore.md").write_text(fix(sec("\n# Esplorazione dichiarata (prima della pre-registrazione)\n", r"\n# ").split("\n", 2)[2]), encoding="utf-8")
(out / "_S1.4e_received.md").write_text(fix(sec("\n# Dubbi ricevuti\n", r"\n# ").split("\n", 2)[2]), encoding="utf-8")

V = pd.read_csv(R / "S1.4e_verdicts.csv"); RN = pd.read_csv(R / "S1.4e_ruleN.csv").set_index("arm"); CP = pd.read_csv(R / "S1.4e_cp.csv"); D = pd.read_csv(R / "S1.4e_desc.csv")
REF = pd.read_csv(R / "S1.4e_refs.csv"); PH = pd.read_csv(R / "S1.4e_PH1_raters.csv"); E12 = pd.read_csv(R / "explore/S1.4e_E12_pairs.csv")
v = lambda i, a, arm: V[(V.id == i) & (V.archetype == a) & (V.arm == arm)].iloc[0]
f2 = lambda x: f"{x:.2f}"
dense = ["A1", "A2", "A3", "A4"]; p3t = ["A1", "A2", "A3", "A4", "A6"]
rng = lambda xs, k=100: f"{round(k * min(xs)):.0f}-{round(k * max(xs)):.0f}"
bA = lambda a, arm: D[(D.desc == "banco A < 6 um") & (D.archetype == a) & (D.arm == arm)].est.iloc[0]
ph = lambda r, a: PH[(PH.what == "contatto sulle coppie comuni") & (PH.rater == r) & (PH.archetype == a) & (PH.cls == "lt6")].value.iloc[0]
cp1 = CP[(CP.id == "CP-1") & (CP.contrast == "ET-K2")]; cp2 = CP[(CP.id == "CP-2") & (CP.contrast == "K2")].set_index("archetype")
er_d = D[(D.desc == "P-1 diff vs K2 (per seme)") & (D.arm == "ER") & D.archetype.isin(dense)].est
fr = D[D.desc == "P-3 frazione (vero - ER) / (ET - ER)"].est
g2 = E12[(E12.archetype == "A2") & (E12.d_manual < 5) & (E12.status == "ok")]
rest = PH[PH.what == "coppie di Luca in finestre senza Federica"].set_index("archetype").n
tms = PH[PH.what == "t_ms mediano per item"].set_index("rater").value
glog = lambda s: re.search(r"\[" + s + r"\] \d+ unita', 0 errori, ([\d.]+) min", (LOGS / f"{s}.log").read_text()).group(1)
rows = [
 ("intestino 0.57 → 0.47", "S1.4e_verdicts.csv, P-1 intestino, K2 ed ET", f"{f2(v('P-1','A1','K2').est)} {f2(v('P-1','A1','ET').est)}", "0.57 0.47"),
 ("tumore del colon 0.69 → 0.59 (vero al più 0.57)", "S1.4e_verdicts.csv, P-1 tumore del colon", f"{f2(v('P-1','A2','K2').est)} {f2(v('P-1','A2','ET').est)} {f2(v('P-1','A2','ET').hi)}", "0.69 0.59 0.57"),
 ("stroma 0.36 → 0.28 (al più 0.28)", "S1.4e_verdicts.csv, P-1 stroma", f"{f2(v('P-1','A3','K2').est)} {f2(v('P-1','A3','ET').est)} {f2(v('P-1','A3','ET').hi)}", "0.36 0.28 0.28"),
 ("linfonodo 0.83 → 0.73 (al più 0.75)", "S1.4e_verdicts.csv, P-1 linfonodo", f"{f2(v('P-1','A4','K2').est)} {f2(v('P-1','A4','ET').est)} {f2(v('P-1','A4','ET').hi)}", "0.83 0.73 0.75"),
 ("tre tessuti passano «al limite»", "S1.4e_verdicts.csv, P-1 ET tumore, stroma, linfonodo", " ".join(v("P-1", a, "ET").verdict for a in ["A2", "A3", "A4"]), "indeterminata indeterminata indeterminata"),
 ("4 misure su 10, il disco 3", "S1.4e_ruleN.csv", f"{RN.loc['ET'].score_P1_P2} {RN.loc['K2'].score_P1_P2} {RN.loc['K2'].winner}", "4 3 K2"),
 ("ellisse a caso: al più 0.02 in più", "S1.4e_desc.csv, P-1 diff vs K2 (per seme), ER, tessuti densi", f"{er_d.min():.3f} {er_d.max():.3f}", "0.006 0.020"),
 ("ellisse orientata toglie da 6 a 15 coppie su 100 (sotto 6 µm, stessi centri)", "S1.4e_cp.csv, CP-1 ET-K2", rng(-cp1.est), "6-15"),
 ("di traverso: 77-87 veri", "S1.4e_refs.csv, R-2 (stima)", rng(REF[(REF.ref == "R-2") & REF.archetype.isin(p3t)].est), "77-87"),
 ("86-94 con l'ellisse orientata, 61-75 a caso", "S1.4e_verdicts.csv, P-3 ET ed ER", f"{rng([v('P-3', a, 'ET').est for a in p3t])} {rng([v('P-3', a, 'ER').est for a in p3t])}", "86-94 61-75"),
 ("circa a metà strada", "S1.4e_desc.csv, P-3 frazione (vero - ER) / (ET - ER)", f"{fr.min():.2f}-{fr.max():.2f}", "0.51-0.69"),
 ("intestino 0.33 contro il tuo 0.35", "S1.4e_desc.csv, banco A < 6 um, intestino", f"{f2(bA('A1','ET'))} {f2(bA('A1','manuale'))}", "0.33 0.35"),
 ("tumore del colon 0.72 contro 0.30; stroma 0.74 contro 0.47", "S1.4e_desc.csv, banco A < 6 um", f"{f2(bA('A2','ET'))} {f2(bA('A2','manuale'))} {f2(bA('A3','ET'))} {f2(bA('A3','manuale'))}", "0.72 0.30 0.74 0.47"),
 ("Cellpose 5 volte su 6 contro le tue 2 su 6", "explore/S1.4e_E12_pairs.csv, tumore < 5 µm, coppie agganciate", f"{int(g2.contact_op.sum())}/{len(g2)} {int((g2.answer == 'touch').sum())}/{len(g2)}", "5/6 2/6"),
 ("Federica 44 su 100 contro le tue 17 (una decina di coppie)", "S1.4e_PH1_raters.csv, tumore < 6 µm", f"{round(100 * ph('Federica', 'A2'))} {round(100 * ph('Luca', 'A2'))}", "44 17"),
 ("il disco resta fuori nel tumore del colon e nel linfonodo (CP-2)", "S1.4e_cp.csv, CP-2 K2", f"{cp2.loc['A2'].verdict} {cp2.loc['A4'].verdict}", "smentita smentita"),
 ("circa 320 coppie (233 di tumore, 87 di stroma)", "S1.4e_PH1_raters.csv", f"{rest['A2'] + rest['A3']} {rest['A2']} {rest['A3']}", "320 233 87"),
 ("circa un secondo per coppia", "S1.4e_PH1_raters.csv, t_ms mediano (Federica)", f"{tms['Federica'] / 1000:.1f}", "1.0"),
 ("300 simulazioni 9.6 minuti, 24 finestre 5.8", "logs genB.log, genA.log", f"{glog('genB')} {glog('genA')}", "9.6 5.8"),
]
bad = [r for r in rows if r[2] != r[3]]
pd.DataFrame(rows, columns=["frase", "fonte", "valore letto", "testo"]).to_csv(R / "S1.4e_verification_log.csv", index=False)
tab = "| frase del riassunto | fonte (tabella o file, riga) | valore letto |\n|---|---|---|\n" + "".join(f"| {a} | {b} | {c} |\n" for a, b, c, _ in rows)
(out / "_S1.4e_fonti.md").write_text(tab, encoding="utf-8")
if bad:
    for r in bad: print("DIVERSO:", r)
    raise SystemExit(1)
print(f"{len(rows)} frasi verificate, tutte coincidono")
