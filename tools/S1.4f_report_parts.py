#!/usr/bin/env python3
"""tools/S1.4f_report_parts.py - estrae dal file dello step (docs/steps/S1.4f.md) le sezioni incluse nel report e costruisce la tabella
«Da dove vengono i numeri del riassunto»: ogni numero di «In breve» e' riletto dai file dei risultati e confrontato con il testo (si ferma
se uno non coincide). Uso (root del repo): .venv/bin/python -B tools/S1.4f_report_parts.py
Uscite: reports/_S1.4f_*.md, results/S1.4f/S1.4f_verification_log.csv"""
import pathlib, re, numpy as np, pandas as pd
R = pathlib.Path("results/S1.4f"); out = pathlib.Path("reports"); LOGS = pathlib.Path("/mnt/micron/geo_spatialtrans/S1.4f/logs")
full = pathlib.Path("docs/steps/S1.4f.md").read_text(encoding="utf-8"); t = "\n" + full.split("\n---\n", 1)[1]
def sec(start, stop_re):
    i = t.index(start); m = re.search(stop_re, t[i + len(start):]); j = i + len(start) + (m.start() if m else len(t)); return t[i:j].strip() + "\n"
fix = lambda s: s.replace("](/decisions/", "](../docs/decisions/").replace("](/steps/", "](../docs/steps/")
(out / "_S1.4f_inbreve.md").write_text(fix(sec("\n# In breve\n", r"\n# ")), encoding="utf-8")
(out / "_S1.4f_rinviati.md").write_text(fix(sec("\n# Dubbi rinviati\n", r"\n# ").replace("# Dubbi rinviati\n", "", 1)), encoding="utf-8")
(out / "_S1.4f_verifiche.md").write_text(sec("\n# Verifiche di Luca\n", r"\n# ").replace("# Verifiche di Luca\n", "", 1), encoding="utf-8")
(out / "_S1.4f_dopo.md").write_text(fix(sec("\n# Dopo il risultato\n", r"\n# In breve\n").replace("# Dopo il risultato\n", "", 1)).replace("\n## ", "\n### "), encoding="utf-8")
(out / "_S1.4f_prereg.md").write_text(fix(sec("\n# Attese pre-registrate\n", r"\n# Dopo il risultato\n").replace("# Attese pre-registrate\n", "", 1)).replace("\n## ", "\n### "), encoding="utf-8")
(out / "_S1.4f_explore.md").write_text(fix(sec("\n# Esplorazione dichiarata (prima della pre-registrazione)\n", r"\n# ").split("\n", 2)[2]), encoding="utf-8")
(out / "_S1.4f_received.md").write_text(fix(sec("\n# Dubbi ricevuti\n", r"\n# ").split("\n", 2)[2]), encoding="utf-8")

V = pd.read_csv(R / "S1.4f_verdicts.csv"); RN = pd.read_csv(R / "S1.4f_ruleN2.csv").set_index("arm"); CP = pd.read_csv(R / "S1.4f_cp.csv"); D = pd.read_csv(R / "S1.4f_desc.csv")
RF = pd.read_csv(R / "S1.4f_refs.csv"); PH = pd.read_csv(R / "S1.4f_PH1_territory_nn.csv"); S1 = pd.read_csv(R / "S1.4f_S1.csv")
PHR = pd.read_csv("results/S1.4e/S1.4e_PH1_raters.csv")
v = lambda i, a, arm, cls=None: V[(V.id == i) & (V.archetype == a) & (V.arm == arm) & ((V.cls == cls) if cls else True)].iloc[0]
f2 = lambda x: f"{x:.2f}"; dense = ["A1", "A2", "A3", "A4"]
quota = D[D.desc == "P-1 quota dell'effetto di ET conservata da EP"].set_index("archetype").est
cp1 = CP[(CP.id == "CP-1")]; ep = cp1[cp1.contrast == "EP-K2"].set_index("archetype").est; et = cp1[cp1.contrast == "ET-K2"].set_index("archetype").est
terr = PH[(PH.what == "territory") & (PH.cls == "nn6-8") & PH.archetype.isin(dense)].groupby(["source", "archetype", "roi_id"]).frac.mean().groupby(["source", "archetype"]).mean()
ph = lambda r, a: PHR[(PHR.what == "contatto sulle coppie comuni") & (PHR.rater == r) & (PHR.archetype == a) & (PHR.cls == "lt6")].value.iloc[0]
RUNLOG = {"calB": R / "S1.4f_chain_run1.log", "genB": R / "S1.4f_chain_run1.log", "genA": R / "S1.4f_chain_run2.log"}   # la ripresa salta le unita' fatte (0.0 min)
glog = lambda s: re.search(r"\[" + s + r"\] \d+ unita', 0 errori, ([\d.]+) min", RUNLOG[s].read_text()).group(1)
rng = lambda xs, k=100: f"{round(k * min(xs)):.0f}-{round(k * max(xs)):.0f}"
rows = [
 ("intestino 0.57 → 0.51", "S1.4f_verdicts.csv, P-1 intestino, K2 ed EP", f"{f2(v('P-1','A1','K2').est)} {f2(v('P-1','A1','EP').est)}", "0.57 0.51"),
 ("tumore del colon 0.69 → 0.64 (vero al più 0.57)", "S1.4f_verdicts.csv, P-1 tumore del colon", f"{f2(v('P-1','A2','K2').est)} {f2(v('P-1','A2','EP').est)} {f2(v('P-1','A2','EP').hi)}", "0.69 0.64 0.57"),
 ("stroma 0.36 → 0.34 (al più 0.28)", "S1.4f_verdicts.csv, P-1 stroma", f"{f2(v('P-1','A3','K2').est)} {f2(v('P-1','A3','EP').est)} {f2(v('P-1','A3','EP').hi)}", "0.36 0.34 0.28"),
 ("linfonodo 0.83 → 0.78 (al più 0.75)", "S1.4f_verdicts.csv, P-1 linfonodo", f"{f2(v('P-1','A4','K2').est)} {f2(v('P-1','A4','EP').est)} {f2(v('P-1','A4','EP').hi)}", "0.83 0.78 0.75"),
 ("circa metà dei contatti tolti dall'orientamento esatto; nello stroma un quarto", "S1.4f_desc.csv, quota dell'effetto di ET conservata da EP (intestino, tumore, linfonodo; stroma)",
  f"{quota[['A1','A2','A4']].min():.2f}-{quota[['A1','A2','A4']].max():.2f} {quota['A3']:.2f}", "0.52-0.56 0.24"),
 ("tumore del colon e linfonodo fuori dal vero, stroma al limite", "S1.4f_verdicts.csv, P-1 EP", " ".join(v("P-1", a, "EP").verdict for a in ["A2", "A4", "A3"]), "smentita smentita indeterminata"),
 ("3 misure di contatto su 10, come il disco; ne servivano almeno 2 in più", "S1.4f_ruleN2.csv", f"{RN.loc['EP'].score_P1_P2} {RN.loc['K2'].score_P1_P2} {RN.loc['K2'].winner}", "3 3 K2"),
 ("allineamento dei nuclei vicini giusto in 5 tessuti su 5", "S1.4f_verdicts.csv, P-3 EP", str(int((V[(V.id == 'P-3') & (V.arm == 'EP')].verdict == 'confermata').sum())), "5"),
 ("con tutte e tre le regole troppo allineati a 5-8 µm", "S1.4f_verdicts.csv, P-5 caselle smentite per V0, V1, V2", " ".join(str(int((V[(V.id == 'P-5') & (V.arm == a)].verdict == 'smentita').sum())) for a in ["V0", "V1", "V2"]), "6 4 4"),
 ("intestino a 6-8 µm: veri metà di traverso, come a caso; simulati 62 su 100", "S1.4f_refs.csv R-4 intestino nn6-8 (stima); S1.4f_verdicts.csv P-5 EP intestino nn6-8",
  f"{RF[(RF.archetype == 'A1') & (RF.cls == 'nn6-8')].est.iloc[0]:.2f} {round(100 * v('P-5', 'A1', 'EP', 'nn6-8').est)}", "0.49 62"),
 ("asse del territorio di traverso circa 8 volte su 10 nei veri e nei simulati (6-8 µm)", "S1.4f_PH1_territory_nn.csv, territori, tessuti densi, nn6-8",
  f"{rng(terr['real'])} {rng(terr['sim'])}", "81-82 80-84"),
 ("sui tuoi punti toglie da 3 a 7 coppie su 100 sotto i 6 µm", "S1.4f_cp.csv, CP-1 EP-K2", rng(-ep), "3-7"),
 ("circa metà dell'orientamento esatto", "S1.4f_cp.csv, CP-1 (EP-K2)/(ET-K2)", f"{(ep / et).min():.2f}-{(ep / et).max():.2f}", "0.44-0.48"),
 ("togliere i nuclei contati due volte non cambia il risultato", "S1.4f_S1.csv, contrasti entro 2 errori standard", f"{int(S1.within_2se.sum())}/{len(S1)}", "8/8"),
 ("tu 17 su 100, Federica 44 (stesse coppie vicine del tumore)", "results/S1.4e/S1.4e_PH1_raters.csv, tumore < 6 µm", f"{round(100 * ph('Luca', 'A2'))} {round(100 * ph('Federica', 'A2'))}", "17 44"),
 ("taratura 24,5 minuti, 300 simulazioni 12,2, 24 finestre 8,3", "S1.4f_chain_run1.log (calB, genB), S1.4f_chain_run2.log (genA)", f"{glog('calB')} {glog('genB')} {glog('genA')}", "24.5 12.2 8.3"),
]
bad = [r for r in rows if r[2] != r[3]]
extra = [("catena circa un'ora con 20 processori", "registro dei job del server: 2703 s (fino all'arresto dello stadio delle coppie) + 794 s (ripresa)", "58 min", "58 min"),
         ("test, mutanti e regressioni circa mezz'ora per passata; tre passate", "registro dei job del server: 2278 s, 2170 s, 2046 s", "34-38 min, 3", "34-38 min, 3")]
pd.DataFrame(rows + extra, columns=["frase", "fonte", "valore letto", "testo"]).to_csv(R / "S1.4f_verification_log.csv", index=False)
tab = "| frase del riassunto | fonte (tabella o file, riga) | valore letto |\n|---|---|---|\n" + "".join(f"| {a} | {b} | {c} |\n" for a, b, c, _ in rows + extra)
(out / "_S1.4f_fonti.md").write_text(tab, encoding="utf-8")
if bad:
    for r in bad: print("DIVERSO:", r)
    raise SystemExit(1)
print(f"{len(rows)} frasi verificate dai file, tutte coincidono (+ {len(extra)} dal registro dei job)")
