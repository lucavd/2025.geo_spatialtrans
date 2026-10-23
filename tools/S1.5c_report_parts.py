#!/usr/bin/env python3
"""tools/S1.5c_report_parts.py - estrae dal file dello step (docs/steps/S1.5c.md) le sezioni incluse nel report e costruisce la tabella
«Da dove vengono i numeri del riassunto» rileggendo ogni valore dai file dei risultati. Uso (root del repo): .venv/bin/python -B tools/S1.5c_report_parts.py"""
import pathlib, re, json, numpy as np, pandas as pd
R = pathlib.Path("results/S1.5c"); out = pathlib.Path("reports")
t = pathlib.Path("docs/steps/S1.5c.md").read_text(encoding="utf-8").split("\n---\n", 1)[1]
def sec(start, stop_re):
    i = t.index(start); m = re.search(stop_re, t[i + len(start):]); j = i + len(start) + (m.start() if m else len(t))
    return t[i:j].strip() + "\n"
fix = lambda s: s.replace("](../decisions/", "](../docs/decisions/")
(out / "_S1.5c_inbreve.md").write_text(fix(sec("\n# In breve\n", r"\n# ")), encoding="utf-8")
(out / "_S1.5c_rinviati.md").write_text(sec("\n# Dubbi rinviati\n", r"\n# ").replace("# Dubbi rinviati\n", "", 1), encoding="utf-8")
(out / "_S1.5c_verifiche.md").write_text(sec("\n# Verifiche di Luca\n", r"\n# ").replace("# Verifiche di Luca\n", "", 1), encoding="utf-8")
(out / "_S1.5c_dopo.md").write_text(fix(sec("\n# Dopo il risultato\n", r"\n# In breve\n").replace("# Dopo il risultato\n", "", 1)), encoding="utf-8")
(out / "_S1.5c_prereg.md").write_text(fix(sec("\n# Attese pre-registrate\n", r"\n# Dopo il risultato\n").replace("# Attese pre-registrate\n", "", 1)), encoding="utf-8")
V = pd.read_csv(R / "S1.5c_verdicts.csv"); CP = pd.read_csv(R / "S1.5c_cp1.csv").set_index("arm"); RU = pd.read_csv(R / "S1.5c_rules.csv")
PH = pd.read_csv(R / "S1.5c_posthoc_terr.csv").groupby(["archetype", "source"]).mean(numeric_only=True); G4 = pd.read_csv(R / "S1.5c_G4_calcheck.csv"); MC = pd.read_csv(R / "S1.5c_G4_mc.csv")
D = pd.read_csv(R / "S1.5c_desc.csv")
def v(i, a, u, arm="N-a"): return V[(V.id == i) & (V.archetype == a) & (V.unit == u) & (V.arm == arm)].iloc[0]
conf = lambda arm, i: int(((V.arm == arm) & (V.id == i) & (V.verdict == "confermata")).sum())
g4a = G4[(G4.part == "a") & G4.checked]; g4b = G4[G4.part == "b"]
tm = D[D.desc.str.startswith("D-desc-4")].pivot(index="archetype", columns="arm", values="est"); ratio = (tm["N-a"] / tm["U"])
rows = [
 ("52 misure su 92, contro 19", "S1.5c_cp1.csv, colonna spatial (N-a, U), n_items", f"{CP.loc['N-a','spatial']} / {CP.loc['N-a','n_items']}; {CP.loc['U','spatial']}"),
 ("raggruppamento 23 su 34 (oggi 4)", "S1.5c_cp1.csv, P-6", f"{CP.loc['N-a','P-6']}; {CP.loc['U','P-6']}; voci {int(((V.id=='P-6')&(V.arm=='N-a')).sum())}"),
 ("distanza dal vicino 10 tipi su 17 (oggi 4)", "S1.5c_cp1.csv, P-5", f"{CP.loc['N-a','P-5']}; {CP.loc['U','P-5']}"),
 ("spazio attorno al nucleo nel tumore del colon e nello stroma", "S1.5c_verdicts.csv, P-3 N-a A2, A3", f"{v('P-3','A2','tessuto').est:.4f} {v('P-3','A2','tessuto').verdict}; {v('P-3','A3','tessuto').est:.4f} {v('P-3','A3','tessuto').verdict}"),
 ("tipi mescolati: le misure tornano a quelle di oggi", "S1.5c_cp1.csv, G contro U: P-4, P-5, P-6", f"G {CP.loc['G','P-4']}, {CP.loc['G','P-5']}, {CP.loc['G','P-6']}; U {CP.loc['U','P-4']}, {CP.loc['U','P-5']}, {CP.loc['U','P-6']}"),
 ("regola: distanza per nucleo; modello del gruppo", "S1.5c_rules.csv, regole 1 e 2", "; ".join(f"{r.a} {r.a_value} contro {r.b} {r.b_value} -> {r.result}" for r in RU.itertuples() if r.rule[0] in "12")),
 ("taglie: 6 tessuti su 6, 14 tipi su 17", "S1.5c_cp1.csv, P-8a, P-8b (N-a)", f"{CP.loc['N-a','P-8a']}; {CP.loc['N-a','P-8b']}"),
 ("territori troppo uniformi in quattro tessuti", "S1.5c_verdicts.csv, P-2 N-a smentite", ", ".join(f"{a} {v('P-2',a,'tessuto').est:.2f}" for a in ["A1","A2","A3","A4","A5","A6"] if v('P-2',a,'tessuto').verdict == "smentita")),
 ("una cellula su venti sopra tre volte il tipico (stroma una su dieci); simulato una su duecento in intestino e tumore, una su venticinque nello stroma", "S1.5c_posthoc_terr.csv, frac_gt3", "; ".join(f"{a} vero {PH.loc[(a,'reale'),'frac_gt3']:.3f}, simulato {PH.loc[(a,'N-a'),'frac_gt3']:.3f}" for a in ["A1","A2","A3"])),
 ("variabilità simulata da circa il 40 % sotto a meno del 20 % sotto", "S1.5c_posthoc_terr.csv, cv e cv_trim", "; ".join(f"{a} {PH.loc[(a,'N-a'),'cv']/PH.loc[(a,'reale'),'cv']-1:+.2f} -> {PH.loc[(a,'N-a'),'cv_trim']/PH.loc[(a,'reale'),'cv_trim']-1:+.2f}" for a in ["A1","A2"])),
 ("linfonodo: 82 su 100 contro al più 75", "S1.5c_verdicts.csv, P-1 N-a A4", f"{v('P-1','A4','tessuto').est:.3f}; {v('P-1','A4','tessuto').hi}"),
 ("cellule quasi senza citoplasma 3,4 % contro al più 2,9 %", "S1.5c_verdicts.csv, P-8c A2", f"N-a {v('P-8c','A2','tessuto').est:.4f}; limite {v('P-8c','A2','tessuto').hi}; U {v('P-8c','A2','tessuto','U').est:.4f}"),
 ("taglie 20 contro 21", "S1.5c_cp1.csv, size", f"{CP.loc['N-a','size']}; {CP.loc['U','size']}"),
 ("fibroblasti dello stroma appena dentro il limite (semina di oggi)", "S1.5c_verdicts.csv, P-8b U A3 fibroblasti", (lambda r: f"{r.est:.3f} + 2 x {r.se:.3f} = {r.est + 2*r.se:.3f} contro {r.hi:.3f}")(v('P-8b','A3','fibroblasti','U'))),
 ("75 confronti su 79 entro la soglia", "S1.5c_G4_calcheck.csv, parte a, bin controllati", f"{int(g4a['pass'].sum())}/{len(g4a)}; parte b {int(g4b['pass'].sum())}/{len(g4b)}"),
 ("con 50 ripetizioni coincidono entro 0,07 nei 79 confronti del controllo", "S1.5c_G4_mc.csv unito a S1.5c_G4_calcheck.csv (bin controllati), |log R_fit_R10 - log py10_mean| massimo", (lambda m: f"{np.nanmax(np.abs(np.log(m.R_fit_R10) - np.log(m.py10_mean))):.3f} su {len(m)} bin; tutti i bin {np.nanmax(np.abs(np.log(MC.R_fit_R10) - np.log(MC.py10_mean))):.3f}")(MC.merge(g4a[["archetype", "cell_type", "bin"]], on=["archetype", "cell_type", "bin"]))),
 ("il revisore si aspetta circa 7 sforamenti per rumore", "S1.5c_review.json, r1, g4_judgement", "7.3 (testo del revisore)"),
 ("tutti i calcoli identici", "S1.5c_review.json, r1, checks", "; ".join(f"{c['n_identical']}/{c['n_compared']}" for c in json.load(open(R / "S1.5c_review.json"))["r1"]["checks"])),
 ("tempo da 1,1 a 4,3 volte; intestino 38 contro 29 s", "S1.5c_desc.csv, D-desc-4", f"{ratio.min():.2f}-{ratio.max():.2f}; A1 {tm.loc['A1','N-a']:.1f} contro {tm.loc['A1','U']:.1f}"),
]
tab = "| frase del riassunto | fonte | valore letto |\n|---|---|---|\n" + "\n".join(f"| {a} | {b} | {c} |" for a, b, c in rows) + "\n"
(out / "_S1.5c_fonti.md").write_text(tab, encoding="utf-8"); print(tab)
