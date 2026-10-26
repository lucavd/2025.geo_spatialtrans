#!/usr/bin/env python3
"""tools/S1.5e_report_parts.py - estrae dal file dello step (docs/steps/S1.5e.md) le sezioni incluse nel report e costruisce la tabella
«Da dove vengono i numeri del riassunto» e il registro delle verifiche: ogni numero scritto in «In breve» viene riletto dai file dei
risultati e confrontato con il testo (si ferma se uno non coincide).
Uso (root del repo): .venv/bin/python -B tools/S1.5e_report_parts.py
Uscite: reports/_S1.5e_*.md, results/S1.5e/S1.5e_verification_log.csv"""
import pathlib, re, numpy as np, pandas as pd
R = pathlib.Path("results/S1.5e"); out = pathlib.Path("reports"); LOGS = pathlib.Path("/mnt/micron/geo_spatialtrans/S1.5e")
full = pathlib.Path("docs/steps/S1.5e.md").read_text(encoding="utf-8"); t = "\n" + full.split("\n---\n", 1)[1]
def sec(start, stop_re):
    i = t.index(start); m = re.search(stop_re, t[i + len(start):]); j = i + len(start) + (m.start() if m else len(t))
    return t[i:j].strip() + "\n"
fix = lambda s: s.replace("](/decisions/", "](../docs/decisions/").replace("](../decisions/", "](../docs/decisions/").replace("](/steps/", "](../docs/steps/")
(out / "_S1.5e_inbreve.md").write_text(fix(sec("\n# In breve\n", r"\n# ")), encoding="utf-8")
(out / "_S1.5e_rinviati.md").write_text(sec("\n# Dubbi rinviati\n", r"\n# ").replace("# Dubbi rinviati\n", "", 1), encoding="utf-8")
(out / "_S1.5e_verifiche.md").write_text(sec("\n# Verifiche di Luca\n", r"\n# ").replace("# Verifiche di Luca\n", "", 1), encoding="utf-8")
(out / "_S1.5e_dopo.md").write_text(fix(sec("\n# Dopo il risultato\n", r"\n# In breve\n").replace("# Dopo il risultato\n", "", 1)).replace("\n## ", "\n### "), encoding="utf-8")
(out / "_S1.5e_prereg.md").write_text(fix(sec("\n# Attese pre-registrate\n", r"\n# Dopo il risultato\n").replace("# Attese pre-registrate\n", "", 1)).replace("\n## ", "\n### "), encoding="utf-8")
(out / "_S1.5e_explore.md").write_text(fix(sec("\n# Esplorazione fatta prima della pre-registrazione (dichiarata)\n", r"\n# ").split("\n", 2)[2]), encoding="utf-8")
(out / "_S1.5e_received.md").write_text(fix(sec("\n# Dubbi ricevuti\n", r"\n# ").split("\n", 2)[2]), encoding="utf-8")

V = pd.read_csv(R / "S1.5e_verdicts.csv"); RA = pd.read_csv(R / "S1.5e_ruleA.csv").iloc[0]; RQ = pd.read_csv(R / "S1.5e_ruleQ.csv").iloc[0]
U = pd.read_csv(R / "S1.5e_panel_units.csv"); U0 = pd.read_csv("results/S1.5d/S1.5d_panel_units.csv"); CP = pd.read_csv(R / "S1.5e_cp.csv")
E1 = pd.read_csv(R / "S1.5e_explore_density_ref.csv"); E3 = pd.read_csv(R / "S1.5e_explore_contacts_summary.csv"); D = pd.read_csv(R / "S1.5e_desc.csv")
PH = pd.read_csv(R / "S1.5e_posthoc_loss.csv"); HX = pd.read_csv(R / "S1.5e_het_XA.csv"); T = (R / "S1.5e_test_S1.5e.log").read_text()
cal = (LOGS / "calpos_gen.log").read_text()
def v(i, a, arm="RV"): return V[(V.id == i) & (V.archetype == a) & (V.unit == "tessuto") & (V.arm == arm)].iloc[0]
NEW = ["X", "XA", "XAD", "XADs"]; pool = lambda a: HX[(HX.archetype == a) & (HX.fold == "pooled")].iloc[0]
nf = U[U.n_failed > 0].iloc[0]; key = lambda d: d[(d.arm == "Q") & (d.archetype == nf.archetype) & (d.roi_id == nf.roi_id) & (d.seed == nf.seed)].md5.iloc[0]
cp3 = CP[CP.cp.str.startswith("CP-3") & CP.item.str.split().str[0].isin(NEW)]; cp3a = cp3[cp3.item.str.endswith("A4")]; cp3o = cp3[~cp3.item.str.endswith("A4")]
cpv = lambda s: float(s.split()[0]); e3 = lambda src: E3[(E3.archetype == "A2") & (E3.source == src)].s1_b0.iloc[0]
rm = lambda a, arms: D[D.desc.str.startswith("D-desc-1b") & (D.archetype == a) & D.arm.isin(arms)].est
ph6 = PH[PH.archetype == "A6"].set_index("which"); r1 = lambda w: ph6.loc[w, "loss_rep1"]
mx = re.search(r"T8b.*senza (\d+), con (\d+)", T).groups(); cm = re.search(r"\[calpos\] (\d+) unita', (\d+) errori, ([\d.]+) min", cal).groups()
rng = lambda xs, f="{:.0f}": f"{f.format(min(xs))}-{f.format(max(xs))}"
rows = [
 ("6 cellule su 8 980 in un solo tessuto simulato dell'intestino", "S1.5e_panel_units.csv, n_failed > 0 (n + n_failed)",
  f"{int(U.n_failed.sum())} {int(nf.n + nf.n_failed)} {int((U.n_failed > 0).sum())} {nf.archetype}", "6 8980 1 A1"),
 ("lo stesso tessuto simulato dello step precedente, identico bit per bit", "S1.5e_panel_units.csv contro results/S1.5d/S1.5d_panel_units.csv, md5 di Q intestino r5 seme 5",
  str(key(U) == key(U0)), "True"),
 ("le varianti nuove non perdono nessuna cellula", "S1.5e_panel_units.csv, n_failed in X, XA, XAD, XADs", f"{int(U[U.arm.isin(NEW)].n_failed.sum())}", "0"),
 ("le varianti nuove 15 o 16 misure su 40, quella di oggi 19", "S1.5e_ruleA.csv", f"{RA.winner} {RA.RV} {RA.X} {RA.XA} {RA.XADs} {RA.sens_XAD} {RA.n_verdicts}", "RV 19 15 15 16 16 40"),
 ("nel caso che perdeva 533 cellule, con il tetto nessuna", "S1.5e_test_S1.5e.log, T8b", f"{mx[0]} {mx[1]}", "533 0"),
 ("vicino a meno di 5 µm: si toccano 96 nuclei simulati su 100 contro 47 veri", "S1.5e_explore_contacts_summary.csv, tumore del colon, classe 0-5 µm, RV e manuale (17 focali)",
  f"{100 * e3('RV'):.0f} {100 * e3('manual'):.0f} {E3[(E3.archetype == 'A2') & (E3.source == 'manual')].n_b0.iloc[0]}", "96 47 17"),
 ("tumore del colon: circa due su tre in ogni variante, contro al più 57 su 100", "S1.5e_verdicts.csv, P-1 tumore del colon, RV e varianti nuove",
  f"{rng([v('P-1', 'A2', a).est for a in ['RV'] + NEW], '{:.2f}')} {100 * v('P-1', 'A2').hi:.0f}", "0.62-0.69 57"),
 ("intestino: spazio vuoto 2 % delle finestre e 8 % del tessuto", "S1.5e_explore_density_ref.csv, intestino, media delle finestre",
  f"{100 * E1[E1.archetype == 'A1'].void_frac_window.mean():.0f} {100 * E1[E1.archetype == 'A1'].void_frac_roi.mean():.0f}", "2 8"),
 ("densità del conteggio manuale circa 12 400, più alta del tessuto intero", "S1.5e_verdicts.csv, P-15 intestino RV (densità del preset) e intervallo vero",
  f"{v('P-15', 'A1').est:.0f} > {v('P-15', 'A1').hi:.0f}", "12419 > 11522"),
 ("vuoti allungati: circa 9 970, appena sotto il minimo vero 9 975", "S1.5e_verdicts.csv, P-15 intestino XA", f"{v('P-15', 'A1', 'XA').est:.0f} {v('P-15', 'A1', 'XA').lo:.0f}", "9974 9975"),
 ("linfonodo: vuoti veri 0.2 % del tessuto", "S1.5e_explore_density_ref.csv, linfonodo", f"{100 * E1[E1.archetype == 'A4'].void_frac_roi.mean():.1f}", "0.2"),
 ("linfonodo: la calibrazione sceglie il 30 % di vuoti, più stretti della distanza fra due cellule", "S1.5e_het_XA.csv, linfonodo pooled (asse corto = ℓ / √κ)",
  f"{pool('A4').void_frac:.2f} {pool('A4').void_ell_um / np.sqrt(pool('A4').void_aspect):.1f}", "0.30 2.0"),
 ("linfonodo: toglie un terzo delle cellule", "S1.5e_desc.csv, D-desc-1b linfonodo, X e XA", rng(rm("A4", ["X", "XA"]), "{:.2f}"), "0.30-0.32"),
 ("linfonodo: circa 19 000 contro almeno 21 500", "S1.5e_verdicts.csv, P-15 linfonodo X, XA, XADs", f"{rng([v('P-15', 'A4', a).est for a in ['X', 'XA', 'XADs']])} {v('P-15', 'A4').lo:.0f}", "18933-19642 21533"),
 ("cervello: circa 800 contro almeno 900", "S1.5e_verdicts.csv, P-15 cervello XA, XADs", f"{rng([v('P-15', 'A5', a).est for a in ['XA', 'XADs']])} {v('P-15', 'A5').lo:.0f}", "795-799 917"),
 ("cuore: 15 % di vuoti dove i veri sono meno dell'1 %", "S1.5e_het_XA.csv, cuore pooled; S1.5e_explore_density_ref.csv, cuore",
  f"{pool('A6').void_frac:.2f} {100 * E1[E1.archetype == 'A6'].void_frac_roi.mean():.1f}", "0.15 0.7"),
 ("cuore: con due ripetizioni le due scelte erano quasi pari (prima ripetizione)", "S1.5e_posthoc_loss.csv, cuore, loss_rep1",
  " ".join(f"{r1(w):.1f}" for w in ph6.index), "52.9 55.0"),
 ("cervello: area in vuoti molto allungati da un sesto a circa un quarto, contro almeno tre decimi", "S1.5e_verdicts.csv, P-16 cervello RV, XA, XAD, XADs",
  f"{v('P-16', 'A5').est:.3f} {rng([v('P-16', 'A5', a).est for a in ['XA', 'XAD', 'XADs']], '{:.3f}')} {v('P-16', 'A5').lo:.3f}", "0.166 0.245-0.286 0.305"),
 ("intestino: 43-48 vuoti contro 96-117", "S1.5e_verdicts.csv, P-17 intestino XA, XADs", f"{rng([v('P-17', 'A1', a).est for a in ['XA', 'XADs']])} {v('P-17', 'A1').lo:.0f}-{v('P-17', 'A1').hi:.0f}", "43-48 96-117"),
 ("intestino: più di metà dell'area vuota in vuoti molto allungati contro al più il 45 %", "S1.5e_verdicts.csv, P-16 intestino XA, XADs",
  f"{rng([v('P-16', 'A1', a).est for a in ['XA', 'XADs']], '{:.3f}')} {v('P-16', 'A1').hi:.3f}", "0.524-0.548 0.450"),
 ("contatti e distanze: tiene in 17 casi su 20 fuori dal linfonodo", "S1.5e_cp.csv, CP-3 sulle varianti nuove, senza il linfonodo", f"{int(cp3o['pass'].sum())}/{len(cp3o)}", "17/20"),
 ("linfonodo: si toccano 6-8 su 100 più di quanto dicono le distanze", "S1.5e_cp.csv, CP-3 linfonodo, varianti nuove", rng([100 * cpv(x) for x in cp3a.value]) + f" {int((~cp3a['pass']).sum())}/{len(cp3a)}", "6-8 4/4"),
 ("senza numero imposto 60 contro 58 e 48 contro 46", "S1.5e_ruleQ.csv", f"{RQ.spatial_noquotas} {RQ.spatial_quotas} {RQ.types_noquotas} {RQ.types_quotas} {RQ.result}", "60 58 48 46 senza quote"),
 ("6 720 semine in 59 minuti, 0 errori", "/mnt/micron/geo_spatialtrans/S1.5e/calpos_gen.log", f"{cm[0]} {cm[1]} {float(cm[2]):.0f}", "6720 0 59"),
]
L = pd.DataFrame(rows, columns=["frase", "fonte", "valore_letto", "valore_scritto"]); L["coincide"] = L.valore_letto == L.valore_scritto
L.to_csv(R / "S1.5e_verification_log.csv", index=False)
bad = L[~L.coincide]; assert bad.empty, bad.to_string()
tab = "| frase del riassunto | fonte | valore letto |\n|---|---|---|\n" + "\n".join(f"| {a} | {b} | {c} |" for a, b, c, _ in rows) + "\n"
(out / "_S1.5e_fonti.md").write_text(tab, encoding="utf-8"); print(f"{len(L)} righe, {int(L.coincide.sum())} coincidono")
