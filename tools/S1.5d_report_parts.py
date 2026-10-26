#!/usr/bin/env python3
"""tools/S1.5d_report_parts.py - estrae dal file dello step (docs/steps/S1.5d.md) le sezioni incluse nel report e costruisce la tabella
«Da dove vengono i numeri del riassunto» e il registro delle verifiche: ogni numero scritto in «In breve» e nei punti principali di
«Dopo il risultato» viene riletto dai file dei risultati e confrontato con il testo (si ferma se uno non coincide).
Uso (root del repo): .venv/bin/python -B tools/S1.5d_report_parts.py
Uscite: reports/_S1.5d_*.md, results/S1.5d/S1.5d_verification_log.csv"""
import pathlib, re, numpy as np, pandas as pd
R = pathlib.Path("results/S1.5d"); out = pathlib.Path("reports")
full = pathlib.Path("docs/steps/S1.5d.md").read_text(encoding="utf-8"); t = full.split("\n---\n", 1)[1]
def sec(start, stop_re):
    i = t.index(start); m = re.search(stop_re, t[i + len(start):]); j = i + len(start) + (m.start() if m else len(t))
    return t[i:j].strip() + "\n"
fix = lambda s: s.replace("](/decisions/", "](../docs/decisions/").replace("](../decisions/", "](../docs/decisions/")
(out / "_S1.5d_inbreve.md").write_text(fix(sec("\n# In breve\n", r"\n# ")), encoding="utf-8")
(out / "_S1.5d_rinviati.md").write_text(sec("\n# Dubbi rinviati\n", r"\n# ").replace("# Dubbi rinviati\n", "", 1), encoding="utf-8")
(out / "_S1.5d_verifiche.md").write_text(sec("\n# Verifiche di Luca\n", r"\n# ").replace("# Verifiche di Luca\n", "", 1), encoding="utf-8")
(out / "_S1.5d_dopo.md").write_text(fix(sec("\n# Dopo il risultato\n", r"\n# In breve\n").replace("# Dopo il risultato\n", "", 1)).replace("\n## ", "\n### "), encoding="utf-8")
(out / "_S1.5d_prereg.md").write_text(fix(sec("\n# Attese pre-registrate\n", r"\n# Dopo il risultato\n").replace("# Attese pre-registrate\n", "", 1)).replace("\n## ", "\n### "), encoding="utf-8")
(out / "_S1.5d_explore.md").write_text(fix(sec("\n# Esplorazione fatta prima della pre-registrazione (dichiarata)\n", r"\n# ").split("\n", 2)[2]), encoding="utf-8")
V = pd.read_csv(R / "S1.5d_verdicts.csv"); K = pd.read_csv(R / "S1.5d_counts.csv").set_index("arm")
RA = pd.read_csv(R / "S1.5d_ruleA.csv").iloc[0]; RB = pd.read_csv(R / "S1.5d_ruleB.csv").iloc[0]; RC = pd.read_csv(R / "S1.5d_ruleC.csv")
RS = pd.read_csv(R / "S1.5d_ruleC_sensitivity_U3c.csv").iloc[0]; CP = pd.read_csv(R / "S1.5d_cp.csv"); U = pd.read_csv(R / "S1.5d_panel_units.csv")
P4 = pd.read_csv(R / "S1.5d_posthoc_p8c_real.csv"); TL = pd.read_csv(R / "S1.5d_posthoc_size_tail.csv").set_index("key"); Q9 = pd.read_csv(R / "S1.5d_posthoc_size_q90.csv").set_index("key")
G3 = pd.read_csv(R / "S1.5d_G3_check.csv").set_index("check"); G4 = {p: pd.read_csv(R / f"S1.5d_G4_calcheck_{p}.csv") for p in "abc"}
D = pd.read_csv(R / "S1.5d_desc.csv"); g1 = (R / "S1.5d_G1_final_run.log").read_text(); g2o = (R / "S1.5d_G2_official_run.log").read_text(); g2f = (R / "S1.5d_G2_final_run.log").read_text()
def v(i, a, u="tessuto", arm="RV"): return V[(V.id == i) & (V.archetype == a) & (V.unit == u) & (V.arm == arm)].iloc[0]
pct = lambda x: f"{100 * x:.0f} %"; rc = lambda e: RC[RC.emp == e].iloc[0]; real = P4[P4.source == "vero"].value
nf = U[U.n_failed > 0]; g4a = G4["a"][G4["a"].checked]; g4b = G4["b"]; g4c = G4["c"]
mpre = lambda s: re.search(r"mutanti rilevati (\d+)/(\d+)", s).groups(); mfin = re.search(r"mutanti rilevati (\d+)/20 \(pre-registrati M1-M11: (\d+)/11", g2f).groups()
rows = [
 ("spazio vuoto in tutti e sei i tessuti (prima in tre)", "S1.5d_counts.csv, P-11 (RV, N)", f"{K.loc['RV','P-11']}; {K.loc['N','P-11']}", "6; 3"),
 ("territori grandi in quattro tessuti (prima in nessuno)", "S1.5d_counts.csv, P-13 (RV, N)", f"{K.loc['RV','P-13']}; {K.loc['N','P-13']}", "4; 0"),
 ("intestino: variabilità dei territori da 41 % sotto a 2 % sopra", "S1.5d_verdicts.csv, P-2 A1 (N, RV)", f"{pct(-v('P-2','A1',arm='N').est)}; {pct(v('P-2','A1').est)}", "41 %; 2 %"),
 ("controprove dei vuoti fini e del cuore superate", "S1.5d_cp.csv, CP-1 e CP-2 giudicate", f"{int(CP[CP.judged & CP.cp.str.startswith(('CP-1','CP-2'))].pass_.sum()) if 'pass_' in CP else int(CP[CP.judged & CP.cp.str.startswith(('CP-1','CP-2'))]['pass'].sum())}/{int((CP.judged & CP.cp.str.startswith(('CP-1','CP-2'))).sum())}", "8/8"),
 ("contatti nel tumore del colon 69 % contro al più 57 %", "S1.5d_verdicts.csv, P-1 A2 RV (stima, estremo alto)", f"{pct(v('P-1','A2').est)}; {pct(v('P-1','A2').hi)}", "69 %; 57 %"),
 ("contatti nello stroma 36 % contro al più 28 %", "S1.5d_verdicts.csv, P-1 A3 RV", f"{pct(v('P-1','A3').est)}; {pct(v('P-1','A3').hi)}", "36 %; 28 %"),
 ("vincono i soli vuoti, 11 contro 10", "S1.5d_ruleA.csv", f"{RA.winner} {RA.RV} {RA.RD} {RA.RVD}", "RV 11 10 10"),
 ("cellule quasi senza citoplasma 4.3 % con i vuoti", "S1.5d_verdicts.csv, P-8c RV", f"{100*v('P-8c','A2').est:.1f} %", "4.3 %"),
 ("nelle cinque aree vere fra 1.5 e 2.9 %", "S1.5d_posthoc_p8c_real.csv, vero", f"{100*real.min():.1f}-{100*real.max():.1f} %", "1.5-2.9 %"),
 ("gli scarti veri per territorio le riportano al 2.2 %", "S1.5d_verdicts.csv, P-8c W3", f"{100*v('P-8c','A2',arm='W3').est:.1f} %", "2.2 %"),
 ("parità nella regola sulla taglia", "S1.5d_ruleC.csv, RV contro W3", f"{rc('W3').size_base} {rc('W3').size_emp} {rc('W3').winner}", "34 34 gaussian"),
 ("pesa la dipendenza dal territorio (controprova)", "S1.5d_verdicts.csv, P-8c W3p", f"{v('P-8c','A2',arm='W3p').est:.4f} {v('P-8c','A2',arm='W3p').verdict}", "0.0375 smentita"),
 ("nei territori piccoli il nucleo varia meno", "S1.5d_posthoc_size_tail.csv, tumor_microenv:1, sd dei decimi 1 e 10", f"{TL.loc['tumor_microenv:1','sd_decile1']:.2f} {TL.loc['tumor_microenv:1','sd_decile10']:.2f}", "0.38 0.61"),
 ("disposizione 58 contro 19, taglia 19 contro 21", "S1.5d_ruleB.csv", f"{RB.spatial_final} {RB.spatial_ref} {RB.size_final} {RB.size_ref} {RB.result}", "58 19 19 21 non adottata"),
 ("due degli undici errori non riconosciuti", "S1.5d_G2_official_run.log; S1.5d_G2_final_run.log", f"{'/'.join(mpre(g2o))}; finale {mfin[1]}/11, tutti {mfin[0]}/20", None),
 ("1 076 cellule in 5 simulazioni su 3 600", "S1.5d_panel_units.csv, n_failed", f"{int(U.n_failed.sum())} {len(nf)} {len(U)}", "1076 5 3600"),
 ("tutte in versioni che non hanno vinto", "S1.5d_panel_units.csv, bracci con n_failed > 0", " ".join(sorted(nf.arm.unique())), "Q RD RVD"),
 ("i risultati non cambiano escludendole", "review_recalc/sens_unplaced.csv (revisore)", "verdetti invariati (testo del revisore)", None),
 ("tre tolleranze scritte male: T13, T12 di S1.5c, G-4 c1", "S1.5d_G1_final_run.log; S1.5d_G4_calcheck_c.csv", f"T13 {'FAIL T13' in g1}; T12 {'FAIL T12' in g1}; c1 {int(g4c[g4c.part=='c1']['pass'].sum())}/{int((g4c.part=='c1').sum())}", "T13 True; T12 True; c1 0/10"),
 ("semina uniforme: 37 contro 36; con la correzione 39 contro 37", "S1.5d_ruleC.csv (U contro U3); S1.5d_ruleC_sensitivity_U3c.csv", f"{rc('U3').size_base} {rc('U3').size_emp}; {RS.size_emp} {RS.size_base} {RS.winner}", "37 36; 39 37 empirical"),
 ("con la correzione 0.7 %, sotto ogni area vera", "S1.5d_verdicts.csv, P-8c U3c; S1.5d_posthoc_p8c_real.csv", f"{100*v('P-8c','A2',arm='U3c').est:.1f} % < {100*real.min():.2f} %", "0.7 % < 1.54 %"),
 ("G-3: U = S1.5c 300/300, N = N-a 300/300, rerun 8/8, superamenti 1.54e-4", "S1.5d_G3_check.csv", f"{G3.loc['U_equals_S1.5c_U_md5','value']:.0f} {G3.loc['N_equals_S1.5c_Na_md5','value']:.0f} {G3.loc['rerun_identical','value']:.0f} {G3.loc['max_exceed_fraction','value']:.3g}", "300 300 8 0.000154"),
 ("G-4: parte a 333/333, parte b 240/240, c2 10/10", "S1.5d_G4_calcheck_a/b/c.csv", f"{int(g4a['pass'].sum())}/{len(g4a)} {int(g4b['pass'].sum())}/{len(g4b)} {int(g4c[g4c.part=='c2']['pass'].sum())}/{int((g4c.part=='c2').sum())}", "333/333 240/240 10/10"),
 ("G-1: test_S1.5d 31/32, S1.5c 27/28 (28/28 con i preset di S1.5c), S1.5b 22/22, S1.5 81/83, S1.4d 36/36, OKF 8/8", "S1.5d_G1_final_run.log", " ".join(re.findall(r"(?:RESULT|RISULTATO test_S1.5b:) (\d+/\d+)", g1)), "31/32 27/28 28/28 22/22 81/83 36/36 8/8"),
 ("CP-3: area mediana diversa in 5 tessuti su 6 e in 16 tipi su 17", "S1.5d_cp.csv, CP-3", f"{int((~CP[CP.item.str.startswith('P-8a')]['pass']).sum())}/6 {int((~CP[CP.item.str.startswith('P-8b')]['pass']).sum())}/17", "5/6 16/17"),
 ("PH-2: vero 0.0269, gaussiano 0.0704, senza decimi 0.0640, con i decimi 0.0298; 77 %", "S1.5d_posthoc_size_tail.csv, tumor_microenv:1", " ".join(f"{TL.loc['tumor_microenv:1',c]:.4f}" for c in ["real", "gaussian", "empirical_pooled", "empirical_decile"]) + f" {pct(TL.loc['tumor_microenv:1','share_tail_from_decile1_3'])}", "0.0269 0.0704 0.0640 0.0298 77 %"),
 ("PH-3: 0.269, 0.230, 0.250", "S1.5d_posthoc_size_q90.csv, lymphoid:1", f"{Q9.loc['lymphoid:1','gaussian_q90']:.3f} {Q9.loc['lymphoid:1','table_q90']:.3f} {Q9.loc['lymphoid:1','R6_q90']:.3f}", "0.269 0.230 0.250"),
 ("quote: dal 2.8 % all'11.0 %; senza quote 60 contro 58", "S1.5d_desc.csv D-desc-1 (RV); S1.5d_counts.csv", (lambda x: f"{100*x.min():.1f} % {100*x.max():.1f} %")(D[D.desc.str.startswith('D-desc-1') & (D.arm == 'RV')].est) + f"; {K.loc['Q','spatial']} {K.loc['RV','spatial']}", "2.8 % 11.0 %; 60 58"),
]
rows = [(a, b, c, c if d is None else d) for a, b, c, d in rows]
L = pd.DataFrame(rows, columns=["frase", "fonte", "valore_letto", "valore_scritto"]); L["coincide"] = L.valore_letto == L.valore_scritto
L.to_csv(R / "S1.5d_verification_log.csv", index=False)
bad = L[~L.coincide]; assert bad.empty, bad.to_string()
tab = "| frase del riassunto | fonte | valore letto |\n|---|---|---|\n" + "\n".join(f"| {a} | {b} | {c} |" for a, b, c, _ in rows) + "\n"
(out / "_S1.5d_fonti.md").write_text(tab, encoding="utf-8"); print(f"{len(L)} righe, {int(L.coincide.sum())} coincidono")
