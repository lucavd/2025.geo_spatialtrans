#!/usr/bin/env python3
"""tools/S1.5_report_parts.py - estrae dal file dello step (docs/steps/S1.5.md) le sezioni incluse nel report (il testo esiste
una volta sola) e costruisce la tabella «Da dove vengono i numeri del riassunto» rileggendo ogni valore dai file dei risultati.
Uso (root del repo): .venv/bin/python -B tools/S1.5_report_parts.py"""
import pathlib, re, json
import pandas as pd
R = pathlib.Path("results/S1.5"); out = pathlib.Path("reports"); T = "/mnt/micron/geo_spatialtrans/S1.5/tests/final2"
t = pathlib.Path("docs/steps/S1.5.md").read_text(encoding="utf-8").split("\n---\n", 1)[1]
def sec(start, stop_re):
    i = t.index(start); m = re.search(stop_re, t[i + len(start):]); j = i + len(start) + (m.start() if m else len(t))
    return t[i:j].strip() + "\n"
(out / "_S1.5_inbreve.md").write_text(sec("\n# In breve\n", r"\n# "), encoding="utf-8")
(out / "_S1.5_rinviati.md").write_text(sec("\n# Dubbi rinviati\n", r"\n# ").replace("# Dubbi rinviati\n", "", 1), encoding="utf-8")
(out / "_S1.5_verifiche.md").write_text(sec("\n# Verifiche di Luca\n", r"\n# ").replace("# Verifiche di Luca\n", "", 1), encoding="utf-8")
h_rev = "\n## Revisione avversariale e correzioni"; h_dopo = "\n## Dopo il risultato"; h_dev = "\n## Note di sviluppo"
rev = sec(h_rev, r"\n#{1,2} "); rev = rev.split("\n", 1)[1]
(out / "_S1.5_review.md").write_text(rev, encoding="utf-8")
(out / "_S1.5_dopo.md").write_text(sec(h_dopo, r"\n#{1,2} ").split("\n", 1)[1], encoding="utf-8")
(out / "_S1.5_sviluppo.md").write_text(sec(h_dev, r"\n#{1,2} ").split("\n", 1)[1], encoding="utf-8")
# ---- fonti dei numeri del riassunto
v = pd.read_csv(R / "S1.5_verdicts.csv"); w = v[v.winner_arm]
g = lambda i, a, q=None: w[(w.id == i) & (w.archetype == a) & ((w.quantity.str.startswith(q)) if q else True)].iloc[0]
log = pathlib.Path(T, "test_S1.5.log").read_text(encoding="utf-8")
e10d = re.search(r"PASS E10d .* ([0-9.]+) vs ([0-9.]+) um", log); res = re.findall(r"RESULT (\d+)/(\d+) PASS", log)[-1]
s1 = re.search(r"PASS S1 .* (\d+)\s*$", log, re.M)
pu = pd.read_csv(R / "S1.5_panel_units.csv"); rc = pd.read_csv(R / "S1.5_G3_recheck.csv"); rp = pd.read_csv(R / "S1.5_G3_repro.csv")
fv = pd.read_csv(R / "S1.5_G4_faults_v2.csv"); nv = pd.read_csv(R / "S1.5_G4_nulls_v2.csv")
rf = pd.read_csv(R / "S1.5_review_faults_after.csv"); rf = rf[(rf.fault != "none") & rf.detected.notna()]
rm = pd.read_csv(R / "S1.5_review_mutants_after.csv"); rmR = rm[rm.id.str.startswith("R")]
sr = pd.read_csv(R / "S1.5_size_rule.csv").set_index("arm")
cal = pd.read_csv("results/S1.4/S1.4_calibration.csv"); cal = cal[cal.method == "cellpose_rgb"].set_index("archetype").density_mm2
r3 = pd.read_csv("results/R3b/R3b_refs.csv"); r3 = r3[r3.metric == "density_mm2"].set_index("archetype")
bio = pathlib.Path("docs/BIO_REFERENCES.md").read_text(encoding="utf-8")
man = {a: float(re.search(rf"^\| {b} \| {a} \|[^|]*\| ([0-9.]+) \|", bio, re.M).group(1)) for a, b in zip(["A1","A2","A3","A4","A5","A6"], ["B-043","B-044","B-045","B-046","B-047","B-048"])}
ph = pd.read_csv(R / "S1.5_posthoc_units.csv"); phr = ph[ph.arm == "realgen"].groupby(["archetype", "roi_id"])[["median_nc", "frac_nc095"]].mean().groupby("archetype").mean()
pw = pd.read_csv(R / "S1.5_posthoc_windows.csv"); pws = pw.groupby(["arm", "archetype", "seed"]).s1.mean().groupby(["arm", "archetype"]).mean()
cp = pd.read_csv(R / "S1.5_CP1_verdicts.csv"); cp1 = cp[(cp.variant == "pipeline") & (cp.scale_true == 1)].set_index("archetype")
cs = pd.read_csv(R / "S1.5_CP1_scale.csv"); csb = pd.read_csv(R / "S1.5_CP1_scale_est_b.csv")
wrong = lambda d: int(((~d.refused) & ((d.ratio - 1).abs() > 0.2) & (d.variant == "pipeline")).sum())
npipe = int((cs.variant == "pipeline").sum())
sens = [c for c in json.load(open(R / "review/review_recalc.json"))["structured_output"]["checks"] if c["id"] == "SENS-P6-A6"][0]
p8 = w[w.id == "P-8"]; tim = p8[p8.quantity.str.startswith("tempo")].set_index("archetype"); mem = p8[p8.quantity.str.startswith("memoria")].set_index("archetype")
f = lambda x, d=3: f"{x:.{d}f}"
rows = [
 ("raggio mediano 4,86 µm con la scala giusta e con quella 10 volte più grande", "tests/final2/test_S1.5.log, riga E10d", f"{e10d.group(1)} vs {e10d.group(2)} µm"),
 ("le altre 81 passano", "tests/final2/test_S1.5.log, RESULT", f"{res[0]}/{res[1]} (FAIL: E10, E10b)"),
 ("metadati esatti in tutti e 5 i dataset", "tests/final2/test_S1.5.log, riga S1", s1.group(1) if s1 else "manca"),
 ("600 simulazioni senza errori", "S1.5_panel_units.csv, n_error", f"{int((pu.n_error == 0).sum())}/{len(pu)}"),
 ("rifatte con lo stesso seme danno file identici", "S1.5_G3_repro.csv; S1.5_G3_recheck.csv md5 = md5_panel", f"{int(rp.md5_rerun_identical.sum())}/{len(rp)}; {int((rc.md5 == rc.md5_panel).sum())}/{len(rc)}"),
 ("tutti i 630 guasti rilevati", "S1.5_G4_faults_v2.csv, status_expected = ERROR", f"{int((fv.status_expected == 'ERROR').sum())}/{len(fv)}"),
 ("rilevati tutti i guasti costruiti dal revisore", "S1.5_review_faults_after.csv, detected_any", f"{int(rf.detected_any.sum())}/{len(rf)}"),
 ("18 delle 19 modifiche sabotate del revisore", "S1.5_review_mutants_after.csv, id R01-R19", f"{int(rmR.detected.sum())}/{len(rmR)}; non rilevato: {', '.join(rmR.id[~rmR.detected])}"),
 ("stessi esiti dei controlli nelle 600 simulazioni", "S1.5_G3_recheck.csv statuses = S1.5_panel_units.csv statuses", f"{int((rc.merge(pu, on=['archetype','roi_id','seed','arm']).pipe(lambda m: (m.statuses_x == m.statuses_y).sum())))}/{len(rc)}"),
 ("linfonodo: 5 minuti e 10 secondi, 8,5 GB", "S1.5_verdicts.csv P-8 A4", f"{f(tim.loc['A4'].estimate, 1)} s; {f(mem.loc['A4'].estimate, 2)} GB"),
 ("altri tessuti da 18 secondi a poco più di 2 minuti", "S1.5_verdicts.csv P-8 tempo", f"{f(tim.drop('A4').estimate.min(), 1)}-{f(tim.drop('A4').estimate.max(), 1)} s"),
 ("area mediana del nucleo dentro la variabilità in sei tessuti; 13, 31, 64 µm²", "S1.5_verdicts.csv P-1", f"{int((w[w.id == 'P-1'].outcome == 'confermata').sum())}/6; {f(g('P-1','A1').estimate,2)}, {f(g('P-1','A2').estimate,2)}, {f(g('P-1','A5').estimate,2)}"),
 ("sei tessuti su sei contro quattro", "S1.5_size_rule.csv", f"preserve_area {sr.loc['preserve_area'].P1_confirmed}, none {sr.loc['none'].P1_confirmed}"),
 ("densità contate a mano diverse da quelle del modello", "BIO_REFERENCES B-043...B-048; S1.4_calibration.csv (cellpose_rgb), R3b_refs.csv (A4, A6)", "manuale / modello: " + "; ".join(f"{a} {man[a]:.0f} / {(cal[a] if a in ['A1','A2','A3','A5'] else r3.loc[a].consensus_sdarea):.0f}" for a in ["A1","A2","A3","A4","A5","A6"])),
 ("frazione del nucleo giusta in quattro tessuti", "S1.5_verdicts.csv P-2", f"{int((w[w.id == 'P-2'].outcome == 'confermata').sum())}/6"),
 ("intestino 0,18 contro 0,22-0,23", "S1.5_verdicts.csv P-2 A1", f"{f(g('P-2','A1').estimate)} [{f(g('P-2','A1').int_lo,4)}, {f(g('P-2','A1').int_hi,4)}]"),
 ("tumore del colon 0,306 contro almeno 0,309", "S1.5_verdicts.csv P-2 A2", f"{f(g('P-2','A2').estimate)} ± {f(2*g('P-2','A2').se,4)}; limite {f(g('P-2','A2').int_lo,4)}"),
 ("con i territori delle cellule vere torna (0,32)", "S1.5_posthoc_units.csv realgen A2 median_nc", f(phr.loc["A2"].median_nc)),
 ("3,6 % quasi senza citoplasma contro al massimo 2,9 %", "S1.5_verdicts.csv P-4 A2", f"{f(g('P-4','A2').estimate,4)}; massimo reale {f(g('P-4','A2').int_hi,4)}"),
 ("tumore del colon 62 su 100 contro 22-57", "S1.5_verdicts.csv P-5 A2", f"{f(g('P-5','A2').estimate)} [{g('P-5','A2').int_lo}, {g('P-5','A2').int_hi}]"),
 ("linfonodo 83 contro 35-75", "S1.5_verdicts.csv P-5 A4", f"{f(g('P-5','A4').estimate)} [{g('P-5','A4').int_lo}, {g('P-5','A4').int_hi}]"),
 ("intestino e cervello tornano; stroma e cuore al limite", "S1.5_verdicts.csv P-5", ", ".join(f"{a} {g('P-5',a).outcome}" for a in ["A1","A3","A5","A6"])),
 ("con i centri veri dentro in cinque tessuti, al limite nel linfonodo", "S1.5_posthoc_windows.csv realgen (media dei semi)", ", ".join(f"{a} {pws.loc[('realgen', a)]:.3f}" for a in ["A1","A2","A3","A4","A5","A6"])),
 ("variabilità delle aree più bassa di un quarto fino a quasi la metà, in cinque tessuti", "S1.5_verdicts.csv P-6", ", ".join(f"{a} {g('P-6',a).estimate:+.3f}" for a in ["A1","A2","A3","A4","A5"])),
 ("nel cuore torna, ma è fragile", "S1.5_verdicts.csv P-6 A6; review/review_recalc.json SENS-P6-A6", f"{g('P-6','A6').estimate:+.3f} [{g('P-6','A6').lo2:+.4f}, {g('P-6','A6').hi2:+.4f}]; {sens['reviewer_value']}"),
 ("a 1 micron per pixel la scala giusta solo nel linfonodo; nel cuore sbaglia di 4 volte", "S1.5_CP1_verdicts.csv pipeline, 1 µm/px", "; ".join(f"{a} {cp1.loc[a].outcome}" for a in ["A1","A2","A3","A4","A5","A6"]) + f"; cuore {cp1.loc['A6'].ratio_mean:.2f}"),
 ("23 casi su 120 un numero sbagliato senza avvisare", "S1.5_CP1_scale.csv pipeline, non rifiutati con abs(rapporto − 1) > 0.2", f"{wrong(cs)}/{npipe}"),
 ("un controllo aggiunto dopo li riduce a 8", "S1.5_CP1_scale_est_b.csv, stessa regola", f"{wrong(csb)}/{npipe}"),
]
tab = "## Da dove vengono i numeri del riassunto\n\n| frase del riassunto | fonte (file, riga) | valore letto |\n|---|---|---|\n" + "".join(f"| {a} | {b} | {c} |\n" for a, b, c in rows)
(out / "_S1.5_fonti.md").write_text(tab, encoding="utf-8"); print(tab)
