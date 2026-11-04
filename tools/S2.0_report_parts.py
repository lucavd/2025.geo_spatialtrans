#!/usr/bin/env python3
"""tools/S2.0_report_parts.py - estrae dal file dello step (docs/steps/S2.0.md) le sezioni incluse nel report, costruisce il riepilogo della
revisione avversariale dai CSV dei revisori e la tabella «Da dove vengono i numeri del riassunto» rileggendo ogni valore dai file dei risultati.
Uso (root del repo): .venv/bin/python -B tools/S2.0_report_parts.py"""
import pathlib, re, json, glob
import numpy as np, pandas as pd
R = pathlib.Path("results/S2.0"); RW = R / "wfix"; out = pathlib.Path("reports"); RUN = "/mnt/micron/geo_spatialtrans/S2.0/run1"; REV = "/mnt/micron/geo_spatialtrans/S2.0/review"
t = pathlib.Path("docs/steps/S2.0.md").read_text(encoding="utf-8").split("\n---\n", 1)[1]
def sec(start, stop_re):
    i = t.index(start); m = re.search(stop_re, t[i + len(start):]); j = i + len(start) + (m.start() if m else len(t))
    return t[i:j].strip() + "\n"
(out / "_S2.0_inbreve.md").write_text(sec("\n# In breve\n", r"\n# ").replace("# In breve\n", "", 1), encoding="utf-8")
(out / "_S2.0_rinviati.md").write_text(sec("\n# Dubbi rinviati\n", r"\n# ").replace("# Dubbi rinviati\n", "", 1), encoding="utf-8")
(out / "_S2.0_verifiche.md").write_text(sec("\n# Verifiche di Luca\n", r"\n# ").replace("# Verifiche di Luca\n", "", 1), encoding="utf-8")
# ---- revisione avversariale
def passes(f):
    d = pd.read_csv(f); p = d["pass"].astype(str).str.strip().str.lower().isin(["true", "1"]); return int(p.sum()), len(d), d[~p]
rev = [("documento", "documento/S2.0_G5_documento_checks.csv", "(e) documento contro dati e decisioni: ogni numero, fonte e decisione"),
       ("nucleo", "nucleo/S2.0_G5_nucleo_comparisons.csv", "(c) controllo interno in 3 ROI (tumore, cervello), (d) ψ osservato e simulato, PH-2, PH-3 in 7 ROI"),
       ("strisce", "strisce/S2.0_G5_strisce_comparisons.csv", "(a) M1 dai dati grezzi in intestino, linfonodo, cervello; 5 prove di falsificazione"),
       ("bordi", "bordi/S2.0_G5_bordi_checks.csv", "(b) quattro regole rivalutate su 40 finestre, ottimizzazione propria su 20, R-A da tutte le 314")]
md = ["# Revisione avversariale (G-5)\n", "Quattro revisori indipendenti, ciascuno con codice proprio e accesso in sola lettura; rapporti e codice in `" + REV + "/<revisore>/`.\n",
      "| revisore | che cosa | controlli superati | non superati |", "|---|---|---|---|"]
for nm, f, what in rev:
    k, n, bad = passes(f"{REV}/{f}"); md.append(f"| {nm} | {what} | {k}/{n} | {', '.join(bad['id'].astype(str).head(8))}{' …' if len(bad) > 8 else ''} |")
nd = len(re.findall(r"^- \*\*D-", open(f"{REV}/documento/S2.0_G5_documento.md").read(), re.M))
md += ["", f"Dubbi del revisore del documento: {nd}. Come sono stati chiusi è scritto qui sotto e nel file dello step.", "",
       "- **Strisce**: righe e colonne con tessuto e zero UMI davano un fattore sentinella che falsava la parità (linfonodo, cervello); la parte lenta dei fattori deriva con le iterazioni. Chiusi con PH-5: righe a zero UMI tolte, per il simulatore la sola parte fine (stabile), passo 1 calibrato su sezioni sintetiche. Gli esiti di P-1 e P-2 non cambiano.",
       "- **Bordi**: soluzioni ferme con la coda spenta (fino a 32 unità di log-verosimiglianza nel campione del revisore). Chiuso con PH-6: 107 coppie finestra × regola riadattate, nessun esito pre-registrato cambia, cambiano lunghezza e massa della coda; i valori del design sono quelli corretti. La tolleranza di 1e-6 sulla log-verosimiglianza era irraggiungibile (modello in precisione singola): giudicata a 0.05, come in R5.",
       "- **Nucleo**: 203/203; l'eccesso di variabilità c'è anche nella corona e nello spazio acellulare (provvisorio): il campo del simulatore copre tutto il tessuto (DEC-040). Nell'intestino v è sottostimato di circa il 20 % (bin profondi sulla maschera non spostata): BACKLOG.",
       "- **Documento**: i numeri sbagliati e le frasi oltre i dati sono corretti nel design v0.2; le voci che restano aperte sono nei limiti (sezione 17) e nel BACKLOG.", ""]
(out / "_S2.0_review.md").write_text("\n".join(md) + "\n", encoding="utf-8")
# ---- fonti dei numeri del riassunto
TT = pd.read_csv(RW / "S2.0_edges_tiles.csv"); TT = TT[~TT.error.astype(bool)]
E = pd.read_csv(RW / "S2.0_edges_sections.csv"); E4 = E[E.key.isin(["A1", "A23", "A4", "A5"])]
def ev(rule, p): x = E4[(E4.rule == rule) & (E4.param == p)]["median"]; return x.min(), x.max()
L = pd.read_csv(RW / "S2.0_CP1_lumens.csv"); J = L[L.judged_CP1 & (L.rule == "C2")]
I = pd.read_csv(RW / "S2.0_inner_summary.csv"); I2 = I[I.rule == "C2"].set_index("archetype")
P5 = pd.read_csv(R / "S2.0_PH5_sections.csv"); P5 = P5[P5.n_iter == 10]
ST4 = json.load(open(f"{RUN}/stripes/A4.json")); PH = pd.read_csv(R / "S2.0_PH23_disp.csv").set_index("archetype")
share = 1 - (PH.ratio - 1) / (PH.ratio_pre - 1)
wf = json.load(open(f"{RUN}/edges_wfix/A23_t2700_1800.json"))
ell_max = E[(E.rule == "C2") & (E.param == "ell") & (E.key != "A6")].sort_values("median").iloc[-1]
f = lambda x, d=2: f"{x:.{d}f}"
B1 = pd.read_csv(f"{REV}/bordi/S2.0_G5_bordi_checks.csv"); B1MAX = B1[B1.id.str.match(r"b1-A\d")].abs_diff.astype(float).max()
SEP = {k: sorted(TT[TT.key == k][r + "_aic"].sum() for r in ("C1", "C2", "C3")) for k in ("A1", "A23", "A4", "A5")}
SEPKEY = min(SEP, key=lambda k: SEP[k][1] - SEP[k][0]); SEPMIN = (SEP[SEPKEY][1] - SEP[SEPKEY][0]) / 2
rows = [
 ("zona di bordo: entro 7-25 µm dal bordo, allargamento di 4.5-6.2 µm", "results/S2.0/wfix/S2.0_edges_sections.csv, C2, D e sigma_e, sezioni giudicate", f"D {f(ev('C2','D')[0])}-{f(ev('C2','D')[1])} µm; σ_e {f(ev('C2','sigma_e')[0])}-{f(ev('C2','sigma_e')[1])} µm"),
 ("vince in tutte e quattro le sezioni giudicate", "results/S2.0/wfix/S2.0_RA.json (somma degli AIC, soluzioni corrette)", json.dumps({k: v["winner"] for k, v in json.load(open(RW / "S2.0_RA.json"))["per_section"].items()})),
 ("bordi dei buchi interni: nessuna delle 8 fuori (6 dentro, 2 incerte)", "results/S2.0/wfix/S2.0_CP1_lumens.csv, C2, celle giudicate", ", ".join(f"{v}: {n}" for v, n in J.verdict.value_counts().items())),
 ("nucleo/citoplasma riprodotto in tumore, stroma, cervello; incerto in intestino e linfonodo; non nel cuore", "results/S2.0/wfix/S2.0_inner_summary.csv, C2", "; ".join(f"{a} {f(I2.loc[a,'ratio_mean'],3)} {I2.loc[a,'verdict']}" for a in ["A1", "A2", "A3", "A4", "A5", "A6"])),
 ("coda fino a 41.5 µm nel linfonodo", "results/S2.0/wfix/S2.0_edges_sections.csv, C2, ell", f"{ell_max.key} {f(ell_max['median'])} µm"),
 ("il difetto c'era già in R5 in una finestra del tumore", "run1/edges_wfix/A23_t2700_1800.json, C0 (nucleo di R5)", f"miglioramento {f(wf['C0']['improvement'])} unità"),
 ("alternanza fino a 1.75 volte fra colonne dispari e pari nell'intestino", "results/S2.0/S2.0_PH5_sections.csv, A1 cols", f"{f(P5[(P5.key=='A1')&(P5.axis=='cols')].odd_over_even.iloc[0],3)}"),
 ("ogni riga e colonna cattura circa il 20 % in più o in meno", "results/S2.0/S2.0_PH5_sections.csv, sd_repro_hp (log)", f"{f(P5.sd_repro_hp.min(),3)}-{f(P5.sd_repro_hp.max(),3)}"),
 ("alternanza del linfonodo da 0.98 a 1.00", "run1/stripes/A4.json (rows) e S2.0_PH5_sections.csv (A4 rows)", f"{f(ST4['rows']['odd_over_even'],4)} → {f(P5[(P5.key=='A4')&(P5.axis=='rows')].odd_over_even.iloc[0],4)}"),
 ("variabilità da 3 a 11 volte quella simulata", "results/S2.0/S2.0_PH23_disp.csv, ratio_pre", f"{f(PH.ratio_pre.min())}-{f(PH.ratio_pre.max())}"),
 ("da un quarto (cuore) a sette decimi (cervello) dalle strisce, circa metà negli altri", "results/S2.0/S2.0_PH23_disp.csv, 1 − (ratio − 1)/(ratio_pre − 1)", "; ".join(f"{a} {f(share[a])}" for a in share.index)),
 ("quadratini vicini si somigliano, a 4 µm molto meno", "results/S2.0/S2.0_PH23_disp.csv, exc_lag1 / exc_lag2", f"2 µm {f(PH.exc_lag1.min())}-{f(PH.exc_lag1.max())}; 4 µm {f(PH.exc_lag2.min())}-{f(PH.exc_lag2.max())}"),
 ("scarti fino a 2 centesimi di unità (precisione singola)", f"{REV}/bordi/S2.0_G5_bordi_checks.csv, righe b1-<sezione>-*, abs_diff massimo", f"{B1MAX:.4f}"),
 ("la somma sulle finestre di una sezione separa le regole di oltre cento unità", "results/S2.0/wfix/S2.0_edges_tiles.csv: differenza minima fra le somme degli AIC della prima e della seconda regola, divisa per 2 (unità di log-verosimiglianza)", f"{SEPMIN:.1f} ({SEPKEY})"),
 ("il revisore del documento ha trovato 38 punti", f"{REV}/documento/S2.0_G5_documento.md, voci D-", str(nd)),
]
for nm, fcsv, _ in rev:
    k, n, _b = passes(f"{REV}/{fcsv}"); rows.append((f"revisore {nm}", f"{REV}/{fcsv}, colonna pass", f"{k}/{n}"))
for g in ("S2.0_G1.txt", "S2.0_G2.txt", "S2.0_G3.txt", "S2.0_G4.txt"):
    p = R / g; rows.append((f"cancello {g[5:7]}", f"results/S2.0/{g}, ultima riga", p.read_text().strip().split("\n")[-1] if p.exists() else "mancante"))
fm = "| frase del riassunto | fonte (tabella o file, riga) | valore letto |\n|---|---|---|\n" + "\n".join(f"| {a} | {b} | {c} |" for a, b, c in rows) + "\n"
(out / "_S2.0_fonti.md").write_text(fm, encoding="utf-8"); print(fm)
