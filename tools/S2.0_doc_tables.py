"""tools/S2.0_doc_tables.py — S2.0: tabella dei parametri del design dello Step 2 e righe nuove di BIO_REFERENCES, generate dai
risultati (nessun numero ricopiato a mano). Uscite: results/S2.0/S2.0_param_table.md, results/S2.0/S2.0_bio_rows.md.
Uso: .venv/bin/python tools/S2.0_doc_tables.py <outdir_run1>"""
import sys, json, ast
import numpy as np, pandas as pd
OUT = sys.argv[1]; RES = "results/S2.0"
SEC = ["A1", "A23", "A4", "A5", "A6"]; ARCH = ["A1", "A2", "A3", "A4", "A5", "A6"]; SOF = {"A1": "A1", "A2": "A23", "A3": "A23", "A4": "A4", "A5": "A5", "A6": "A6"}
DS = {"A1": "Visium_HD_Mouse_Small_Intestine", "A23": "Visium_HD_Human_Colon_Cancer", "A4": "Visium_HD_Human_Lymph_Node_FFPE", "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}
BG = {"A1": ("B-342", 0.161), "A23": ("B-348", 0.197), "A4": ("B-354", 0.042), "A5": ("B-360", 0.309), "A6": ("B-366", 0.117)}   # R5, gia' in BIO_REFERENCES
E = pd.read_csv(f"{RES}/S2.0_edges_sections.csv"); I = pd.read_csv(f"{RES}/S2.0_inner_summary.csv"); IR = pd.read_csv(f"{RES}/S2.0_inner_roi.csv")
P = pd.read_csv(f"{RES}/S2.0_PH23_disp.csv").set_index("archetype"); L = pd.read_csv(f"{RES}/S2.0_CP1_lumens.csv")
ST = {k: json.load(open(f"{OUT}/stripes/{k}.json")) for k in SEC}
def e(k, rl, p): x = E[(E.key == k) & (E.rule == rl) & (E.param == p)].iloc[0]; return x["median"], x.lo, x.hi
HEART = {p: float(np.median([e(k, "C2", p)[0] for k in ("A1", "A23", "A4", "A5")])) for p in ("sigma_e", "D", "w", "ell")}
f2 = lambda x: f"{x:.2f}"; f3 = lambda x: f"{x:.3f}"
bio, nb = [], 374
def add(arch, metr, val, iv, unit, src, note):
    global nb; bio.append(f"| B-{nb} | {arch} | {metr} | {val} | {iv} | {unit} | {src} | S2.0 | {note} |"); nb += 1; return f"B-{nb-1}"
SRC_ST = "tools/S2.0_stripes.py@S2.0 (results/S2.0/S2.0_predictions.csv, run1/stripes/<sezione>.json)"
bst = {}
for k in SEC:
    r, c_ = ST[k]["rows"], ST[k]["cols"]; t = ST[k]["tiles"]
    arch = "A2, A3" if k == "A23" else k
    bst[k] = add(arch, "strisce a 2 µm: parità dispari/pari (righe; colonne), SD riproducibile dei fattori non periodici per riga e per colonna (log), correlazione a passo 1",
                 f"parità {f3(r['odd_over_even'])}; {f3(c_['odd_over_even'])}. SD {f3(r['sd_repro'])}; {f3(c_['sd_repro'])}",
                 f"SD righe [{f3(r['sd_lo'])}, {f3(r['sd_hi'])}], colonne [{f3(c_['sd_lo'])}, {f3(c_['sd_hi'])}]; passo 1 {f3(r['lag1_repro'])}; {f3(c_['lag1_repro'])}; r fra metà {f2(r['corr_halves'])}; {f2(c_['corr_halves'])}; SD fra tessere della parità {f3(t['rows']['sd_between_repro'])}; {f3(t['cols']['sd_between_repro'])}",
                 "rapporto; log", f"{SRC_ST} + {DS[k]}", "P-1 smentita (non periodica presente); P-2 " + ("smentita nelle righe e indeterminata nelle colonne: parità variabile nello spazio" if k == "A4" else "confermata") + "; PH-1: i fattori si ripetono sul vetro lontano (cattura tecnica per riga e colonna)")
bc2 = {}
for k in SEC:
    vals = {p: e(k, "C2", p) for p in ("sigma_e", "D", "w", "ell")}
    arch = "A2, A3" if k == "A23" else k
    nt = int(E[(E.key == k) & (E.rule == "C2") & (E.param == "sigma_e")].n_tiles.iloc[0])
    dA = E[(E.key == k) & (E.rule == "C2") & (E.param == "dAIC_vs_C0_per_window")]["median"].iloc[0]
    note = "regola d'uscita scelta da R-A (zona di bordo, vince in 4 sezioni su 4 sulla somma degli AIC)" if k != "A6" else "5 finestre: D non identificata (intervallo fino al limite 10 000 µm); nel preset si usa la mediana delle altre 4 sezioni (DEC-037)"
    bc2[k] = add(arch, "spostamento ai bordi, regola C2 (zona di bordo): σ_e, profondità D, massa della coda w, lunghezza della coda ℓ (mediana delle finestre)",
                 f"σ_e {f2(vals['sigma_e'][0])} µm; D {f2(vals['D'][0])} µm; w {f3(vals['w'][0])}; ℓ {f2(vals['ell'][0])} µm",
                 "; ".join(f"{n} [{f2(v[1])}, {f2(v[2])}]" for n, v in (("σ_e", vals["sigma_e"]), ("D", vals["D"]), ("w", vals["w"]), ("ℓ", vals["ell"]))) + f"; {nt} finestre; ΔAIC mediano per finestra contro C0 {dA:.1f}",
                 "µm; frazione", f"tools/S2.0_edges.py, S2.0_analysis.py@S2.0 (results/S2.0/S2.0_edges_sections.csv, S2.0_RA.json) + {DS[k]}", note)
br, bv = {}, {}
for A_ in ARCH:
    x = I[(I.archetype == A_) & (I.rule == "C2")].iloc[0]; rr = IR[(IR.archetype == A_) & (IR.rule == "C2")].r.dropna()
    br[A_] = add(A_, "rapporto nucleo/citoplasma prima dello spostamento con la regola C2 (invertito sulla statistica di R4 PH-1) e controllo interno (previsto/osservato senza i bin entro 2 µm dal bordo del nucleo)",
                 f"r {f2(x.r_median)} (mediana di 5 ROI)", f"ROI {f2(rr.min())}-{f2(rr.max())}; controllo interno {f3(x.ratio_mean)} [{f3(x.lo)}, {f3(x.hi)}], {x.verdict}", "adimensionale",
                 f"tools/S2.0_inner.py, S2.0_analysis.py@S2.0 (results/S2.0/S2.0_inner_summary.csv) + 5 ROI di R2, nuclei di R4, {DS[SOF[A_]]}",
                 "il modello in avanti non descrive l'interno del nucleo: valore da non usare senza riserva" if x.verdict == "smentita" else ("da confermare" if x.verdict != "confermata" else "regola C2, σ interno 0.7 µm"))
    p = P.loc[A_]
    bv[A_] = add(A_, "variabilità dentro il nucleo oltre il caso a 2 µm: varianza moltiplicativa v in più (dopo strisce e spostamento) e correlazione in più dei residui fra bin adiacenti e a 4 µm",
                 f"v {f3(p.v)}", f"[{f3(p.v_lo)}, {f3(p.v_hi)}]; senza i fattori di riga e colonna {f3(p.v_pre)}; correlazione in più a 2 µm {f3(p.exc_lag1)} [{f3(p.exc_lag1_lo)}, {f3(p.exc_lag1_hi)}], a 4 µm {f3(p.exc_lag2)} [{f3(p.exc_lag2_lo)}, {f3(p.exc_lag2_hi)}]; dispersione oss/sim {f2(p.ratio)}",
                 "adimensionale", f"tools/S2.0_disp.py, S2.0_ph23.py@S2.0 (results/S2.0/S2.0_PH23_disp.csv) + 5 ROI di R2, {DS[SOF[A_]]}",
                 "post hoc PH-2/PH-3 (P-4 smentita per eccesso come scritta); struttura corta (1-2 µm), biologica o tecnica non distinguibile" + ("; controllo interno del cuore smentito" if A_ == "A6" else ""))
j = L[L.judged_CP1 & (L.rule == "C2")]
bl = add("A1, A2-A3, A4, A6", "lumi interni: UMI osservati/previsti con la regola C2 a 2-4 e 4-8 µm dal bordo (CP-1)",
         "; ".join(f"{k} {b} {f2(o)}" for k, b, o in zip(j.key, j.band, j.obs_over_pred)), "0 smentite su 8 rispetto a [0.8, 1.25]; " + "; ".join(f"{k} {b} [{f2(lo)}, {f2(hi)}]" for k, b, lo, hi in zip(j.key, j.band, j.lo, j.hi)),
         "rapporto", "tools/S2.0_lumens.py, S2.0_analysis.py@S2.0 (results/S2.0/S2.0_CP1_lumens.csv) + lumi di R5 CP-2", "bootstrap per tessere come R5; C0 riproduce R5 CP-2 esattamente; cervello (2 lumi) descrittivo")
open(f"{RES}/S2.0_bio_rows.md", "w").write("\n".join(bio) + "\n")
# tabella dei parametri
cols = ["intestino (A1)", "tumore (A2)", "stroma (A3)", "linfonodo (A4)", "cervello (A5)", "cuore (A6)"]
def row(name, vals, src): return f"| {name} | " + " | ".join(vals) + f" | {src} |"
T = ["| parametro | " + " | ".join(cols) + " | fonte |", "|---|" + "---|" * 6 + "---|"]
T.append(row("trascritti per cellula", ["per tipo, mediana e dispersione di R4"] * 6, "B-275...B-335; [DEC-035](/decisions/DEC-035.md)"))
T.append(row("rapporto nucleo/citoplasma r (regola C2)", [f2(I[(I.archetype == A_) & (I.rule == 'C2')].r_median.iloc[0]) + ("*" if A_ == "A6" else "") for A_ in ARCH], ", ".join(br[A_] for A_ in ARCH) + "; *controllo interno smentito"))
T.append(row("sfocatura stretta per tutti σ_in (µm)", ["0.7"] * 6, "valore unico per scelta ([DEC-038](/decisions/DEC-038.md)): R5 PH-5 `results/R5/R5_PH5_narrowcore.csv` (sigma_fixed 0.7) e pre-registrazione di M2, dentro 0.5-1 µm di [DEC-037](/decisions/DEC-037.md); valori per tessuto in B-367...B-372"))
for p, nm, f_ in (("sigma_e", "zona di bordo: sfocatura larga σ_e (µm)", f2), ("D", "zona di bordo: profondità D (µm)", f2), ("w", "coda esponenziale: massa w", f3), ("ell", "coda esponenziale: lunghezza ℓ (µm)", f2)):
    v = [f_(e(SOF[A_], "C2", p)[0]) for A_ in ARCH[:5]] + [f_(HEART[p]) + "†"]
    T.append(row(nm, v, ", ".join(bc2[k] for k in SEC) + "; †mediana delle altre 4 sezioni ([DEC-037](/decisions/DEC-037.md)), calcolata da `results/S2.0/S2.0_edges_sections.csv`"))
T.append(row("fondo sul vetro b (UMI per bin di 2 µm)", [f3(BG[SOF[A_]][1]) for A_ in ARCH], ", ".join(BG[k][0] for k in SEC)))
T.append(row("strisce: parità righe dispari/pari", [f3(ST[SOF[A_]]["rows"]["odd_over_even"]) + ("‡" if A_ == "A4" else "") for A_ in ARCH], ", ".join(bst[k] for k in SEC) + "; ‡variabile nello spazio"))
T.append(row("strisce: parità colonne dispari/pari", [f3(ST[SOF[A_]]["cols"]["odd_over_even"]) + ("‡" if A_ == "A4" else "") for A_ in ARCH], ", ".join(bst[k] for k in SEC)))
T.append(row("strisce: SD dei fattori per riga α (log)", [f3(ST[SOF[A_]]["rows"]["sd_repro"]) for A_ in ARCH], ", ".join(bst[k] for k in SEC)))
T.append(row("strisce: SD dei fattori per colonna β (log)", [f3(ST[SOF[A_]]["cols"]["sd_repro"]) for A_ in ARCH], ", ".join(bst[k] for k in SEC)))
T.append(row("strisce: correlazione a passo 1 di α; β (misurata fra le due metà; lo stimatore aggiunge circa −0.03)", [f"{f2(ST[SOF[A_]]['rows']['lag1_repro'])}; {f2(ST[SOF[A_]]['cols']['lag1_repro'])}" for A_ in ARCH], ", ".join(bst[k] for k in SEC) + "; distorsione dallo stimatore: post hoc PH-4 `results/S2.0/S2.0_PH4_lag.csv`"))
T.append(row("variabilità dentro il nucleo v, come scritta da R-C (prototipo senza fattori di riga e colonna)", [f3(P.loc[A_].v_pre) for A_ in ARCH], ", ".join(bv[A_] for A_ in ARCH) + "; `results/S2.0/S2.0_RC.json`"))
T.append(row("variabilità dentro il nucleo v, al netto delle strisce (post hoc PH-2; proposta, scelta di Luca)", [f3(P.loc[A_].v) + ("*" if A_ == "A6" else "") for A_ in ARCH], ", ".join(bv[A_] for A_ in ARCH) + "; post hoc PH-2"))
T.append(row("scala della variabilità dentro il nucleo", ["1-2 µm"] * 6, ", ".join(bv[A_] for A_ in ARCH) + "; post hoc PH-3 (correlazione in più a 2 µm, al più 0.09 a 4 µm)"))
import re
BR = {l.split("|")[1].strip(): [x.strip() for x in l.split("|")[1:-1]] for l in open("docs/BIO_REFERENCES.md").read().split("\n") if l.startswith("| B-")}
REF = {"A1": (270, 271, 272), "A2": (280, 281, 282), "A3": (290, 291, 292), "A4": (302, 303, 304), "A5": (315, 316, 317), "A6": (327, 328, 329)}
def sp8(b):
    m = re.search(r"8 µm ([0-9.]+)$", BR[b][4].split("sparsità")[1].strip().rstrip(";")) if "sparsità" in BR[b][4] else None
    return m.group(1) if m else "—"
def sp2(b):
    m = re.search(r"sparsità 2 µm ([0-9.]+)", BR[b][4]); return m.group(1) if m else "—"
T.append(row("riferimento (controllo, non parametro): UMI per bin di 2 µm in tessuto, media", [BR[f"B-{REF[A_][0]}"][3] for A_ in ARCH], ", ".join(f"B-{REF[A_][0]}" for A_ in ARCH)))
T.append(row("riferimento (controllo, non parametro): UMI per bin di 8 µm in tessuto, media", [BR[f"B-{REF[A_][1]}"][3] for A_ in ARCH], ", ".join(f"B-{REF[A_][1]}" for A_ in ARCH)))
T.append(row("riferimento (controllo, non parametro): sparsità bin × gene a 2 µm; 8 µm", [f"{sp2('B-%d' % REF[A_][2])}; {sp8('B-%d' % REF[A_][2])}" for A_ in ARCH], ", ".join(f"B-{REF[A_][2]}" for A_ in ARCH)))
open(f"{RES}/S2.0_param_table.md", "w").write("\n".join(T) + "\n")
print("\n".join(T)); print(); print("\n".join(bio)); print("heart medians", HEART)
