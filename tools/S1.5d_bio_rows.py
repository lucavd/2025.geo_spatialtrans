"""tools/S1.5d_bio_rows.py - S1.5d: righe nuove di docs/BIO_REFERENCES.md (B-232...) e valori dei preset in R/04b0_tissue_presets.R, generati dai file dei risultati.
- compartimenti ricalibrati (piega «pooled», results/S1.5d/S1.5d_compartments_<A>.csv), tipi principali: .TP_COMPART (ell_um, s, src);
- tabella dei residui del modello di taglia (inst/extdata/size_resid.csv, md5) per gruppo: una riga per gruppo;
- eterogeneita' della famiglia vincente della regola A (results/S1.5d/S1.5d_ruleA.csv; S1.5d_het_<famiglia>.csv, piega «pooled»): .TP_HET.
Uso: .venv/bin/python tools/S1.5d_bio_rows.py <sha della calibrazione>"""
import sys, re, hashlib, pandas as pd, numpy as np
SHA = sys.argv[1]; P = {"A1": "epithelial", "A2": "tumor_microenv", "A3": "stromal_rich", "A4": "lymphoid", "A5": "brain_cortex", "A6": "columnar_muscle"}
bio = open("docs/BIO_REFERENCES.md").read(); last = max(int(x) for x in re.findall(r"^\| B-(\d{3}) \|", bio, re.M)); assert last == 231, last
rows, nid = [], last + 1; comp_src, het_src = {}, {}
for A in P:
    d = pd.read_csv(f"results/S1.5d/S1.5d_compartments_{A}.csv")
    for t, g in d.groupby("cell_type", sort=False):
        p = g[g.fold == "pooled"].iloc[0]; l = g[g.fold != "pooled"]; bid = f"B-{nid}"; nid += 1; comp_src[(P[A], t)] = (bid, p.ell_um, p.s)
        edge = "al bordo della griglia" if (p.edge_ell or p.edge_s) else "dentro la griglia"
        rows.append(f"| {bid} | {A} | compartimenti per tipo ricalibrati (S1.5d, BL-159): **{t}** (ampiezza s del campo latente; lunghezza di correlazione ell) | {p.s:g} | ell={p.ell_um:g}; pieghe lascia-fuori: s {l.s.min():g}-{l.s.max():g}, ell {l.ell_um.min():g}-{l.ell_um.max():g} µm; perdita {p.loss:.1f} contro {p.loss_null:.1f} senza campi | s adimensionale; ell µm | tools/S1.5d_calibration_types.py@{SHA} (piega «pooled», 5 ROI, 20 repliche) + etichette R6 PH-1 (DEC-013), 5 dataset Visium HD | S1.5d | {edge} (griglia s 0-8, ell 4-1024 µm); sostituisce la riga di S1.5c dello stesso tipo nel preset |")
md5 = hashlib.md5(open("inst/extdata/size_resid.csv", "rb").read()).hexdigest(); ch = pd.read_csv("results/S1.5d/S1.5d_size_resid_check.csv")
for _, r in ch.iterrows():
    bid = f"B-{nid}"; nid += 1
    rows.append(f"| {bid} | {r.archetype} | residui veri del modello di taglia, gruppo {r.size_group} ({r.types.replace('|', ', ')}): 10 decimi di log A_T x 101 quantili, tabella inst/extdata/size_resid.csv (md5 {md5}) | {r.n} celle | media dei residui {r.resid_mean:.2g}; deviazione standard {r.resid_sd:.4f} (sigma del gruppo {r.sigma_g:.4f}) | log µm² | tools/S1.5d_size_resid_table.py@{SHA} + {'nuclei di forma StarDist della maschera CO di R3b (DEC-006)' if r.archetype in ('A4', 'A6') else 'nuclei veri con territorio di R6 (insieme di S1.5b)'}, 5 ROI | S1.5d | candidato M3 (regola C); include l'errore di segmentazione, come sigma |")
win = pd.read_csv("results/S1.5d/S1.5d_ruleA.csv").winner.iloc[0]; H = pd.read_csv(f"results/S1.5d/S1.5d_het_{win}.csv")
for A in P:
    h = H[(H.archetype == A) & (H.fold == "pooled")].iloc[0]; l = H[(H.archetype == A) & (H.fold != "pooled")]; bid = f"B-{nid}"; nid += 1
    fmt = lambda v: "NA" if not np.isfinite(v) else f"{v:g}"
    val = f"{h.void_frac:g} / {fmt(h.void_ell_um)} / {h.dens_sigma:g} / {fmt(h.dens_ell_um)}"; het_src[P[A]] = (bid, h)
    rows.append(f"| {bid} | {A} | eterogeneita' della semina (S1.5d, famiglia {win} della regola A): void_frac / void_ell_um / dens_sigma / dens_ell_um | {val} | pieghe lascia-fuori: void_frac {l.void_frac.min():g}-{l.void_frac.max():g}; perdita {h.loss:.1f} contro {h.loss_null:.1f} senza eterogeneita'; bordo: {h.edge if isinstance(h.edge, str) else 'no'} | frazione; µm; adimensionale; µm | tools/S1.5d_calibration_pos.R + tools/S1.5d_calibration_pos.py@{SHA} (piega «pooled», 30 ROI x 2 repliche, griglia estesa una volta), maschere di R2 | S1.5d | bersagli: spazio vuoto H_rel 1.5-4, CV2 in eccesso 25-100 um, CV della distanza dal vicino |")
k = bio.index("| B-231 |"); k = bio.index("\n", k) + 1
open("docs/BIO_REFERENCES.md", "w").write(bio[:k] + "\n".join(rows) + "\n" + bio[k:])
tp = open("R/04b0_tissue_presets.R").read()
m = re.search(r"\.TP_COMPART <- data\.frame\((.*?)\n\)\n", tp, re.S); block = m.group(0)
def vec(name):
    mm = re.search(rf"  {name} = c\((.*?)\),?\n", block); return mm
pres = re.findall(r'"([^"]+)"', vec("preset").group(1)); ctyp = re.findall(r'"([^"]+)"', vec("cell_type").group(1))
ell = [x.strip() for x in vec("ell_um").group(1).split(",")]; s = [x.strip() for x in vec("s").group(1).split(",")]; src = re.findall(r'"([^"]+)"', vec("src").group(1))
for i, (pp, tt) in enumerate(zip(pres, ctyp)):
    if (pp, tt) in comp_src: b_, e_, s_ = comp_src[(pp, tt)]; ell[i] = f"{e_:g}"; s[i] = f"{s_:g}"; src[i] = b_
nb = block
for name, vals, q in [("ell_um", ell, False), ("s", s, False), ("src", src, True)]:
    mm = vec(name); new = ", ".join(f'"{v}"' if q else v for v in vals); nb = nb.replace(mm.group(0), mm.group(0).replace(mm.group(1), new))
tp = tp.replace(block, nb)
old = re.search(r"\.TP_HET <- data\.frame\(.*?stringsAsFactors = FALSE\)\n", tp, re.S).group(0)
f = lambda v: "NA" if not np.isfinite(v) else repr(float(v))
newh = (".TP_HET <- data.frame(preset = c(" + ", ".join(f'"{p}"' for p in P.values()) + "),\n"
        "                      void_frac = c(" + ", ".join(f(het_src[p][1].void_frac) for p in P.values()) + "), void_ell_um = c(" + ", ".join(f(het_src[p][1].void_ell_um) for p in P.values()) + "),\n"
        "                      dens_sigma = c(" + ", ".join(f(het_src[p][1].dens_sigma) for p in P.values()) + "), dens_ell_um = c(" + ", ".join(f(het_src[p][1].dens_ell_um) for p in P.values()) + "),\n"
        "                      src = c(" + ", ".join(f'"{het_src[p][0]}"' for p in P.values()) + "), stringsAsFactors = FALSE)\n")
tp = tp.replace(old, newh)
tp = tp.replace("# Eterogeneita': calibrazione su tutte le ROI della famiglia vincente della regola A (docs/steps/S1.5d.md); prima della calibrazione tutto 0\n# (nessun effetto: la semina a compartimenti resta quella di S1.5c bit per bit).",
                f"# Eterogeneita': calibrazione su tutte le ROI della famiglia vincente della regola A ({win}, docs/steps/S1.5d.md; tools/S1.5d_bio_rows.py).")
tp = tp.replace("# ---- S1.5c: compartimenti per tipo (docs/steps/S1.5c.md; riferimenti B-209...B-231) ----", "# ---- S1.5c: compartimenti per tipo (docs/steps/S1.5c.md; riferimenti B-209...B-231; tipi principali ricalibrati in S1.5d, B-232...B-248) ----")
open("R/04b0_tissue_presets.R", "w").write(tp); print("righe", nid - last - 1, "ultima", f"B-{nid-1}", "vincente", win)
