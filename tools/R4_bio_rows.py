"""tools/R4_bio_rows.py — righe di docs/BIO_REFERENCES.md da results/R4 (stampate su stdout, una per riga, da B-<primo>).
Uso: .venv/bin/python tools/R4_bio_rows.py <primo_id_numerico>"""
import sys
import pandas as pd
first = int(sys.argv[1]); RES = "results/R4"
R = pd.read_csv(f"{RES}/R4_refs.csv"); T = pd.read_csv(f"{RES}/R4_types.csv"); PH1 = pd.read_csv(f"{RES}/R4_PH1.csv")
P = pd.read_csv(f"{RES}/R4_predictions.csv"); W = pd.read_csv(f"{RES}/R4_PH4_parity.csv"); D4 = pd.read_csv(f"{RES}/R4_D4.csv")
DS = {"A1": "Visium_HD_Mouse_Small_Intestine", "A2": "Visium_HD_Human_Colon_Cancer", "A3": "Visium_HD_Human_Colon_Cancer", "A4": "Visium_HD_Human_Lymph_Node_FFPE",
      "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}
SRC = "tools/R4_analysis.py@R4"
def g(a, q):
    r = R[(R.archetype == a) & (R.quantity == q)].iloc[0]; return r
def ci(r, d=3): return f"{r.estimate:.{d}f} [{r.lo:.{d}f}, {r.hi:.{d}f}]"
rows = []; k = first
def add(a, metric, value, interval, unit, src, note):
    global k
    rows.append(f"| B-{k} | {a} | {metric} | {value} | {interval} | {unit} | {src} | R4 | {note} |"); k += 1
for a in ["A1", "A2", "A3", "A4", "A5", "A6"]:
    base = f"{SRC} (results/R4/R4_refs.csv) + {DS[a]}, 5 ROI di R2, bin in tessuto di Space Ranger"
    for res in ("2", "8"):
        add(a, f"UMI per bin a {res} µm in tessuto (media; mediana, 90° e 99° percentile, CV, frazione di bin vuoti)", f"{g(a, f'mean_umi{res}').estimate:.2f}",
            f"media {ci(g(a, f'mean_umi{res}'), 2)}; mediana {ci(g(a, f'median_umi{res}'), 1)}; q90 {g(a, f'q90_umi{res}').estimate:.1f}; q99 {g(a, f'q99_umi{res}').estimate:.1f}; "
            f"CV {g(a, f'cv_umi{res}').estimate:.3f}; vuoti {ci(g(a, f'frac_empty{res}'), 4)}", "UMI per bin; frazione", base,
            "media delle 5 ROI, bootstrap sulle ROI (B 2000, seme 20261030); a 2 µm include l'alternanza di righe e colonne (B per parità sotto)")
    add(a, "geni rilevati per bin a 2 µm e a 8 µm; sparsità (frazione di valori nulli bin × gene)", f"{g(a, 'mean_genes2').estimate:.2f} / {g(a, 'mean_genes8').estimate:.1f}",
        f"2 µm {ci(g(a, 'mean_genes2'), 2)}; 8 µm {ci(g(a, 'mean_genes8'), 1)}; sparsità 2 µm {g(a, 'sparsity2').estimate:.4f}, 8 µm {g(a, 'sparsity8').estimate:.4f}", "geni per bin; frazione", base,
        "tutte le feature Gene Expression della matrice filtrata")
    add(a, "UMI per bin a 2 µm per strato: nucleo, corona (≤ 8 µm dal bordo del nucleo), acellulare (> 8 µm); quota di bin per strato",
        f"{g(a, 'mean_umi2_nucleus').estimate:.2f} / {g(a, 'mean_umi2_ring').estimate:.2f} / {g(a, 'mean_umi2_acellular').estimate:.2f}",
        f"quote {g(a, 'share_nucleus').estimate:.3f} / {g(a, 'share_ring').estimate:.3f} / {g(a, 'share_acellular').estimate:.3f}; bin vuoti {g(a, 'frac_empty2_nucleus').estimate:.4f} / {g(a, 'frac_empty2_ring').estimate:.4f} / {g(a, 'frac_empty2_acellular').estimate:.4f}",
        "UMI per bin; frazione", base + "; nuclei di R6 (S1.5 refs), impronta 8 µm", "media delle 5 ROI")
    cw = PH1[(PH1.archetype == a) & (PH1.quantity == "cellweighted_star")].iloc[0]; p1 = P[(P.id == "P-1") & (P.archetype == a)].iloc[0]
    add(a, "rapporto fra UMI per bin nel nucleo e nella sua corona di 8 µm, pesato per cellula, diviso per lo stesso rapporto con la maschera spostata di 20 µm",
        f"{cw.estimate:.3f}", f"[{cw.lo:.3f}, {cw.hi:.3f}]; pre-registrato (bin raggruppati) {p1.estimate:.3f} [{p1.lo:.3f}, {p1.hi:.3f}], regola R-1: {p1.rule_R1}", "adimensionale",
        "tools/R4_posthoc.py@R4 (results/R4/R4_PH1.csv) + " + DS[a], "riferimento per lo Step 2 (DEC-034); osservato dopo diffusione e registrazione: senza i bin entro 2 µm dal bordo il rapporto e' piu' lontano da 1 (results/R4/R4_sensitivity.csv)")
    add(a, "UMI per cellula (nucleo + corona di 8 µm), cellule interne: mediana; 1° e 99° percentile; CV; massimo",
        f"{g(a, 'cell_umi_median').estimate:.1f}", f"mediana {ci(g(a, 'cell_umi_median'), 1)}; q01 {g(a, 'cell_umi_q01').estimate:.1f}; q99 {g(a, 'cell_umi_q99').estimate:.1f}; "
        f"CV {ci(g(a, 'cell_umi_cv'), 3)}; massimo {g(a, 'cell_umi_max').estimate:.0f}; UMI per bin della cellula {ci(g(a, 'cell_umi_per_bin_median'), 2)}", "UMI per cellula",
        base + "; nuclei di R6", f"bootstrap a due stadi (ROI, cellule); n = {int(g(a, 'cell_umi_median').n_cells)} cellule")
    for _, t in T[T.archetype == a].sort_values("cell_type").iterrows():
        add(a, f"UMI per cellula per tipo: **{t.cell_type}** (mediana; UMI per bin della cellula)", f"{t.median_umi_cell:.1f}",
            f"[{t.lo:.1f}, {t.hi:.1f}]; per bin {t.median_umi_per_bin:.2f} [{t.lo_per_bin:.2f}, {t.hi_per_bin:.2f}]; quota dei nuclei interni {t.share:.3f}", "UMI per cellula; UMI per bin",
            base + "; etichette di R6 (PH-1, DEC-013)", f"n = {t.n_cells} cellule interne; {'nel confronto fra tipi (quota >= 5 %)' if t.in_P4 else 'quota < 5 %'}")
    sec = W[(W.unit == "section") & (W.dataset == DS[a])].iloc[0]; d4 = D4[(D4.archetype == a) & (D4.unit == "square_002um")].iloc[0]
    add(a, "alternanza della cattura a 2 µm: UMI per bin delle colonne dispari / pari e delle righe dispari / pari (sezione intera, bin in tessuto)",
        f"{sec.col_odd_over_even:.3f} / {sec.row_odd_over_even:.3f}", f"ROI: colonne {W[(W.unit.str.startswith(a + '_'))].col_odd_over_even.min():.3f}-{W[(W.unit.str.startswith(a + '_'))].col_odd_over_even.max():.3f}, righe "
        f"{W[(W.unit.str.startswith(a + '_'))].row_odd_over_even.min():.3f}-{W[(W.unit.str.startswith(a + '_'))].row_odd_over_even.max():.3f}; media della sezione {d4.section_mean_umi:.2f} UMI per bin a 2 µm, ROI / sezione {d4.ratio_roi_to_section:.2f}",
        "adimensionale", "tools/R4_posthoc.py@R4 (results/R4/R4_PH4_parity.csv) + " + DS[a], "effetto tecnico della griglia (strisce per riga e colonna, Polanski et al. 2024 doi:10.1093/bioinformatics/btae546); si annulla a 8 µm")
print("\n".join(rows))
print(f"# last B-{k-1}", file=sys.stderr)
