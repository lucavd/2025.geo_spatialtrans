#!/usr/bin/env python3
"""tools/R6_report_parts.py - estrae dal file dello step (docs/steps/R6.md) le sezioni incluse nel report e costruisce la tabella
«Da dove vengono i numeri del riassunto» rileggendo ogni valore dai file dei risultati. Uso (root del repo): .venv/bin/python -B tools/R6_report_parts.py"""
import pathlib, re, pandas as pd
R = pathlib.Path("results/R6"); out = pathlib.Path("reports")
t = pathlib.Path("docs/steps/R6.md").read_text(encoding="utf-8").split("\n---\n", 1)[1]
def sec(start, stop_re):
    i = t.index(start); m = re.search(stop_re, t[i + len(start):]); j = i + len(start) + (m.start() if m else len(t))
    return t[i:j].strip() + "\n"
(out / "_R6_inbreve.md").write_text(sec("\n# In breve\n", r"\n# "), encoding="utf-8")
(out / "_R6_rinviati.md").write_text(sec("\n# Dubbi rinviati\n", r"\n# ").replace("# Dubbi rinviati\n", "", 1), encoding="utf-8")
(out / "_R6_verifiche.md").write_text(sec("\n# Verifiche di Luca\n", r"\n# ").replace("# Verifiche di Luca\n", "", 1), encoding="utf-8")
# ---- fonti dei numeri del riassunto
P0 = pd.read_csv(R / "R6_predictions.csv"); P1 = pd.read_csv(R / "R6_predictions_ph1.csv")
T1 = pd.read_csv(R / "R6_refs_type_ph1.csv").set_index(["archetype", "cell_type"]); T0 = pd.read_csv(R / "R6_refs_type.csv").set_index(["archetype", "cell_type"])
PH = pd.read_csv(R / "R6_posthoc.csv"); K = pd.read_csv(R / "R6_posthoc_kappa.csv"); RT = pd.read_csv(R / "R6_refs_roi_type_ph1.csv")
cl1 = pd.read_csv("/mnt/micron/geo_spatialtrans/R6/labels/A1_exp_umi20_seed20261022_clusters.csv")
cl6 = pd.read_csv("/mnt/micron/geo_spatialtrans/R6/labels/A6_exp_umi20_seed20261022_clusters.csv")
def p(df, i, a): return df[(df.id == i) & (df.archetype == a)].iloc[0]
def ph(q, a=None):
    d = PH[PH.quantity.str.startswith(q)]; d = d[d.archetype == a] if a else d; return d.iloc[0]
f = lambda x, d=2: f"{x:.{d}f}"
rows = [
 ("240 000 nuclei veri delle 30 aree", "R6_refs_type_ph1.csv, somma di n", f"{int(T1.n.sum()):,}".replace(",", " ")),
 ("cellule di Paneth chiamate «endoteliali» (5 260 nuclei)", "labels/A1_exp_..._clusters.csv, cluster 13 (Defa21, Lyz1, Itln1)", f"{int(cl1.loc[cl1.cluster == 13, 'n_nuclei'].iloc[0])} nuclei, type = {cl1.loc[cl1.cluster == 13, 'type'].iloc[0]}"),
 ("tre gruppi di cardiomiociti senza nome", "labels/A6_exp_..._clusters.csv, type = non risolto", f"{int((cl6.type == 'non risolto').sum())} cluster, {int(cl6.loc[cl6.type == 'non risolto', 'n_nuclei'].sum())} nuclei"),
 ("accordo con Space Ranger 91 % e 79 % (corretto), 88 % e 77 % (pre-registrato)", "R6_predictions_ph1.csv e R6_predictions.csv, P-6 A5/A6", f"{p(P1,'P-6','A5').estimate:.3f}, {p(P1,'P-6','A6').estimate:.3f}; {p(P0,'P-6','A5').estimate:.3f}, {p(P0,'P-6','A6').estimate:.3f}"),
 ("impronte spostate di 20 µm: accordo al livello del caso", "R6_posthoc_kappa.csv (PH-2), kappa primaria contro spostata", "; ".join(f"{r.archetype} {r.footprint}: {r.kappa:.2f}" for _, r in K[K.variant == 'ph1'].iterrows())),
 ("neuroni circa 3 volte la glia (3,2 corretto; 2,5 pre-registrato, al limite)", "R6_predictions.csv / _ph1.csv, P-1", f"{p(P0,'P-1','A5').estimate:.2f} [{p(P0,'P-1','A5').lo:.2f}, {p(P0,'P-1','A5').hi:.2f}]; {p(P1,'P-1','A5').estimate:.2f} [{p(P1,'P-1','A5').lo:.2f}, {p(P1,'P-1','A5').hi:.2f}]"),
 ("neuroni 111 µm², glia 34 (astrociti 50, oligodendrociti 31, microglia 23)", "R6_refs_type_ph1.csv A5; R6_posthoc.csv glia pooled", f"{T1.loc[('A5','neuroni'),'median_area_nuc']:.0f}; {ph('median nuclear area glia').estimate:.0f}; {T1.loc[('A5','astrociti'),'median_area_nuc']:.0f}, {T1.loc[('A5','oligodendrociti'),'median_area_nuc']:.0f}, {T1.loc[('A5','microglia'),'median_area_nuc']:.0f}"),
 ("neuroni due nuclei su tre", "R6_refs_type_ph1.csv A5 neuroni fraction", f"{T1.loc[('A5','neuroni'),'fraction']:.2f}"),
 ("oligodendrociti dall'8 % al 36 % fra le aree", "R6_refs_roi_type_ph1.csv A5 oligodendrociti", f"{RT[(RT.archetype=='A5')&(RT.cell_type=='oligodendrociti')].fraction.min():.2f} - {RT[(RT.archetype=='A5')&(RT.cell_type=='oligodendrociti')].fraction.max():.2f}"),
 ("differenze fra tipi del 24-42 % in intestino, tumore, stroma", "R6_predictions_ph1.csv P-2 A1, A2, A3", ", ".join(f"{a} {p(P1,'P-2',a).estimate:.2f}" for a in ['A1','A2','A3'])),
 ("linfonodo: B e T differiscono del 5 % (16,8 contro 16,0 µm²)", "R6_refs_type_ph1.csv A4; P-2 A4", f"{T1.loc[('A4','linfociti B'),'median_area_nuc']:.1f}, {T1.loc[('A4','linfociti T'),'median_area_nuc']:.1f}; P-2 {p(P1,'P-2','A4').estimate:.3f}"),
 ("cuore: differenza fra tipi 11 %; cardiomiociti da 16 a 28 µm² fra le aree", "P-2 A6 (_ph1); R6_refs_roi_type_ph1.csv A6 cardiomiociti", f"{p(P1,'P-2','A6').estimate:.2f}; {RT[(RT.archetype=='A6')&(RT.cell_type=='cardiomiociti')].median_area_nuc.min():.1f} - {RT[(RT.archetype=='A6')&(RT.cell_type=='cardiomiociti')].median_area_nuc.max():.1f}"),
 ("tumore del colon: dispersione dentro il solo epitelio 0,56 (tutte le cellule 0,55)", "R6_posthoc.csv PH-6 A2 epiteliale normale; BL-138 (sigma 0.55)", f"{ph('epiteliale normale', 'A2').estimate:.2f}"),
 ("cardiomiociti non più allungati: 0,82 contro 0,83", "R6_refs_type_ph1.csv A6 median_ecc_nuc; P-3", f"{T1.loc[('A6','cardiomiociti'),'median_ecc_nuc']:.3f}, fibroblasti {T1.loc[('A6','fibroblasti'),'median_ecc_nuc']:.3f}, immunitari {T1.loc[('A6','immunitari'),'median_ecc_nuc']:.3f}; P-3 {p(P1,'P-3','A6').estimate:.3f} [{p(P1,'P-3','A6').lo:.3f}, {p(P1,'P-3','A6').hi:.3f}]"),
 ("intestino: epitelio 68 %, stroma e muscolo 21 %, immunitarie 10 %", "R6_refs_type_ph1.csv A1", f"{T1.loc[('A1','epiteliale'),'fraction']:.2f}, {T1.loc[('A1','stromale e muscolo liscio'),'fraction']:.2f}, {T1.loc[('A1','immunitario'),'fraction']:.2f}"),
 ("tumore del colon: epitelio 71 %, fibroblasti 20 %", "R6_posthoc.csv PH-4 A2; R6_refs_type_ph1.csv A2", f"{ph('fraction epithelial', 'A2').estimate:.2f} [{ph('fraction epithelial', 'A2').lo:.2f}, {ph('fraction epithelial', 'A2').hi:.2f}]; {T1.loc[('A2','fibroblasti'),'fraction']:.2f}"),
 ("stroma: fibroblasti 32 %, epitelio 41 %, immunitarie 19 %", "R6_refs_type_ph1.csv A3; PH-4 A3", f"{T1.loc[('A3','fibroblasti'),'fraction']:.2f}, {ph('fraction epithelial', 'A3').estimate:.2f}, {T1.loc[('A3','immunitario'),'fraction']:.2f}"),
 ("linfonodo: B 61 % nei follicoli, 42 % nella paracorticale; 46 % senza tipo; 84 trascritti per nucleo", "R6_posthoc.csv PH-5; R6_refs_type_ph1.csv A4 non risolto; median_umi", f"{ph('fraction linfociti B among resolved, follicular').estimate:.2f}, {ph('fraction linfociti B among resolved, paracortical').estimate:.2f}; {T1.loc[('A4','non risolto'),'fraction']:.2f}; umi mediano A4 {RT[RT.archetype=='A4'].median_umi.median():.0f}"),
 ("cuore: cardiomiociti 76 %, fibroblasti 13 %, immunitarie e vascolari 11 %", "R6_refs_type_ph1.csv A6", f"{T1.loc[('A6','cardiomiociti'),'fraction']:.2f}, {T1.loc[('A6','fibroblasti'),'fraction']:.2f}, {T1.loc[('A6','immunitari'),'fraction']:.2f}"),
]
md = "| frase del riassunto | fonte (tabella o file, riga) | valore letto |\n|---|---|---|\n" + "\n".join(f"| {a} | {b} | {c} |" for a, b, c in rows) + "\n"
(out / "_R6_fonti.md").write_text(md, encoding="utf-8"); print(md)
