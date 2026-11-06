"""S2.2: righe nuove di docs/BIO_REFERENCES.md dalle tabelle di S2.2 (numeri letti dai file). Inserite dopo l'ultima riga B-nnn.
Uso: .venv/bin/python tools/S2.2_bio_rows.py [--write]"""
import sys, re, pandas as pd
R = "results/S2.2"; sha = sys.argv[sys.argv.index("--sha") + 1] if "--sha" in sys.argv else "HEAD"
e = pd.read_csv(f"{R}/S2.2_E1_real.csv"); f = pd.read_csv(f"{R}/S2.2_depth_factor.csv").set_index("archetype")
DS = {"A1": "Visium_HD_Mouse_Small_Intestine", "A2": "Visium_HD_Human_Colon_Cancer", "A3": "Visium_HD_Human_Colon_Cancer", "A4": "Visium_HD_Human_Lymph_Node_FFPE",
      "A5": "Visium_HD_6p5mm_Mouse_Brain", "A6": "Visium_HD_6p5mm_Human_Heart"}
txt = open("docs/BIO_REFERENCES.md").read()
last = max(int(x) for x in re.findall(r"^\| B-(\d+) ", txt, flags=re.M)); n = last; rows = []
for A, g in e.groupby("archetype"):
    for col, name, unit in (("m3", "somiglianza fra blocchi di 8 µm adiacenti nei quadratini dei nuclei con tipo (coseno medio)", "—"),
                            ("m4_ratio", "misura di Daniele: coseno fra le due metà dei quadratini della stessa cellula / fra cellule diverse dello stesso tipo (> 30 µm)", "—")):
        n += 1
        rows.append(f"| B-{n:03d} | {A} | {name} | {g[col].mean():.3f} | min–max 5 ROI: {g[col].min():.3f}–{g[col].max():.3f}; n quadratini {int(g.n_bins.sum())} | {unit} | tools/S2.2_real.py@{sha} + {DS[A]} | S2.2 | impronta nucleo + 8 µm (R6), tipi del preset; ROI sovrapposte (BL-215): intervallo ottimistico |")
for A in f.index:
    n += 1; r = f.loc[A]
    rows.append(f"| B-{n:03d} | {A} | fattore di profondità del simulatore (UMI veri per bin di 2 µm in tessuto / UMI attesi dalle cellule simulate dei preset con i totali per tipo di R6) | {r.f:.3f} | IC 95 % bootstrap 5 ROI: {r.lo:.3f}–{r.hi:.3f}; min–max {r.f_min:.3f}–{r.f_max:.3f} | — | tools/S2.2_depth.py@{sha} + {DS[A]} | S2.2 | regola R-3; scomposizione in results/S2.2/S2.2_PH2.csv; da ricalcolare con le densità nuove (BL-218); ROI sovrapposte (BL-215) |")
print("\n".join(rows))
if "--write" in sys.argv:
    lines = txt.split("\n"); idx = max(i for i, l in enumerate(lines) if l.startswith(f"| B-{last:03d} "))
    lines[idx + 1:idx + 1] = rows; open("docs/BIO_REFERENCES.md", "w").write("\n".join(lines))
    print(f"scritte {len(rows)} righe: B-{last+1:03d}...B-{n:03d}")
