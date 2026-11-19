"""tools/F3.3_docs.py - F3.3, documentazione senza previsione (solo dati di taratura e tabelle dei passi precedenti).
D-1: per tessuto, valori del simulatore stimati sulle 5 aree di taratura di R2 o sulle loro finestre (censimento di F3.2, results/F3.2/F3.2_census.csv)
contro le unita' indipendenti (5 aree, 1 sezione). D-2: giudizi confermati di F3.1 «costruiti dalla taratura» (un parametro stimato sulle stesse aree
per far tornare la stessa misura o una sua parte): M-2, UMI per cellula, costruita dal fattore di profondita' di S2.5 (regola R-2: UMI medio per
quadratino in tessuto delle stesse aree); M-1, UMI assegnati, «in parte» (densita' e taglia dei nuclei stimate sulle stesse aree, non tarate su M-1);
le altre no. Uso: python3 tools/F3.3_docs.py"""
import pandas as pd
C = pd.read_csv("results/F3.2/F3.2_census.csv"); C["n_values"] = pd.to_numeric(C.n_values, errors="coerce").fillna(1)
rows = []
for A in ["A1", "A2", "A3", "A4", "A5", "A6"]:
    c = C[C.archetype == A]; areas = c.area_set.astype(str).str.contains("R2_ROI|R2b_WINDOW|R2c")
    rows.append({"archetype": A, "rows_total": len(c), "rows_from_calibration_areas": int(areas.sum()), "values_from_calibration_areas": int(c.n_values[areas].sum()),
                 "rows_scalar_from_areas": int((areas & (c.kind == "scalar")).sum()), "rows_table_from_areas": int((areas & (c.kind == "table")).sum()),
                 "independent_areas": 5, "independent_sections": 1})
D1 = pd.DataFrame(rows); D1.to_csv("results/F3.3/F3.3_D1.csv", index=False); print(D1.to_string(index=False))
V = pd.read_csv("results/F3.1/F3.1_verdicts.csv"); V = V[V.verdict == "confermata"].copy()
V["built"] = V.id.map({"M-2": "costruita (fattore di profondita' di S2.5 tarato sull'UMI per quadratino delle stesse aree)",
                       "M-1": "in parte (densita' e taglia dei nuclei dalle stesse aree)"}).fillna("no")
V[["id", "archetype", "tissue", "verdict", "built"]].to_csv("results/F3.3/F3.3_D2.csv", index=False); print(V[["id", "archetype", "built"]].to_string(index=False))
