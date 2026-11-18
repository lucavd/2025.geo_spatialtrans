#!/usr/bin/env python3
"""tools/F3.2_census.py - F3.2, cancello G-1 (seconda parte): provenienza di ogni valore del censimento.
Legge results/F3.2/F3.2_census_values.csv (tools/F3.2_census.R) e docs/BIO_REFERENCES.md; scrive results/F3.2/F3.2_census.csv con, per riga:
fonte (riga B-nnn o file), sessione di stima, sezione, insieme di aree (R2_ROI, R2b_WINDOW, R5_EDGE_TILE, SECTION, DESIGN) e livelli di verifica
per cui il valore e' fuori campione (L1 aree nuove della stessa sezione; L2 altre sezioni dello stesso tessuto; L3 tessuti nuovi).
Le famiglie senza riga B-nnn nel codice prendono la fonte dalla documentazione roxygen delle funzioni (tabella SRC qui sotto, con i file citati).
Uso: python3 tools/F3.2_census.py
"""
import csv, re, sys
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results/F3.2"
SECTION_OF = {"A1": "A1", "A2": "A23", "A3": "A23", "A4": "A4", "A5": "A5", "A6": "A6", "all": "all"}
TECH = ("step2.transcripts", "step2.capture")                                # chimica e cattura: comuni ai tessuti per natura (L3 li verifica)
# famiglia, parametro (regex) -> fonte, sessione, insieme di aree, nota (dalla documentazione di transcript_params(), capture_params(),
# visium_hd_depth_factor(), load_cell_profiles() e dalle righe B-nnn che citano)
SRC = [
    ("step2.transcripts", r"^r$", "S2.0/wfix/S2.0_inner_summary.csv; B-384...B-394", "S2.0", "R2_ROI", "mediana delle 5 aree di R2"),
    ("step2.transcripts", r"^(sigma_e|D|w|ell)$", "S2.0/run1_wfix/params_flat.json, regola C2; B-379...B-383", "S2.0", "R5_EDGE_TILE",
     "tessere di bordo di R5 della sezione; cuore: mediana delle altre sezioni (DEC-037)"),
    ("step2.transcripts", r"^b$", "R5, fondo del vetro lontano; B-342...B-366", "R5", "SECTION", "conservato per confronto, non usato dal simulatore"),
    ("step2.transcripts", r"^sigma_in$", "DEC-038 (sfocatura stretta unica)", "S2.0", "DESIGN", "valore unico per tutti i tessuti"),
    ("step2.transcripts", r"^field_(s2|L)$", "results/S2.3/S2.3_R1_choice.csv (regola R-1, taratura sul nucleo)", "S2.3", "R2_ROI", "30 aree di R2"),
    ("step2.capture", r"^(parity_row|parity_col|sd_alpha|sd_beta)$", "results/S2.0/S2.0_PH5_sections.csv; B-397...B-401", "S2.0", "SECTION", "sezione intera"),
    ("step2.capture", r"^(rho_alpha|rho_beta)$", "results/S2.4/S2.4_R1_generator.csv (regola R-1)", "S2.4", "SECTION", "dalle strisce della sezione intera"),
    ("step2.capture", r"^(b|w2|ell2)$", "results/S2.4/S2.4_R3_fits.csv (regola R-3, DEC-053)", "S2.4", "SECTION", "vetro vicino al tessuto, sezione intera"),
    ("step2.capture", r"^(sd_row|sd_col|corr_um)$", "documento di progetto dello Step 2, sezione 7 (campo della parita' del linfonodo)", "S2.0", "SECTION", "tessere di 300 um della sezione"),
    ("step2.depth", r"^depth_factor$", "results/S2.5/S2.5_R2_depth.csv (regola R-2); B-518...B-523", "S2.5", "R2_ROI", "30 aree di R2"),
    ("step2.profiles", r".*", "tools/S2.2_build_profiles.py (nuclei veri delle aree di R2); results/S2.2/S2.2_profiles_md5.txt", "S2.2", "R2_ROI", "nuclei delle 5 aree di R2 del tessuto"),
    ("cell_layer.size_resid(option)", r".*", "tools/S1.5d_size_resid_table.py (residui veri del modello di taglia)", "S1.5d", "R2_ROI", "celle interne delle 5 aree di R2; opzione, non predefinita"),
    ("cell_layer.defaults", r".*", "scelte del design (DEC-019, DEC-022, DEC-024)", "S1.5b-S1.5e", "DESIGN", "impostazione, non stima"),
    ("step2.grid", r".*", "documento di progetto dello Step 2 (margini di 280 um, DEC-038; griglia di 2 um)", "S2.1", "DESIGN", "impostazione, non stima"),
    ("step2.pieces", r".*", "DEC-074 (tessere fisse di 56 um)", "S2.7", "DESIGN", "impostazione, non stima"),
]

def brows():
    rows = {}
    for ln in open(ROOT / "docs/BIO_REFERENCES.md", encoding="utf-8"):
        if ln.startswith("| B-"):
            c = [x.strip() for x in ln.strip().strip("|").split("|")]
            rows[c[0]] = dict(id=c[0], metric=c[2], value=c[3], source=c[6], session=c[7], note=c[8] if len(c) > 8 else "")
    return rows

def area_set_of_b(b):
    s = b["source"] + " " + b["note"] + " " + b["session"]
    if b["session"] == "R2b" or "finestre" in s and "R2b" in s: return "R2b_WINDOW"
    if re.search(r"ROI R2|5 ROI|30 ROI|maschere native R2|R2_rois_checked|R6|pooled|celle interne", s): return "R2_ROI"
    return "UNKNOWN"

def main():
    v = pd.read_csv(OUT / "F3.2_census_values.csv", dtype=str, keep_default_na=False)
    B = brows(); out = []
    for _, r in v.iterrows():
        rec = r.to_dict(); src = r["src_in_code"]
        if src.startswith("B-"):
            b = B.get(src)
            if b is None: rec.update(source_id=src, source_session="", source_text="RIGA ASSENTE", area_set="UNKNOWN", note="")
            else: rec.update(source_id=src, source_session=b["session"], source_text=b["source"][:160], area_set=area_set_of_b(b), note=b["note"][:120])
        elif src.endswith(".csv"):        # density_model: S1.4_calibration.csv / R3b_refs.csv
            rec.update(source_id=src, source_session="S1.4" if "S1.4" in src else "R3b", source_text=f"results/{'S1.4' if 'S1.4' in src else 'R3b'}/{src}",
                       area_set="R2_ROI", note="densita' con cui e' stato stimato il modello di taglia")
        else:
            hit = [s for s in SRC if s[0] == r["family"] and re.search(s[1], r["parameter"].split(":")[-1])]
            if hit:
                s = hit[0]; rec.update(source_id=s[2].split(";")[-1].strip(), source_session=s[3], source_text=s[2], area_set=s[4], note=s[5])
            else:
                rec.update(source_id="", source_session="", source_text="", area_set="UNKNOWN", note="")
        rec["section"] = SECTION_OF[r["archetype"]]
        a = rec["area_set"]
        rec["L1_out_of_sample"] = a in ("R2_ROI", "R2b_WINDOW", "R5_EDGE_TILE")      # le aree nuove stanno a >= 280 um da queste
        rec["L2_out_of_sample"] = a != "DESIGN"
        rec["L3_out_of_sample"] = r["family"] in TECH and a != "DESIGN"
        out.append(rec)
    d = pd.DataFrame(out); d.to_csv(OUT / "F3.2_census.csv", index=False)
    n_unk = int((d.area_set == "UNKNOWN").sum())
    print("righe", len(d), "senza insieme di aree", n_unk)
    print(d.groupby(["area_set"]).size().to_string())
    print(d.groupby(["family", "area_set"]).size().to_string())
    if n_unk: print(d[d.area_set == "UNKNOWN"][["family", "archetype", "parameter", "src_in_code"]].to_string())

if __name__ == "__main__":
    main()
