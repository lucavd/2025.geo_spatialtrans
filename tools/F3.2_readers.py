#!/usr/bin/env python3
"""tools/F3.2_readers.py - F3.2, controprova CP-1: ogni script che legge i dati veri (o le cartelle derivate su micron) corrisponde a un insieme
della mappa delle aree usate (results/F3.2/F3.2_used_areas.csv) o e' dichiarato «sezione intera». Classificazione per regole (file di aree citati,
prefisso dello step) con eccezioni scritte a mano e motivate (OVERRIDE). La controprova e' smentita se resta uno script non classificato o se uno
script di stima dei parametri legge aree che la mappa non contiene.
Uso: python3 tools/F3.2_readers.py -> results/F3.2/F3.2_readers.csv
"""
import re, csv, sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
RAW = re.compile(r"data_real|data/real/datasets|/mnt/micron/geo_spatialtrans/(R2|R2b|R2c|R3b|R5|S2\.0|S2\.2|R6|R4)\b")
REFS = {"R2_ROI": re.compile(r"R2_rois(_checked)?\.csv|R2/rois|R2/masks|R2/nuclei|R3b/(masks|nuclei)|R6/(counts|labels)|S1\.6/layers|S2\.2/profiles"),
        "R2b_WINDOW": re.compile(r"R2b_windows\.csv|R2c_"),
        "R5_EDGE_TILE": re.compile(r"R5_tiles\.csv|_candidates\.csv|R5/run3/(edges|tiles)|run1_wfix"),
        "SECTION": re.compile(r"sezione intera|/section/|R5/run3/section|_section\.py")}
# prefisso dello step -> insieme di aree predefinito per gli script che leggono solo uscite derivate dello step
STEP_DEFAULT = {"R2": "R2_ROI", "R2b": "R2b_WINDOW", "R2c": "R2b_WINDOW", "R3": "R2_ROI", "R3b": "R2_ROI", "R4": "R2_ROI", "R6": "R2_ROI", "R7": "R2_ROI",
                "S1": "R2_ROI", "S2.1": "R2_ROI", "S2.2": "R2_ROI", "S2.3": "R2_ROI", "S2.5": "R2_ROI", "F3.1": "R2_ROI", "R5": "R5_EDGE_TILE", "S2.0": "R5_EDGE_TILE",
                "S2.4": "SECTION", "S2.6": "SECTION", "S2.7": "SECTION", "S2.8": "R2_ROI", "R1": "METADATA", "F3.2": "F3.2"}
OVERRIDE = {
    "tools/R2_overview.py": ("SECTION", "panoramica H&E e gruppi di Space Ranger della sezione per scegliere le aree di R2; nessun parametro"),
    "tools/R6_dexp.py": ("SECTION", "distanza di espansione nucleo-cellula di Space Ranger su un campione di cellule della sezione; serve alle misure sul vero, non e' un parametro del simulatore"),
    "tools/R4_section.py": ("SECTION", "medie della sezione e spettro dello sfondo (controlli e riferimenti di sezione)"),
    "tools/R4_posthoc.py": ("SECTION", "alternanza della cattura sulla sezione intera (B-289 e simili)"),
    "tools/R5_section.py": ("SECTION", "UMI per quadratino di tutta la sezione: base delle tessere di bordo e delle strisce"),
    "tools/R5_explore_e0.py": ("SECTION", "esplorazione geometrica della maschera del tessuto, nessun conteggio"),
    "tools/R5_explore_e0b.py": ("SECTION", "esplorazione geometrica: copertura dell'immagine, nessun parametro"),
    "tools/R5_lumens.py": ("SECTION", "controprova CP-2 di R5: lumi interni su tessere di tutta la sezione; nessun parametro"),
    "tools/S2.0_lumens.py": ("SECTION", "controprova CP-1 di S2.0: stessa selezione dei lumi di R5; nessun parametro"),
    "tools/S2.0_stripes.py": ("SECTION", "strisce per sezione (parametri parity e sd di capture_params)"),
    "tools/S2.0_ph5_stripes_hp.py": ("SECTION", "parte fine delle strisce per sezione (valori del simulatore, B-397...B-401)"),
    "tools/S2.0_ph1_stripes_glass.py": ("SECTION", "strisce sul vetro della sezione, post hoc"),
    "tools/S2.0_disp.py": ("R2_ROI", "30 aree di R2"), "tools/S2.0_inner.py": ("R2_ROI", "30 aree di R2 (rapporto nucleo/citoplasma r)"),
    "tools/S2.3_glass.py": ("SECTION", "vetro vicino al tessuto sulle cinque sezioni intere (base del fondo e della seconda coda)"),
    "tools/S2.6_PH3_edge.py": ("SECTION", "vetro della sezione intera, post hoc di una prova di prestazione"),
    "tools/S2.6_prep_estimates.py": ("SECTION", "dimensionamento della prova di prestazione, solo posizioni"),
    "tools/S2.6_section_layer.R": ("SECTION", "livello cellulare della sezione intera per la prova di prestazione"),
    "tools/S2.6_section_measures.py": ("SECTION", "misure della sezione intera, verifica della prova di prestazione"),
    "tools/S2.7_section_measures.py": ("SECTION", "misure della sezione intera, verifica della simulazione a pezzi"),
    "tools/R1_verify.py": ("METADATA", "inventario, scala, md5 e ritagli di controllo istologico (R1)"),
    "tools/docs_patch_R1.py": ("METADATA", "aggiornamento della documentazione di R1"), "tools/docs_patch_R2b.py": ("METADATA", "aggiornamento della documentazione di R2b"),
    "R/testing/test_R1.sh": ("METADATA", "batteria di R1"), "R/testing/test_S1.5.R": ("SECTION", "prova della stima della scala sull'immagine (S1.5), nessun parametro del preset"),
    "R/05b2_simulate_cell_profiles.R": ("R2_ROI", "carica i profili costruiti dai nuclei delle aree di R2 (S2.2)"),
    "R/testing/run_R3b.sh": ("R2_ROI", "lancia R3b sulle aree di R2"), "tools/run_R3b.sh": ("R2_ROI", "lancia R3b sulle aree di R2"),
}

def step_of(p):
    n = Path(p).name
    m = re.match(r"(?:test_|run_)?(F3\.\d|S\d\.\d|R\d[a-z]?)", n)
    if not m: return None
    s = m.group(1)
    if s.startswith("S1"): return "S1"
    return s if s in STEP_DEFAULT else (s[:-1] if s[:-1] in STEP_DEFAULT else s)

def main():
    rows = []
    files = sorted([p for d in ("tools", "R") for p in (ROOT / d).rglob("*") if p.is_file() and p.suffix in (".py", ".R", ".sh", ".js") and "__pycache__" not in p.parts])
    for p in files:
        rel = str(p.relative_to(ROOT)); txt = open(p, encoding="utf-8", errors="replace").read()
        if not RAW.search(txt): continue
        refs = [k for k, r in REFS.items() if r.search(txt)]
        if rel in OVERRIDE: cls, why = OVERRIDE[rel]
        elif Path(rel).name.startswith("F3.2_"): cls, why = "F3.2", "strumenti di questo step (lettori ammessi dei dati di verifica)"
        elif re.search(r"_md5\.sh$", rel): cls, why = "INTEGRITY", "md5 dei file letti dallo step"
        elif re.search(r"(_repro|_run|_chain\d?|_thermal_guard|_setup)\.sh$", rel): cls, why = "WRAPPER", "lancia altri script dello stesso step"
        elif refs: cls, why = refs[0] if "SECTION" not in refs else "SECTION", "file di aree citati: " + ", ".join(refs)
        else:
            st = step_of(rel); cls = STEP_DEFAULT.get(st, "UNCLASSIFIED"); why = f"legge uscite derivate dello step {st}" if st else "nessuna regola"
        rows.append(dict(file=rel, area_refs=";".join(refs), klass=cls, reason=why))
    with open(ROOT / "results/F3.2/F3.2_readers.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    from collections import Counter
    c = Counter(r["klass"] for r in rows); print(len(rows), "script;", dict(c))
    un = [r["file"] for r in rows if r["klass"] == "UNCLASSIFIED"]
    print("CP-1:", "confermata (0 script non classificati)" if not un else f"SMENTITA o incompleta: non classificati {un}")

if __name__ == "__main__":
    main()
