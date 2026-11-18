#!/usr/bin/env python3
"""tools/F3.2_guard.py - F3.2, cancello G-4 (c): nessun codice fuori dai lettori ammessi nomina i dati di verifica.
Cerca in R/ e tools/ (ricorsivo) i nomi dei dataset di verifica, le loro etichette e i percorsi data_holdout e F3.2_sealed.
Lettori ammessi: data/holdout/ALLOWED_READERS.tsv (file, step, scopo). Un lettore nuovo si aggiunge li', con lo step che lo usa.
Uso: python3 tools/F3.2_guard.py [cartella ...]  -> stampa le occorrenze non ammesse; esce con 1 se ce ne sono.
"""
import re, sys, csv
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
NAMES = ["data_holdout", "F3.2_sealed", "Visium_HD_Human_Colon_Cancer_P1", "Visium_HD_Human_Colon_Cancer_P5", "Visium_HD_Human_Colon_Normal_P3",
         "Visium_HD_Human_Colon_Normal_P5", "Visium_HD_6p5mm_Human_Colon_Cancer", "Visium_HD_Mouse_Brain", "Visium_HD_Human_Kidney_FFPE",
         "Visium_HD_6p5mm_Rat_Liver", "CRC_P1_SR300", "CRC_P5_SR300", "CRC_6p5mm_SR410", "MouseBrain_FFPE_SR401", "HumanKidney_FFPE_SR401",
         "RatLiver_6p5mm_SR410"]
PAT = re.compile(r"(?<![A-Za-z0-9_])(" + "|".join(re.escape(n) for n in NAMES) + r")(?![A-Za-z0-9])")
EXT = {".py", ".R", ".r", ".sh", ".js", ".qmd", ".Rmd"}

def allowed():
    f = ROOT / "data/holdout/ALLOWED_READERS.tsv"
    return {r["file"] for r in csv.DictReader(open(f), delimiter="\t")} if f.exists() else set()

def scan(dirs):
    ok = allowed(); hits = []
    for d in dirs:
        for p in sorted(Path(d).rglob("*")):
            if not p.is_file() or p.suffix not in EXT or "__pycache__" in p.parts: continue
            rel = str(p.resolve().relative_to(ROOT)) if p.resolve().is_relative_to(ROOT) else str(p)
            if rel in ok: continue
            for i, ln in enumerate(open(p, encoding="utf-8", errors="replace"), 1):
                m = PAT.search(ln)
                if m: hits.append((rel, i, m.group(1), ln.strip()[:120]))
    return hits

if __name__ == "__main__":
    dirs = sys.argv[1:] or [ROOT / "R", ROOT / "tools"]
    h = scan(dirs)
    for x in h: print("FUGA", *x)
    print(f"occorrenze non ammesse: {len(h)}"); sys.exit(1 if h else 0)
