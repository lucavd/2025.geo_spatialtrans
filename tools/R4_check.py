"""tools/R4_check.py — cancello G-3 di R4 (comando pre-registrato; revisione F04). Uso: .venv/bin/python tools/R4_check.py [outdir]
Legge i controlli scritti dagli stadi (G-3a e G-3b in R4_section.py, G-3c in R4_roi.py) e stampa PASS/FAIL per sezione e per ROI."""
import sys, json, glob
OUT = sys.argv[1] if len(sys.argv) > 1 else "/mnt/micron/geo_spatialtrans/R4/run1"
n = fail = 0
for f in sorted(glob.glob(f"{OUT}/section/*.json")):
    S = json.load(open(f))
    for sc, v in S["G3a"].items():
        n += 1; fail += not v["pass"]
        print(f"{'PASS' if v['pass'] else 'FAIL'}  G-3a {S['dataset']} {sc}: n bin {v['ours']['n_bins']} = {v['ref']['n_bins']}, diff. rel. UMI {v['rel_diff_umi']:.1e}, geni {v['rel_diff_genes']:.1e}")
    g = S["G3b"]; n += 1; fail += not g["pass"]
    print(f"{'PASS' if g['pass'] else 'FAIL'}  G-3b {S['dataset']}: bin a 8 um diversi dalla somma dei 16 bin a 2 um = {g['n_diff']} su {g['n_bins8']}")
for f in sorted(glob.glob(f"{OUT}/roi/*_summary.json")):
    S = json.load(open(f)); ok = all(v["matrix_identical"] and v["table_identical"] for v in S["G3c"].values()); n += 1; fail += not ok
    print(f"{'PASS' if ok else 'FAIL'}  G-3c {S['archetype']}_{S['roi_id']}: impronte nucleo, nucleo + corona, maschera spostata = R6")
print(f"R4_check: {n - fail}/{n} PASS"); sys.exit(1 if fail else 0)
