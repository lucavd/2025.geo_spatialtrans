#!/usr/bin/env bash
# tools/S2.0_repro.sh <run1> <run2> — S2.0, G-4: confronto fra due esecuzioni complete con lo stesso seme.
# (1) md5 delle tabelle di risultato (results/S2.0 contro <run2>/results); (2) file per unita' (finestre, controllo interno, lumi,
# dispersione, strisce) identici a meno del solo campo dei tempi "sec". Esce con 0 se tutto coincide.
set -u
R1=$1; R2=$2; RES1=${S20_RES1:-results/S2.0}; RES2=${S20_RES2:-$R2/results}; FAIL=0
for f in S2.0_edges_tiles.csv S2.0_edges_sections.csv S2.0_inner_roi.csv S2.0_inner_summary.csv S2.0_RA.json S2.0_RB.json S2.0_RC.json S2.0_predictions.csv S2.0_CP1_lumens.csv S2.0_disp_roi.csv; do
  a=$(md5sum "$RES1/$f" 2>/dev/null | cut -c1-32); b=$(md5sum "$RES2/$f" 2>/dev/null | cut -c1-32)
  if [ -n "$a" ] && [ "$a" = "$b" ]; then echo "PASS table $f $a"; else echo "FAIL table $f run1=$a run2=$b"; FAIL=1; fi
done
.venv/bin/python - "$R1" "$R2" <<'EOF' || FAIL=1
import sys, os, json, glob, hashlib
import numpy as np
R1, R2 = sys.argv[1], sys.argv[2]; bad = 0; n = 0
def strip(o):
    if isinstance(o, dict): return {k: strip(v) for k, v in o.items() if k != "sec"}
    if isinstance(o, list): return [strip(v) for v in o]
    return o
for sub in ("edges", "inner", "disp", "stripes"):
    for f1 in sorted(glob.glob(f"{R1}/{sub}/*.json")):
        if os.path.basename(f1).startswith(("PH1_", "e3")): continue      # analisi post hoc ed E3: non fanno parte della catena
        f2 = f1.replace(R1, R2, 1); n += 1
        if not os.path.exists(f2): print("FAIL missing", f2); bad += 1; continue
        if strip(json.load(open(f1))) != strip(json.load(open(f2))): print("FAIL differs", f1); bad += 1
for f1 in sorted(glob.glob(f"{R1}/lumens/*.csv")) + sorted(glob.glob(f"{R1}/stripes/*.npz")):
    if os.path.basename(f1).startswith(("PH1_", "e3")): continue
    f2 = f1.replace(R1, R2, 1); n += 1
    if not os.path.exists(f2): print("FAIL missing", f2); bad += 1; continue
    if f1.endswith(".csv"):
        same = open(f1, "rb").read() == open(f2, "rb").read()
    else:
        a, b = np.load(f1), np.load(f2); same = set(a.files) == set(b.files) and all(np.array_equal(a[k], b[k], equal_nan=True) for k in a.files)
    if not same: print("FAIL differs", f1); bad += 1
print(f"unit files compared {n}, differing or missing {bad}")
sys.exit(1 if bad else 0)
EOF
[ $FAIL -eq 0 ] && echo "G-4 PASS" || echo "G-4 FAIL"
