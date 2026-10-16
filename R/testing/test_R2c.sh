#!/usr/bin/env bash
# R/testing/test_R2c.sh — R2c: tutti i check C eseguibili con un comando (PASS/FAIL per asserzione). Dalla root del repo.
#   C-1 copertura (results/R2c/R2c_coverage.csv), C-2 annotatore (node + R2c_check_annotator.py), C-3 identita' con cp2,
#   C-4 mutanti dell'analisi (MU1-MU3 devono fallire), test delle funzioni dell'analisi.
set -u
cd "$(dirname "$0")/../.."
PY=.venv/bin/python; fail=0
ok(){ if [ "$1" = 0 ]; then echo "PASS $2"; else echo "FAIL $2"; fail=$((fail+1)); fi; }
$PY - <<'E'; ok $? "C-1 copertura pre-registrata PASS in tutti gli archetipi x segmentatore"
import pandas as pd, sys
c = pd.read_csv("results/R2c/R2c_coverage.csv", keep_default_na=False); c = c[c.role == "pre-registrato"]
print(c[["archetype", "seg", "coverage", "verdict"]].to_string(index=False)); sys.exit(0 if c.verdict.isin(["PASS", "NA"]).all() else 1)
E
$PY - <<'E'; ok $? "C-3 S-1 dei segmentatori = cp2 di S1.4c (48 ROI x segmentatore, 0 differenze)"
import pandas as pd, sys
c = pd.read_csv("results/R2c/R2c_c3.csv"); sys.exit(0 if (c.n_diff.sum() == 0 and c.labels_equal.all() and len(c) == 48) else 1)
E
(cd tools && node R2c_annotator_core.test.js > /tmp/r2c_node.log); ok $? "C-2 core JS ($(grep -c ^PASS /tmp/r2c_node.log) test)"
for r in Luca Federica; do $PY tools/R2c_check_annotator.py $r > /tmp/r2c_c2_$r.log; ok $? "C-2 annotatore $r ($(grep -c ^PASS /tmp/r2c_c2_$r.log) PASS)"; done
$PY tools/R2c_analysis_test.py > /tmp/r2c_ana.log 2>/dev/null; ok $? "analisi: test delle funzioni ($(grep -c ^PASS /tmp/r2c_ana.log) PASS)"
for m in MU1 MU2 MU3; do R2C_MUTANT=$m $PY tools/R2c_analysis_test.py > /dev/null 2>&1; [ $? -ne 0 ]; ok $? "C-4 mutante $m rilevato"; done
echo "== $fail FAIL"; exit $fail
