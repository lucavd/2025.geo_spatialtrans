#!/usr/bin/env bash
# tools/run_S1.3d.sh — S1.3d (BL-084): run ufficiale, stadi in ordine e rilanciabili (ogni stadio salta ciò che è già scritto).
# Uso (root del repo, dal commit del run, senza modificare gli script durante l'esecuzione): bash tools/run_S1.3d.sh <fase> [outdir]
#   fase = main (tutto tranne H) | perf | b1 | analysis | holdout
set -u
F=${1:-main}; O=${2:-/mnt/micron/geo_spatialtrans/S1.3d}; mkdir -p "$O/logs" "$O/checks"
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent PYTHONDONTWRITEBYTECODE=1 S13D_CORES=${S13D_CORES:-40}
PY=.venv/bin/python; RS="Rscript --vanilla"; L="$O/logs"
step() { echo "[$(date '+%F %T')] $*" | tee -a "$L/run.log"; }
if [ "$F" = main ]; then
  step "md5 ingressi R3";            $RS tools/S1.3d_run.R md5 - "$O" >> "$L/md5.log" 2>&1 || exit 1
  step "C-1 test";                   $RS R/testing/test_S1.3d.R > "$L/test.log" 2>&1 || exit 1
  step "test exactgeom";             $PY results/S1.3b/review/exact/test_exactgeom.py > "$L/test_exactgeom.log" 2>&1 || exit 1
  step "export N";                   $RS tools/S1.3d_run.R export N "$O" > "$L/export_N.log" 2>&1 || exit 1
  for S in U CP1 S; do step "export $S"; $RS tools/S1.3d_run.R export_built $S "$O" > "$L/export_$S.log" 2>&1 || exit 1; done
  $RS tools/S1.3d_run.R list N "$O" > "$O/keys_N.txt"; $RS tools/S1.3d_run.R list O2 "$O" > "$O/keys_O2.txt"; $RS tools/S1.3d_run.R list K "$O" > "$O/keys_K.txt"
  grep -E '^(U1|U3|U4|U5)' "$O/keys_U.txt" > "$O/keys_U_o2.txt"; cat "$O/keys_O2.txt" "$O/keys_U_o2.txt" "$O/keys_CP1.txt" > "$O/keys_o2_all.txt"
  step "arbitro N, U, CP1";          for S in N U CP1; do $PY tools/S1.3d_arbiter_batch.py "$O" "$O/keys_$S.txt" --o2-keys "$O/keys_o2_all.txt" --cores 40 > "$L/arb_$S.log" 2>&1 || exit 1; done
  step "arbitro S (40 blocchi)";     seq 0 39 | xargs -P 40 -I{} $PY tools/S1.3d_arbiter.py "$O/geom" "$O/truth" S_cperf --o2 --chunk {}/40 > "$L/arb_S.log" 2>&1 || exit 1
  $PY tools/S1.3d_arbiter_merge.py "$O/truth" S_cperf 40 >> "$L/arb_S.log" 2>&1 || exit 1
  step "candidati N";                $RS tools/S1.3d_run.R cands N "$O" > "$L/cands_N.log" 2>&1 || exit 1
  for S in U CP1; do step "candidati $S"; $RS tools/S1.3d_run.R cands_built $S "$O" > "$L/cands_$S.log" 2>&1 || exit 1; done
  step "copertura N, U, CP1";        for S in N U CP1; do $PY tools/S1.3d_coverage.py "$O" "$O/keys_$S.txt" --cores 40 > "$L/cov_$S.log" 2>&1 || exit 1; done
  step "mutanti";                    $RS tools/S1.3d_run.R mutants A3_r3_CSR_03 "$O" > "$L/mutants.log" 2>&1 || exit 1
  step "C-2 / CP-4 arbitro";         shuf -n 25 --random-source=<(yes 20261014) "$O/keys_N.txt" > "$O/keys_chk_rand.txt"
  cat "$O/keys_K.txt" "$O/keys_chk_rand.txt" > "$O/keys_chk_n.txt"
  $PY tools/S1.3d_arbiter_check.py "$O" "$O/keys_chk_n.txt" 20 "$O/checks/arbiter_check_N.json" --cores 40 > "$L/chk_N.log" 2>&1 || exit 1
  $PY tools/S1.3d_arbiter_check.py "$O" "$O/keys_U.txt" 10 "$O/checks/arbiter_check_U.json" --cores 40 > "$L/chk_U.log" 2>&1 || exit 1
  step "main fatto"
fi
if [ "$F" = perf ]; then step "prestazioni S"; bash tools/S1.3d_perf.sh "$O" 4 > "$L/perf_driver.log" 2>&1
  echo S_cperf > "$O/keys_S.txt"; $PY tools/S1.3d_coverage.py "$O" "$O/keys_S.txt" --cores 1 > "$L/cov_S.log" 2>&1; step "perf fatto"; fi
if [ "$F" = b1 ]; then step "B-1"; $RS tools/S1.3d_b1.R "$O" > "$L/b1.log" 2>&1; step "b1 fatto"; fi
if [ "$F" = analysis ]; then step "analisi"; $RS tools/S1.3d_analysis.R "$O" > "$L/analysis.log" 2>&1; step "analisi fatta"; fi
if [ "$F" = holdout ]; then
  step "export H"; $RS tools/S1.3d_run.R export H "$O" > "$L/export_H.log" 2>&1 || exit 1; $RS tools/S1.3d_run.R list H "$O" > "$O/keys_H.txt"
  step "arbitro H"; $PY tools/S1.3d_arbiter_batch.py "$O" "$O/keys_H.txt" --o2-keys "$O/keys_H.txt" --cores 40 > "$L/arb_H.log" 2>&1 || exit 1
  step "candidati H"; S13D_O2=1 $RS tools/S1.3d_run.R cands H "$O" > "$L/cands_H.log" 2>&1 || exit 1
  $PY tools/S1.3d_coverage.py "$O" "$O/keys_H.txt" --cores 40 > "$L/cov_H.log" 2>&1
  step "analisi con H"; $RS tools/S1.3d_analysis.R "$O" --holdout > "$L/analysis_holdout.log" 2>&1; step "holdout fatto"; fi
