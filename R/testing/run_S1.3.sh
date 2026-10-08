#!/usr/bin/env bash
# S1.3 — rifa' tutti i numeri del report di tessellate_voronoi() (dalla root del repo)
#   bash R/testing/run_S1.3.sh [outdir]   (default /mnt/micron/geo_spatialtrans/S1.3; stadi rilanciabili)
# Ordine: check C (geos, deldir) → mutanti → ROI (real, null, cp3) → CP-1/CP-2/D-1 → prestazioni (per ultime, macchina meno carica) → analisi
set -euo pipefail
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent
OUT=${1:-/mnt/micron/geo_spatialtrans/S1.3}; LOG=$OUT/logs; mkdir -p "$LOG"
for b in geos deldir; do
  S13_BACKEND=$b Rscript --vanilla R/testing/test_S1.3.R > "$LOG/test_$b.log" 2>&1
  tail -n 1 "$LOG/test_$b.log"
done
for m in M1 M2 M3 M4 M5 M6; do
  S13_BACKEND=geos S13_MUTANT=$m Rscript --vanilla R/testing/test_S1.3.R > "$LOG/test_geos_$m.log" 2>&1 || true
  tail -n 1 "$LOG/test_geos_$m.log"
done
for st in real null cp3; do
  Rscript --vanilla tools/S1.3_roi.R $st "$OUT" > "$LOG/roi_$st.log" 2>&1
  tail -n 1 "$LOG/roi_$st.log"
done
Rscript --vanilla tools/S1.3_cp.R all > "$LOG/cp.log" 2>&1; tail -n 3 "$LOG/cp.log"
P=results/S1.3/S1.3_perf.csv
echo "backend,n,elapsed_s,c2_rel,n_multipart,n_snapped,max_rss_mb,status" > "$P"
for b in geos deldir; do
  for n in 10000 20000 40000 80000 160000 449536; do
    lim=7200
    st=ok
    /usr/bin/time -v timeout $lim Rscript --vanilla tools/S1.3_perf.R $b $n > "$LOG/perf_${b}_$n.out" 2> "$LOG/perf_${b}_$n.time" || st=timeout
    line=$(grep '^PERF,' "$LOG/perf_${b}_$n.out" | sed 's/^PERF,//' || true)
    kb=$(grep 'Maximum resident set size' "$LOG/perf_${b}_$n.time" | awk '{print $NF}' || true)
    if [ -z "$line" ]; then line="$b,$n,NA,NA,NA,NA"; fi
    echo "$line,$(( ${kb:-0} / 1024 )),$st" >> "$P"
    tail -n 1 "$P"
  done
done
Rscript --vanilla tools/S1.3_analysis.R "$OUT" > "$LOG/analysis.log" 2>&1; tail -n 40 "$LOG/analysis.log"
