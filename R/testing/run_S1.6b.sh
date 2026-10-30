#!/usr/bin/env bash
# R/testing/run_S1.6b.sh - S1.6b: cancelli G-2 (avvio rapido del README), G-3 (batteria nel clone nuovo), G-4 (merge). docs/steps/S1.6b.md.
# Uso (da qualsiasi cartella):
#   bash R/testing/run_S1.6b.sh clone <sha>   clone nuovo da GitHub in $S16B/clone al commit <sha>, libreria dal lockfile (tools/restore_library.R),
#                                             input dei test fuori da git collegati alla copia di lavoro del server (tools/link_server_inputs.sh)
#   bash R/testing/run_S1.6b.sh quickstart    G-2: blocco dell'avvio rapido estratto dal README del clone, due esecuzioni isolate, stesso md5 di cell_df
#   bash R/testing/run_S1.6b.sh battery       G-3: bash R/testing/run_S1.6.sh test nel clone (FAIL ammessi come in S1.6)
#   bash R/testing/run_S1.6b.sh merge         G-4: dev = step1-cell-layer su GitHub, fast-forward da 7836829, main invariato, simulatore invariato
set -u
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
S16B=${S16B:-/mnt/micron/geo_spatialtrans/S1.6b}; CL="$S16B/clone"; mkdir -p "$S16B"
URL=https://github.com/lucavd/2025.geo_spatialtrans.git
OLD_DEV=7836829c850cb3c1dd2f512470e52d7a851344f3; S16_HEAD=775eec3
MAIN_AT_OPEN=${MAIN_AT_OPEN:-}
mode=${1:-}; fail=0
if [ "$mode" = clone ]; then
  sha=${2:?serve il commit}; rm -rf "$CL"; git clone -q "$URL" "$CL" && git -C "$CL" checkout -q "$sha" || { echo "FAIL clone"; exit 1; }
  echo "clone: $(git -C "$CL" rev-parse HEAD)"
  ( cd "$CL" && Rscript --vanilla tools/restore_library.R ) > "$S16B/restore.log" 2>&1 || fail=1
  tail -3 "$S16B/restore.log"; grep -q '^PASS restore' "$S16B/restore.log" || fail=1
  ( cd "$CL" && bash tools/link_server_inputs.sh ) > "$S16B/link.log" 2>&1 || fail=1; tail -1 "$S16B/link.log"   # passo del README (S1.6b, esplorazione E1)
  [ $fail = 0 ] && echo "PASS clone e libreria" || echo "FAIL clone e libreria"
elif [ "$mode" = quickstart ]; then
  cd "$CL" || exit 1
  python3 - README.md > "$S16B/quickstart.R" <<'PY'
import re, sys
t = open(sys.argv[1]).read()
m = re.search(r"<!-- quickstart:start -->\s*```r\n(.*?)```\s*<!-- quickstart:end -->", t, re.S)
sys.stdout.write(m.group(1) if m else "")
PY
  [ -s "$S16B/quickstart.R" ] || { echo "FAIL blocco dell'avvio rapido non trovato nel README"; exit 1; }
  echo "INFO blocco estratto: $(wc -l < "$S16B/quickstart.R") righe, md5 $(md5sum < "$S16B/quickstart.R" | cut -c1-32)"
  for k in 1 2; do
    Rscript --vanilla -e "source('$S16B/quickstart.R'); saveRDS(cl\$cell_df, '$S16B/qs_cell_df_$k.rds', compress = FALSE)" > "$S16B/quickstart_$k.log" 2>&1 || { echo "FAIL esecuzione $k (vedi quickstart_$k.log)"; fail=1; }
    grep '^\[place_cells\]' "$S16B/quickstart_$k.log"
    grep -q '^\[place_cells\].* 0 ERROR' "$S16B/quickstart_$k.log" || { echo "FAIL esecuzione $k: controlli con errori"; fail=1; }
  done
  m1=$(md5sum < "$S16B/qs_cell_df_1.rds" | cut -c1-32); m2=$(md5sum < "$S16B/qs_cell_df_2.rds" | cut -c1-32)
  [ "$m1" = "$m2" ] && echo "PASS stesso md5 di cell_df nelle due esecuzioni ($m1)" || { echo "FAIL md5 diversi: $m1 $m2"; fail=1; }
  git status --porcelain --untracked-files=no | grep -q . && { echo "FAIL l'avvio rapido ha modificato file versionati"; fail=1; }
  [ $fail = 0 ] && echo "GATE G-2 PASS" || echo "GATE G-2 FAIL"
elif [ "$mode" = battery ]; then
  cd "$CL" || exit 1
  S16_LOG="$S16B/battery" bash R/testing/run_S1.6.sh test > "$S16B/battery.log" 2>&1
  grep -E '^(GATE|FAIL|PASS nessun)' "$S16B/battery.log" | tail -12
  grep -q '^GATE G-1 PASS' "$S16B/battery.log" && echo "GATE G-3 PASS" || echo "GATE G-3 FAIL"
elif [ "$mode" = merge ]; then
  R=$(cd "$(dirname "$0")/../.." && pwd); cd "$R" && git fetch -q origin
  d=$(git ls-remote origin refs/heads/dev | cut -f1); s=$(git ls-remote origin refs/heads/step1-cell-layer | cut -f1); m=$(git ls-remote origin refs/heads/main | cut -f1)
  [ -n "$d" ] && [ "$d" = "$s" ] && echo "PASS dev = step1-cell-layer su GitHub ($d)" || { echo "FAIL dev $d, step1-cell-layer $s"; fail=1; }
  git merge-base --is-ancestor "$OLD_DEV" "$d" && echo "PASS il vecchio dev ($OLD_DEV) e' antenato di dev (fast-forward)" || { echo "FAIL vecchio dev non antenato"; fail=1; }
  if [ -n "$MAIN_AT_OPEN" ]; then [ "$m" = "$MAIN_AT_OPEN" ] && echo "PASS main invariato ($m)" || { echo "FAIL main cambiato: $m"; fail=1; }; else echo "FAIL MAIN_AT_OPEN non dato"; fail=1; fi
  ch=$(git diff --name-only "$S16_HEAD" "$d" -- R/ inst/ renv.lock | grep -v '^R/testing/run_S1.6b.sh$' || true)
  [ -z "$ch" ] && echo "PASS R/, inst/, renv.lock identici a S1.6 ($S16_HEAD), salvo R/testing/run_S1.6b.sh" || { echo "FAIL file del simulatore cambiati:"; echo "$ch"; fail=1; }
  if [ -n "${TESTED:-}" ]; then
    after=$(git diff --name-only "$TESTED" "$d" | grep -vE '^(docs/|reports/|CHANGELOG.md$)' || true)
    [ -z "$after" ] && echo "PASS dopo il commit provato ($TESTED) solo docs/, reports/, CHANGELOG.md" || { echo "FAIL dopo il commit provato cambiano:"; echo "$after"; fail=1; }
  else echo "FAIL TESTED non dato"; fail=1; fi
  [ $fail = 0 ] && echo "GATE G-4 PASS" || echo "GATE G-4 FAIL"
else echo "uso: run_S1.6b.sh clone <sha> | quickstart | battery | merge"; exit 2; fi
exit $fail
