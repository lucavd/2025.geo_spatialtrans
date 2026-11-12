#!/usr/bin/env bash
# R/testing/run_S2.8b.sh - S2.8b: cancelli G-2 (esempi del README e guida d'uso), G-3 (batterie nel clone nuovo), G-4 (unione) e controprova
# CP-2. Derivato da R/testing/run_S1.6b.sh. docs/steps/S2.8b.md. Uso (da qualsiasi cartella):
#   bash R/testing/run_S2.8b.sh clone <sha>   clone nuovo da GitHub in $S28B/clone al commit <sha>, libreria dal lockfile
#                                             (tools/restore_library.R), input dei test fuori da git collegati (tools/link_server_inputs.sh)
#   bash R/testing/run_S2.8b.sh quickstart    G-2 (a, b): i due blocchi dell'avvio rapido del README, due esecuzioni isolate ciascuno
#   bash R/testing/run_S2.8b.sh guide         G-2 (c): due rese della guida d'uso nel clone, impronte uguali, controlli TRUE
#   bash R/testing/run_S2.8b.sh cp2           CP-2: una copia della guida con piece_um = 500 non deve rendersi
#   bash R/testing/run_S2.8b.sh battery       G-3: batteria dello Step 1 (run_S1.6.sh test) e test_S2.1 ... test_S2.8 nel clone
#   bash R/testing/run_S2.8b.sh merge         G-4: dev = step2-visium-grid su GitHub, avanzamento da 9250749, main invariato, simulatore invariato
set -u
LIBREL=renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
S28B=${S28B:-/mnt/micron/geo_spatialtrans/S2.8b}; CL="$S28B/clone"; mkdir -p "$S28B"
URL=https://github.com/lucavd/2025.geo_spatialtrans.git
OLD_DEV=9250749f8dc43aae8555a523be2d69051dc3b312; S28_HEAD=48a5f18; S1_MD5=9f361ea9e1a56d8f46eef8a0a7fba842
GUIDE=vignettes/visium_hd_simulation.qmd
MAIN_AT_OPEN=${MAIN_AT_OPEN:-}
mode=${1:-}; fail=0
extract() {  # extract <file> <marker> -> blocco R fra <!-- marker:start --> e <!-- marker:end -->
  python3 - "$1" "$2" <<'PY'
import re, sys
t = open(sys.argv[1]).read(); k = sys.argv[2]
m = re.search(r"<!-- " + k + r":start -->\s*```r\n(.*?)```\s*<!-- " + k + r":end -->", t, re.S)
sys.stdout.write(m.group(1) if m else "")
PY
}
if [ "$mode" = clone ]; then
  sha=${2:?serve il commit}; rm -rf "$CL"; git clone -q "$URL" "$CL" && git -C "$CL" checkout -q "$sha" || { echo "FAIL clone"; exit 1; }
  echo "clone: $(git -C "$CL" rev-parse HEAD)"
  ( cd "$CL" && Rscript --vanilla tools/restore_library.R ) > "$S28B/restore.log" 2>&1 || fail=1
  tail -2 "$S28B/restore.log"; grep -q '^PASS restore' "$S28B/restore.log" || fail=1
  ( cd "$CL" && bash tools/link_server_inputs.sh ) > "$S28B/link.log" 2>&1 || fail=1; tail -1 "$S28B/link.log"
  [ $fail = 0 ] && echo "PASS clone e libreria" || echo "FAIL clone e libreria"
elif [ "$mode" = quickstart ]; then
  cd "$CL" || exit 1
  extract README.md quickstart > "$S28B/qs1.R"; extract README.md quickstart2 > "$S28B/qs2.R"
  for q in qs1 qs2; do [ -s "$S28B/$q.R" ] || { echo "FAIL blocco $q non trovato nel README"; fail=1; }
    echo "INFO blocco $q: $(wc -l < "$S28B/$q.R") righe, md5 $(md5sum < "$S28B/$q.R" | cut -c1-32)"; done
  for k in 1 2; do
    Rscript --vanilla -e "source('$S28B/qs1.R'); saveRDS(cl\$cell_df, '$S28B/qs1_cell_df_$k.rds', compress = FALSE)" > "$S28B/qs1_$k.log" 2>&1 || { echo "FAIL Step 1 esecuzione $k"; fail=1; }
    grep -q '^\[place_cells\].* 0 ERROR' "$S28B/qs1_$k.log" || { echo "FAIL Step 1 esecuzione $k: controlli con errori"; fail=1; }
    m=$(md5sum < "$S28B/qs1_cell_df_$k.rds" | cut -c1-32)
    [ "$m" = "$S1_MD5" ] && echo "PASS Step 1 esecuzione $k: md5 di cell_df = S1.6b ($m)" || { echo "FAIL Step 1 esecuzione $k: md5 $m"; fail=1; }
    # informativo (dopo la prima esecuzione di S2.8b): l'oggetto e' uguale a quello di S1.6b? byte diversi del file e versione di R che lo scrive
    Rscript --vanilla -e "a <- readRDS('/mnt/micron/geo_spatialtrans/S1.6b/official2/qs_cell_df_1.rds'); b <- readRDS('$S28B/qs1_cell_df_$k.rds'); cat('INFO Step 1 esecuzione $k: cell_df identico a S1.6b:', identical(a, b), '; R', R.version\$major, R.version\$minor, '\n')"
    echo "INFO Step 1 esecuzione $k: byte diversi dal file di S1.6b: $(cmp -l /mnt/micron/geo_spatialtrans/S1.6b/official2/qs_cell_df_1.rds "$S28B/qs1_cell_df_$k.rds" | wc -l) (posizioni: $(cmp -l /mnt/micron/geo_spatialtrans/S1.6b/official2/qs_cell_df_1.rds "$S28B/qs1_cell_df_$k.rds" | head -3 | awk '{print $1}' | tr '\n' ' '))"
    Rscript --vanilla -e "source('$S28B/qs2.R'); saveRDS(list(hd\$counts_2um, hd\$counts_8um), '$S28B/qs2_counts_$k.rds', compress = FALSE); cat('CHECKS', all(unlist(hd\$metadata\$checks)), '\n')" > "$S28B/qs2_$k.log" 2>&1 || { echo "FAIL Step 2 esecuzione $k (vedi qs2_$k.log)"; fail=1; }
    grep -q '^CHECKS TRUE' "$S28B/qs2_$k.log" || { echo "FAIL Step 2 esecuzione $k: controlli interni"; fail=1; }
  done
  a=$(md5sum < "$S28B/qs2_counts_1.rds" | cut -c1-32); b=$(md5sum < "$S28B/qs2_counts_2.rds" | cut -c1-32)
  [ "$a" = "$b" ] && echo "PASS Step 2: stesso md5 di counts_2um e counts_8um nelle due esecuzioni ($a)" || { echo "FAIL Step 2 md5 diversi: $a $b"; fail=1; }
  git status --porcelain --untracked-files=no | grep -q . && { echo "FAIL l'avvio rapido ha modificato file versionati"; fail=1; }
  [ $fail = 0 ] && echo "GATE G-2ab PASS" || echo "GATE G-2ab FAIL"
elif [ "$mode" = guide ]; then
  cd "$CL" || exit 1
  for k in 1 2; do
    # copia non versionata nella stessa cartella della guida: stessa cartella di lavoro degli esempi, risorse incorporate come nella resa
    # versionata (con --output Quarto scrive fuori da vignettes/ e non incorpora le risorse: prima esecuzione di S2.8b)
    cp "$GUIDE" "vignettes/g2_render_$k.qmd"
    R_LIBS="$CL/$LIBREL" quarto render "vignettes/g2_render_$k.qmd" > "$S28B/guide_$k.log" 2>&1 || { echo "FAIL resa $k (vedi guide_$k.log)"; fail=1; }
    [ -s "vignettes/g2_render_$k.html" ] || { echo "FAIL resa $k: vignettes/g2_render_$k.html assente"; fail=1; }
    grep -oE '(FINGERPRINT|CHECK) [A-Za-z0-9_]+ [A-Za-z0-9]+' "vignettes/g2_render_$k.html" | sort > "$S28B/guide_fp_$k.txt"
    echo "INFO resa $k: $(grep -c FINGERPRINT "$S28B/guide_fp_$k.txt") impronte, $(grep -c CHECK "$S28B/guide_fp_$k.txt") controlli"
    cp "vignettes/g2_render_$k.html" "$S28B/"; rm -rf "vignettes/g2_render_$k.qmd" "vignettes/g2_render_$k.html" "vignettes/g2_render_${k}_files"
  done
  [ -s "$S28B/guide_fp_1.txt" ] && cmp -s "$S28B/guide_fp_1.txt" "$S28B/guide_fp_2.txt" && echo "PASS impronte e controlli uguali nelle due rese" || { echo "FAIL impronte diverse o assenti"; fail=1; }
  [ "$(grep -c '^FINGERPRINT' "$S28B/guide_fp_1.txt")" -ge 5 ] || { echo "FAIL meno di 5 impronte"; fail=1; }
  [ "$(grep -c '^CHECK' "$S28B/guide_fp_1.txt")" -ge 5 ] || { echo "FAIL meno di 5 controlli"; fail=1; }
  # ogni controllo deve finire con TRUE (revisione T05: un NA non passa)
  grep '^CHECK' "$S28B/guide_fp_1.txt" | grep -vq ' TRUE$' && { echo "FAIL controlli non TRUE:"; grep '^CHECK' "$S28B/guide_fp_1.txt" | grep -v ' TRUE$'; fail=1; }
  # la guida resa nel commit deve mostrare le stesse impronte (revisione T06: era solo informativo)
  grep -oE '(FINGERPRINT|CHECK) [A-Za-z0-9_]+ [A-Za-z0-9]+' "${GUIDE%.qmd}.html" | sort > "$S28B/guide_fp_committed.txt"
  cmp -s "$S28B/guide_fp_1.txt" "$S28B/guide_fp_committed.txt" && echo "PASS impronte uguali alla guida resa nel commit" || { echo "FAIL impronte diverse dalla guida resa nel commit"; fail=1; }
  git status --porcelain --untracked-files=no | grep -q . && { echo "FAIL la resa ha modificato file versionati"; fail=1; }
  [ $fail = 0 ] && echo "GATE G-2c PASS" || echo "GATE G-2c FAIL"
elif [ "$mode" = cp2 ]; then
  cd "$CL" || exit 1
  sed '0,/piece_um = Inf)/s//piece_um = 500)/' "$GUIDE" > vignettes/cp2_guide.qmd
  grep -c 'piece_um = 500)' vignettes/cp2_guide.qmd
  if R_LIBS="$CL/$LIBREL" quarto render vignettes/cp2_guide.qmd > "$S28B/cp2.log" 2>&1; then echo "FAIL CP-2: la guida con piece_um = 500 si e' resa"; fail=1
  elif grep -q 'multiplo di 56' "$S28B/cp2.log"; then echo "PASS CP-2: la resa si ferma sull'errore del lato dei pezzi (multiplo di 56)"
  else echo "FAIL CP-2: la resa si ferma, ma non per il lato dei pezzi (vedi cp2.log)"; fail=1; fi   # revisione T07
  rm -f vignettes/cp2_guide.qmd vignettes/cp2_guide.html; rm -rf vignettes/cp2_guide_files
elif [ "$mode" = battery ]; then
  cd "$CL" || exit 1
  ( S16_LOG="$S28B/battery1" bash R/testing/run_S1.6.sh test > "$S28B/battery1.log" 2>&1 ) &
  declare -A EXP=( [2.1]=21 [2.2]=36 [2.3]=39 [2.4]=37 [2.5]=40 [2.6]=33 [2.7]=78 [2.8]=46 )
  for s in 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8; do
    ( Rscript --vanilla "R/testing/test_S$s.R" > "$S28B/test_S$s.log" 2>&1; echo "EXIT $?" >> "$S28B/test_S$s.log" ) &
  done
  wait
  grep -E '^(GATE|PASS nessun)' "$S28B/battery1.log" | tail -4
  grep -q '^GATE G-1 PASS' "$S28B/battery1.log" && echo "PASS batteria dello Step 1" || { echo "FAIL batteria dello Step 1"; fail=1; }
  for s in 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8; do
    np=$(grep -cE '^PASS' "$S28B/test_S$s.log" || true); nf=$(grep -cE '^FAIL' "$S28B/test_S$s.log" || true); ex=$(grep -oE '^EXIT [0-9]+' "$S28B/test_S$s.log" | cut -d' ' -f2)
    if [ "$np" = "${EXP[$s]}" ] && [ "$nf" = 0 ] && [ "$ex" = 0 ]; then echo "PASS test_S$s: $np PASS, 0 FAIL (atteso ${EXP[$s]})"
    else echo "FAIL test_S$s: $np PASS, $nf FAIL, uscita $ex (atteso ${EXP[$s]})"; fail=1; fi
  done
  git status --porcelain --untracked-files=no | grep -q . && { echo "FAIL le batterie hanno modificato file versionati"; fail=1; }
  [ $fail = 0 ] && echo "GATE G-3 PASS" || echo "GATE G-3 FAIL"
elif [ "$mode" = merge ]; then
  R=$(cd "$(dirname "$0")/../.." && pwd); cd "$R" && git fetch -q origin
  d=$(git ls-remote origin refs/heads/dev | cut -f1); s=$(git ls-remote origin refs/heads/step2-visium-grid | cut -f1); m=$(git ls-remote origin refs/heads/main | cut -f1)
  [ -n "$d" ] && [ "$d" = "$s" ] && echo "PASS dev = step2-visium-grid su GitHub ($d)" || { echo "FAIL dev $d, step2-visium-grid $s"; fail=1; }
  git merge-base --is-ancestor "$OLD_DEV" "$d" && echo "PASS il vecchio dev ($OLD_DEV) e' antenato di dev (avanzamento)" || { echo "FAIL vecchio dev non antenato"; fail=1; }
  if [ -n "$MAIN_AT_OPEN" ]; then [ "$m" = "$MAIN_AT_OPEN" ] && echo "PASS main invariato ($m)" || { echo "FAIL main cambiato: $m"; fail=1; }; else echo "FAIL MAIN_AT_OPEN non dato"; fail=1; fi
  ch=$(git diff --name-only "$S28_HEAD" "$d" -- R/ inst/ renv.lock renv/activate.R renv/settings.json .Rprofile tools/restore_library.R tools/link_server_inputs.sh | grep -v '^R/testing/run_S2.8b.sh$' || true)
  [ -z "$ch" ] && echo "PASS R/, inst/, renv.lock, renv/, .Rprofile, strumenti d'installazione identici a S2.8 ($S28_HEAD), salvo R/testing/run_S2.8b.sh" || { echo "FAIL file del simulatore o dell'ambiente cambiati:"; echo "$ch"; fail=1; }
  v=$(git show "$d:DESCRIPTION" | sed -n 's/^Version: *//p'); [ "$v" = "0.2.0" ] && echo "PASS Version 0.2.0" || { echo "FAIL Version $v"; fail=1; }   # revisione T08
  python3 - "$S28_HEAD" "$d" <<'PY' || fail=1
import subprocess, sys, re
def dcf(rev):
    t = subprocess.run(["git", "show", f"{rev}:DESCRIPTION"], capture_output=True, text=True, check=True).stdout
    out, k = {}, None
    for line in t.splitlines():
        if line[:1] in (" ", "\t") and k: out[k] += "\n" + line
        elif ":" in line: k, v = line.split(":", 1); out[k] = v
    return out
a, b = dcf(sys.argv[1]), dcf(sys.argv[2])
diff = sorted(k for k in set(a) | set(b) if a.get(k) != b.get(k))
ok = set(diff) <= {"Version", "Date", "Description"}
print(f"{'PASS' if ok else 'FAIL'} DESCRIPTION cambiato solo in Version, Date, Description (campi diversi: {', '.join(diff) or 'nessuno'})")
sys.exit(0 if ok else 1)
PY
  if [ -n "${TESTED:-}" ]; then
    after=$(git diff --name-only "$TESTED" "$d" | grep -vE '^(docs/|reports/|CHANGELOG.md$)' || true)
    [ -z "$after" ] && echo "PASS dopo il commit provato ($TESTED) solo docs/, reports/, CHANGELOG.md" || { echo "FAIL dopo il commit provato cambiano:"; echo "$after"; fail=1; }
  else echo "FAIL TESTED non dato"; fail=1; fi
  [ $fail = 0 ] && echo "GATE G-4 PASS" || echo "GATE G-4 FAIL"
else echo "uso: run_S2.8b.sh clone <sha> | quickstart | guide | cp2 | battery | merge"; exit 2; fi
exit $fail
