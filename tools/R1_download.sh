#!/usr/bin/env bash
# R1_download.sh — download riprendibile dei dataset Visium HD selezionati (manifest CSV) su /mnt/micron
# Uso: bash R1_download.sh <manifest.csv> <dest_root> <status.tsv>
set -uo pipefail
MAN="$1"; ROOT="$2"; STATUS="$3"
mkdir -p "$ROOT" "$(dirname "$STATUS")"
[ -f "$STATUS" ] || printf "archetype\tdataset_id\tfilename\texpected_md5\tobserved_md5\tbytes\tstatus\tattempts\tfinished_utc\n" > "$STATUS"
python3 - "$MAN" <<'EOF' > /tmp/R1_manifest_rows.tsv
import csv,sys
seen=set()
for r in csv.DictReader(open(sys.argv[1])):
    if r['url'] in seen: continue
    seen.add(r['url'])
    print("\t".join([r['archetype'],r['dataset_id'],r['filename'],r['md5'],r['url']]))
EOF
n_ok=0; n_fail=0
while IFS=$'\t' read -r ARCH DS FN MD5 URL; do
  D="$ROOT/$DS"; mkdir -p "$D"; F="$D/$FN"
  # già verificato in un run precedente?
  if grep -qP "\t${FN}\t.*\tOK\t" "$STATUS" 2>/dev/null && [ -f "$F" ]; then echo "SKIP $FN (già OK)"; n_ok=$((n_ok+1)); continue; fi
  ok=0
  for attempt in 1 2 3; do
    echo "[$(date -u +%FT%TZ)] GET $FN (tentativo $attempt)"
    curl -sS -L --retry 5 --retry-delay 10 -C - -o "$F" "$URL" ; rc=$?
    if [ $rc -ne 0 ] && [ $rc -ne 33 ]; then echo "curl rc=$rc"; sleep 20; continue; fi   # 33 = range non supportato (file già completo)
    OBS=$(md5sum "$F" | cut -d' ' -f1)
    if [ "$OBS" = "$MD5" ]; then ok=1; break; else echo "MD5 MISMATCH $FN: atteso $MD5 osservato $OBS — riscarico"; rm -f "$F"; fi
  done
  BYTES=$( [ -f "$F" ] && stat -c %s "$F" || echo 0 )
  if [ $ok -eq 1 ]; then ST=OK; n_ok=$((n_ok+1)); else ST=FAIL; n_fail=$((n_fail+1)); OBS=${OBS:-NA}; fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$ARCH" "$DS" "$FN" "$MD5" "$OBS" "$BYTES" "$ST" "$attempt" "$(date -u +%FT%TZ)" >> "$STATUS"
  echo "[$(date -u +%FT%TZ)] $ST $FN ($BYTES B)"
done < /tmp/R1_manifest_rows.tsv
echo "RIEPILOGO: OK=$n_ok FAIL=$n_fail"
[ $n_fail -eq 0 ]
