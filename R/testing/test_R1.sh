#!/usr/bin/env bash
# test_R1.sh — un comando per tutte le verifiche R1 (PASS/FAIL per asserzione).
#   bash R/testing/test_R1.sh            # completo: md5 di 80 GB (~15 min) + verifiche
#   R1_FAST=1 bash R/testing/test_R1.sh  # salta il ricalcolo md5
set -uo pipefail
cd "$(dirname "$0")/../.."
OUT=results/R1; mkdir -p "$OUT" logs
LOG=logs/R1_test_$(date -u +%Y%m%dT%H%M%SZ).log
fail=0
echo "== R1 test — $(date -u +%FT%TZ) — commit $(git rev-parse --short HEAD 2>/dev/null)" | tee "$LOG"

# 1. dati presenti
[ -L data/real/datasets ] && [ -d data/real/datasets ] && echo "PASS symlink data/real/datasets -> $(readlink -f data/real/datasets)" | tee -a "$LOG" || { echo "FAIL symlink data/real/datasets" | tee -a "$LOG"; fail=1; }

# 2. md5 (checksums.md5 tracciato in git = md5 pubblicati da 10x e ricalcolati in R1)
if [ "${R1_FAST:-0}" = "1" ]; then
  echo "SKIP md5sum -c (R1_FAST=1)" | tee -a "$LOG"; MD5FLAG=--skip-md5
else
  if (cd data/real && md5sum -c --quiet checksums.md5) >>"$LOG" 2>&1; then echo "PASS md5sum -c checksums.md5 (38 file)" | tee -a "$LOG"; else echo "FAIL md5sum -c checksums.md5 (vedi $LOG)" | tee -a "$LOG"; fail=1; fi
  MD5FLAG=--skip-md5   # già verificato sopra; evita il doppio ricalcolo
fi

# 3. verifiche C-R1.2, B-R1.3, CP-1, CP-2 (ritagli)
[ -x .venv/bin/python ] || bash tools/setup_python_env.sh >>"$LOG" 2>&1
.venv/bin/python tools/R1_verify.py --data data/real/datasets --manifest data/real/inventory/R1_download_manifest.csv --out "$OUT" $MD5FLAG >>"$LOG" 2>&1
python3 - "$OUT/R1_checks.csv" <<'EOF' | tee -a "$LOG"
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
bad = 0
for r in rows:
    if r["check"] == "C-R1.1" and r["value"].startswith("0/"): continue   # md5 saltato: coperto dal passo 2
    st = r["status"]; print(f"{st:7s} {r['check']:8s} {r['dataset_id']:36s} {r['value'][:90]}")
    bad += st == "FAIL"
print(f"RIEPILOGO verifiche: {len(rows)} righe, FAIL={bad}, PENDING(verifica visiva)={sum(r['status']=='PENDING' for r in rows)}")
sys.exit(1 if bad else 0)
EOF
[ ${PIPESTATUS[0]} -eq 0 ] || fail=1

# 4. inventario tracciato
for f in data/real/README.md data/real/checksums.md5 data/real/inventory/R1_download_manifest.csv data/real/inventory/R1_archetype_map.csv; do
  [ -s "$f" ] && echo "PASS presente $f" | tee -a "$LOG" || { echo "FAIL manca $f" | tee -a "$LOG"; fail=1; }
done
echo "== esito: $([ $fail -eq 0 ] && echo PASS || echo FAIL) — log $LOG" | tee -a "$LOG"
exit $fail
