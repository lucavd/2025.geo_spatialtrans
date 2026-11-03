#!/usr/bin/env bash
# tools/S2.0_chain.sh <rundir> [NW] — S2.0: catena completa di un'esecuzione (M2 finestre con due esecutori, riepilogo, controllo interno, R-A,
# lumi con C0 e con la regola vincente, dispersione con la regola vincente; M1 strisce se mancano). Rilanciabile (un file per unita').
# Le tabelle di analisi vanno in $S20_RES (default results/S2.0). Pre-registrazione: docs/steps/S2.0.md.
set -u
cd "$(dirname "$0")/.."
RUN=$1; NW=${2:-24}; H=$((NW / 2))
PY=.venv/bin/python; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p "$RUN"; LOG=$RUN/chain.log; echo "$(date -Is) start commit $(git rev-parse --short HEAD) NW=$NW RES=${S20_RES:-results/S2.0}" >> $LOG
[ -f "$RUN/stripes/A6.json" ] || { $PY tools/S2.0_stripes.py "$RUN" >> "$RUN/m1.log" 2>&1; echo "$(date -Is) M1 done" >> $LOG; }
S20_NW=$H $PY tools/S2.0_edges.py "$RUN" real >> "$RUN/m2_fwd.log" 2>&1 &
P1=$!
S20_REVERSE=1 S20_NW=$H $PY tools/S2.0_edges.py "$RUN" real >> "$RUN/m2_rev.log" 2>&1 &
P2=$!
wait $P1 $P2
S20_NW=$NW $PY tools/S2.0_edges.py "$RUN" real >> "$RUN/m2_fwd.log" 2>&1
echo "$(date -Is) M2 done: $(ls $RUN/edges | wc -l) windows" >> $LOG
$PY tools/S2.0_analysis.py "$RUN" edges >> $LOG 2>&1
echo $NW > "$RUN/WORKERS"
$PY tools/S2.0_inner.py "$RUN" "$RUN/params_flat.json" >> "$RUN/inner.log" 2>&1
echo "$(date -Is) inner done: $(ls $RUN/inner | wc -l)" >> $LOG
$PY tools/S2.0_analysis.py "$RUN" inner >> $LOG 2>&1
RES=${S20_RES:-results/S2.0}
WIN=$($PY -c "import json; print(json.load(open('$RES/S2.0_RA.json'))['winner'])")
echo "$(date -Is) R-A winner: $WIN" >> $LOG
for rule in C0 $WIN; do
  [ "$rule" = "None" ] && continue
  for k in A1 A23 A4 A5 A6; do $PY tools/S2.0_lumens.py "$RUN" "$RUN/params_flat.json" $rule $k >> "$RUN/lumens.log" 2>&1 & done
done
wait
echo "$(date -Is) lumens done" >> $LOG
[ "$WIN" != "None" ] && $PY tools/S2.0_disp.py "$RUN" "$RUN/params_flat.json" $WIN >> "$RUN/disp.log" 2>&1
echo "$(date -Is) disp done: $(ls $RUN/disp 2>/dev/null | wc -l)" >> $LOG
echo "$(date -Is) chain done" >> $LOG
