#!/usr/bin/env bash
# tools/S1.3d_perf.sh — S1.3d: tempo e memoria su S, ogni (candidato, operazione) in un processo separato, 3 ripetizioni, mediana in analisi.
# Ripetizioni in ordine (prima tutte le rep 1, che salvano il confronto con la verità); P processi in parallelo (dichiarato in dev_notes).
# Uso (root del repo): bash tools/S1.3d_perf.sh [outdir] [P]
OUTD=${1:-/mnt/micron/geo_spatialtrans/S1.3d}; P=${2:-4}; mkdir -p "$OUTD/logs"
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent
for rep in 1 2 3; do
  for c in V0 F1 F2 F3 F4 F5; do echo "$c O1 $rep"; done
  for op in O2s2 O2s8; do for c in V0 F1 F2 F3 F4 F5 F6; do echo "$c $op $rep"; done; done
done | xargs -P "$P" -L 1 bash -c 'log="'"$OUTD"'/logs/perf_$0_$1_$2.log"; [ -s "$log" ] && grep -q "^PERF," "$log" && exit 0; /usr/bin/time -v Rscript --vanilla tools/S1.3d_perf.R $0 $1 $2 '"$OUTD"' > "$log" 2>&1 || echo "perf FAIL $0 $1 $2"'
grep -h "^PERF," "$OUTD"/logs/perf_*.log | sort | wc -l
