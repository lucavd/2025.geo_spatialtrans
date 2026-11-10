#!/usr/bin/env bash
# tools/S2.6_thermal_guard.sh - S2.6: regola il file WORKERS dello scheduler con la temperatura (Tctl > 88 C: -2, min 4; < 80 C: +1, fino a MAX).
# Non sospende i processi (i tempi misurati restano tempi di calcolo). Registro: <dir>/guard.log (temperatura, carico, WORKERS ogni minuto).
D=${1:-/mnt/micron/geo_spatialtrans/S2.6}; MAX=${2:-16}; L=$D/guard.log
while true; do
  T=$(sensors | awk '/Tctl/ {gsub(/[+°C]/,"",$2); print int($2)}'); W=$(cat $D/WORKERS 2>/dev/null || echo $MAX)
  if [ "$T" -gt 88 ] && [ "$W" -gt 4 ]; then W=$((W - 2)); elif [ "$T" -lt 80 ] && [ "$W" -lt "$MAX" ]; then W=$((W + 1)); fi
  echo "$W" > $D/WORKERS.tmp && mv $D/WORKERS.tmp $D/WORKERS
  echo "$(date +%m-%d_%H:%M) T=$T load=$(cut -d' ' -f1 /proc/loadavg) W=$W" >> $L
  sleep 60
done
