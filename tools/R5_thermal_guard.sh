#!/usr/bin/env bash
# tools/R5_thermal_guard.sh — R5: sospende 8 processi di calcolo (SIGSTOP) se Tctl > 88 °C, li riprende (SIGCONT) quando scende sotto 80 °C.
# Processi sorvegliati: tools/R5_nc.py, R5_fit.py, R5_ph*.py (post hoc), R5_test.py (G-1). I processi sospesi non perdono lavoro.
# Registro in /mnt/micron/geo_spatialtrans/R5/guard.log. (Seconda versione, 2026-11-02: prima copriva solo nc e fit.)
F=/mnt/micron/geo_spatialtrans/R5/run1/stopped_pids.txt; L=/mnt/micron/geo_spatialtrans/R5/guard.log
while true; do
  T=$(sensors | awk '/Tctl/ {gsub(/[+°C]/,"",$2); print int($2)}')
  if [ -s "$F" ] && [ "$T" -lt 80 ]; then kill -CONT $(cat "$F") 2>/dev/null; echo "$(date +%H:%M) T=$T resume $(cat $F)" >> "$L"; : > "$F"
  elif [ ! -s "$F" ] && [ "$T" -gt 88 ]; then
    P=$(pgrep -f 'tools/R5_(nc|fit|ph[a-z0-9_]*|test)\.py' | sort -n | tail -n +2 | head -8 | tr '\n' ' ')
    if [ -n "$P" ]; then kill -STOP $P; echo "$P" > "$F"; echo "$(date +%H:%M) T=$T stop $P" >> "$L"; fi
  fi
  sleep 60
done
