#!/usr/bin/env bash
# tools/S2.0_thermal_guard.sh — S2.0: sospende fino a 8 processi di calcolo dello step (SIGSTOP) se Tctl > 88 °C, li riprende (SIGCONT) sotto 80 °C.
# Come tools/R5_thermal_guard.sh; processi sorvegliati: tools/S2.0_(edges|inner|disp|lumens|stripes|test).py. Registro: <dir>/guard.log.
D=${1:-/mnt/micron/geo_spatialtrans/S2.0}; F=$D/stopped_pids.txt; L=$D/guard.log; : > "$F"
while true; do
  T=$(sensors | awk '/Tctl/ {gsub(/[+°C]/,"",$2); print int($2)}')
  if [ -s "$F" ] && [ "$T" -lt 80 ]; then kill -CONT $(cat "$F") 2>/dev/null; echo "$(date +%m-%d_%H:%M) T=$T resume $(cat $F)" >> "$L"; : > "$F"
  elif [ ! -s "$F" ] && [ "$T" -gt 88 ]; then
    P=$(pgrep -f 'tools/S2\.0_(edges|inner|disp|lumens|stripes|test)\.py' | sort -n | tail -n +2 | head -8 | tr '\n' ' ')
    if [ -n "$P" ]; then kill -STOP $P; echo "$P" > "$F"; echo "$(date +%m-%d_%H:%M) T=$T stop $P" >> "$L"; fi
  fi
  sleep 60
done
