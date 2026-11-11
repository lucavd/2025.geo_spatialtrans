#!/usr/bin/env bash
# tools/S2.8_mutants.sh - S2.8, cancello G-1: mutanti della batteria test_S2.8.R. Ogni mutante e' una copia di R/ con una modifica;
# e' rilevato se la batteria (gruppi a,b) ha almeno un FAIL. Uscita: results/S2.8/S2.8_mutants.txt
set -uo pipefail
cd "$(dirname "$0")/.."
OUT=results/S2.8/S2.8_mutants.txt; mkdir -p results/S2.8; : > "$OUT"
F=R/05b6_simulate_visium_hd_pieces.R
declare -A SED
SED[M1_fingerprint_without_code]='s/^             code = .geo_code_hash("simulate_visium_hd_pieces"), profiles = .geo_object_hash(profiles),$/             profiles = .geo_object_hash(profiles),/'
SED[M2_fingerprint_without_profiles]='s/code = .geo_code_hash("simulate_visium_hd_pieces"), profiles = .geo_object_hash(profiles),$/code = .geo_code_hash("simulate_visium_hd_pieces"),/'
SED[M4_fingerprint_without_layer]='/^             layer = .geo_object_hash/,/unclass(sf::st_as_binary(R)))))$/d; s/code = .geo_code_hash("simulate_visium_hd_pieces"), profiles = .geo_object_hash(profiles),$/code = .geo_code_hash("simulate_visium_hd_pieces"), profiles = .geo_object_hash(profiles))/'
SED[M3_M_left_raw]='s/list(M = Ms, M_raw = St\$M)/list(M = St$M, M_raw = St$M)/'
det=0; tot=0
for m in M1_fingerprint_without_code M2_fingerprint_without_profiles M3_M_left_raw M4_fingerprint_without_layer; do
  D=$(mktemp -d /tmp/s28_mut_XXXX); cp R/*.R "$D"/; sed -i "${SED[$m]}" "$D/$(basename $F)"
  if cmp -s "$F" "$D/$(basename $F)"; then echo "$m NOT_APPLIED" | tee -a "$OUT"; rm -rf "$D"; tot=$((tot+1)); continue; fi
  log=$(S28_RDIR="$D" R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla R/testing/test_S2.8.R a,b 2>&1); nfail=$(printf '%s\n' "$log" | grep -c '^FAIL' || true)
  tot=$((tot+1)); if [ "${nfail:-0}" -gt 0 ]; then det=$((det+1)); r=DETECTED; else r=MISSED; fi
  echo "$m $r ($nfail FAIL: $(printf '%s\n' "$log" | grep '^FAIL' | cut -c7-60 | tr '\n' ';'))" | tee -a "$OUT"
  rm -rf "$D"
done
echo "mutanti rilevati $det/$tot" | tee -a "$OUT"
