#!/usr/bin/env bash
# F3.1: RCTD con il seme dichiarato, due volte per area da cartelle separate; i JSON devono coincidere (senza il tempo)
R=/mnt/micron/geo_spatialtrans/F3.1/main; T=/mnt/micron/geo_spatialtrans/F3.1/determinism; rm -rf $T; mkdir -p $T
for a in real/A6_r2 sim/A1_r1_k0; do for i in 1 2; do d=$T/$(echo $a | tr / _)_$i; mkdir -p $d
  for f in blocks.csv blocks_counts.mtx genes.tsv; do cp $R/$a/$f $d/; done; done; done
ls -d $T/* | xargs -P 4 -I{} bash -c 'b=$(basename {}); s=sim; case $b in real_*) s=real;; esac; A=$(echo $b | sed "s/^[a-z]*_//" | cut -c1-2); Rscript --vanilla tools/F3.1_rctd.R {} '$R'/ref $A $s > {}/rctd.log 2>&1 || echo FAIL {}'
python3 - <<'EOF'
import json, glob
T = "/mnt/micron/geo_spatialtrans/F3.1/determinism"
for a in ["real_A6_r2", "sim_A1_r1_k0"]:
    x, y = [json.load(open(f"{T}/{a}_{i}/rctd.json")) for i in (1, 2)]
    for z in (x, y): z.pop("sec", None)
    print(a, "identical" if x == y else "DIFFERENT", {k: (x.get(k), y.get(k)) for k in x if x.get(k) != y.get(k)})
    print("  ", {k: x[k] for k in ["frac_singlet", "agree_nucleus", "neighbour_same"]})
EOF
