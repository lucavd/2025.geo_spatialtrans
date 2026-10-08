# Revisione avversariale S1.2 — track "code"

Revisore indipendente (sotto-agente), 2026-10-08. Nessun file del repo e' stato modificato; nessun commit.
Tutto cio' che segue e' dentro `results/S1.2/review/code/`.

## Comandi per rifare il ricalcolo (dalla root del repo, su lesexp-server)
```bash
cd ~/2025.geo_spatialtrans
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent
RV=results/S1.2/review/code
# 1) test e C-perf in copia isolata (non sovrascrive results/S1.2/*.csv)
mkdir -p $RV/rerun/results/S1.2 $RV/rerun/renv && cd $RV/rerun && ln -sfn ~/2025.geo_spatialtrans/renv/library renv/library \
  && cp -r ~/2025.geo_spatialtrans/R ~/2025.geo_spatialtrans/tools . \
  && Rscript --vanilla R/testing/test_S1.2.R > test_rerun.log 2>&1 \
  && for i in 1 2 3; do /usr/bin/time -v Rscript --vanilla tools/S1.2_perf.R > perf_run$i.log 2>&1; done; cd ~/2025.geo_spatialtrans
diff results/S1.2/S1.2_test_results.csv $RV/rerun/results/S1.2/S1.2_test_results.csv
diff results/S1.2/S1.2_mutants.csv      $RV/rerun/results/S1.2/S1.2_mutants.csv
# 2) input avversari (seed_centroids usata solo per generare gli output)
Rscript --vanilla $RV/rc_adversarial.R main                              # ~26 min (A12b da sola ~20 min)
for L in C en_GB.UTF-8; do LC_ALL=$L Rscript --vanilla $RV/rc_adversarial.R locale; done
( ulimit -v 20000000; /usr/bin/time -v Rscript --vanilla $RV/rc_adversarial.R mem )
for M in perf_ref perf_bigd; do /usr/bin/time -v Rscript --vanilla $RV/rc_adversarial.R $M; done
Rscript --vanilla $RV/rc_adversarial.R a11big                            # ~1 min, 20 core
# 3) verifiche indipendenti (numpy/scipy: cKDTree, ray casting proprio, densita' armonica, Madow, C1b)
.venv/bin/python $RV/rc_verify.py        # -> verify/V_cases.csv, verify/V_special.json
.venv/bin/python $RV/rc_verify_a11.py    # -> verify/V_A11_boundary_200seeds.csv
```
Uso dei core: mclapply a 20 core (test del progetto: 24), mai in parallelo fra loro.

## File
- `S1.2_adversarial_review_code.csv` — rilievi RA-code-01..23
- `rc_adversarial.R`, `rc_verify.py`, `rc_verify_a11.py` — codice del revisore
- `out/<caso>/` — input (in_*.csv, poligoni), output (centroids/region_df/type_df), meta.csv; `out/A09_rng.csv`, `out/A13b_area_consistency.csv`, `out/A16*`, `out/A17*`
- `verify/` — esiti delle verifiche; `log_*.txt`, `verify_*stdout*.txt` — log
- `rerun/` — copia isolata del codice al commit 6316790 (md5 in `rerun/provenance_md5.txt`) e log del rilancio

Nota (chiusura S1.2): out/ e rerun/ spostati in /mnt/micron/geo_spatialtrans/S1.2/review/code/ (fuori git, 138 MB).
