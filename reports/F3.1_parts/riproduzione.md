
```bash
cd ~/2025.geo_spatialtrans && git checkout f3-methods
bash tools/F3.1_md5.sh before
.venv/bin/python tools/F3.1_reference.py /mnt/micron/geo_spatialtrans/F3.1/main/ref 20261113 A1 A2 A3 A4 A5 A6
.venv/bin/python tools/F3.1_export_real.py /mnt/micron/geo_spatialtrans/F3.1/main/real <aree>   # poi tools/F3.1_b2c.py (label, label_dimg, label_shift) e tools/F3.1_rctd.R
.venv/bin/python tools/F3.1_intervals.py /mnt/micron/geo_spatialtrans/F3.1/main /mnt/micron/geo_spatialtrans/F3.1/explore/explore_measures_v3.csv results/F3.1
.venv/bin/python tools/F3.1_profiles_nuc.py /mnt/micron/geo_spatialtrans/F3.1/profiles_nuc A1 A2 A4 A5 A6
bash tools/F3.1_run_B.sh 6 && bash tools/F3.1_run_PH3.sh 5
.venv/bin/python tools/F3.1_analysis.py /mnt/micron/geo_spatialtrans/F3.1/main results/F3.1
.venv/bin/python tools/F3.1_PH1_mixing.py /mnt/micron/geo_spatialtrans/F3.1/main results/F3.1/F3.1_PH1_mixing.csv
.venv/bin/python tools/F3.1_PH3_analysis.py /mnt/micron/geo_spatialtrans/F3.1/main results/F3.1
.venv-b2c/bin/python tools/F3.1_test.py /mnt/micron/geo_spatialtrans/F3.1/main      # G-1
.venv-b2c/bin/python tools/F3.1_check.py /mnt/micron/geo_spatialtrans/F3.1/main     # G-2
bash tools/F3.1_repro.sh                                                             # G-3
bash tools/F3.1_md5.sh after                                                         # G-4
.venv/bin/python tools/F3.1_figures.py results/F3.1 reports/F3.1_parts && .venv/bin/python tools/F3.1_report_parts.py && quarto render reports/F3.1.qmd
```

Semi: pannello base 20261114 (k = 0, 1; tabella Mersenne-Twister sulle 30 aree ordinate), esplorazione 20261112, riferimento 20261113, intervalli del vero 20261115, intervalli del simulato 20261117, PH-3 20261119, RCTD 20261118, prove 20261116. Uscite in `/mnt/micron/geo_spatialtrans/F3.1/` (25G, non nel repository).

[^bin2cell]: Polański et al. 2024, Bin2cell reconstructs cells from high resolution Visium HD data
[^rctd]: Cable et al. 2022, Robust decomposition of cell type mixtures in spatial transcriptomics (RCTD)
