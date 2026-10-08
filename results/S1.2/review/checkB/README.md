# Revisione avversariale S1.2, track checkB

Esito: `S1.2_adversarial_review_checkB.csv` (RA-checkB-01..20). Tutti i comandi si lanciano dalla root del repo (`~/2025.geo_spatialtrans`) su lesexp-server.
Lo script 4 usa 34 core e gira in circa 70 minuti.

```
R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.2/review/checkB/rv_export.R            # rds -> CSV (data/)
.venv/bin/python results/S1.2/review/checkB/rv_task1_real_pcf.py                                                     # compito 1
.venv/bin/python results/S1.2/review/checkB/rv_task2_recompute.py                                                    # compito 2
.venv/bin/python results/S1.2/review/checkB/rv_task3_b1b.py                                                          # compito 3
R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.2/review/checkB/rv_task3b_b1b_bias.R    # compito 3, diagnosi seed 42
R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.2/review/checkB/rv_task4_sens.R run 34  # compito 4 (rilanciabile)
R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.2/review/checkB/rv_task4_analysis.R
R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.2/review/checkB/rv_task5_extra.R
```

Nota (chiusura S1.2): data/ e sens/ spostati in /mnt/micron/geo_spatialtrans/S1.2/review/checkB/ (fuori git); gli script li cercano qui: creare symlink o adattare il percorso.
