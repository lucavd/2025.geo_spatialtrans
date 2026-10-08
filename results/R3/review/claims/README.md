# Revisione avversariale R3 — "claims" (verdetti vs pre-registrazione)

Data: 2026-10-08. Revisore: sotto-agente "claims". Repo `~/2025.geo_spatialtrans`, branch `step1-cell-layer` @ 508a74a
(working tree con modifiche non mie: `results/R3/R3_selftest_py.csv`, `tools/R3_figures.R`). Nessun file tracciato modificato, nessun commit.
Pre-registrazione `results/R3/R3_preregistration.md` identica al commit 3fac4e7 (`git diff --stat 3fac4e7` vuoto).

## Esito in breve
- 55/55 verdetti di `R3_verdicts.csv` riprodotti (valore ed esito) da codice indipendente; regole di decisione conformi al testo.
- Nulli: n CSR = n osservato 30/30, RSA entro ±3, n_failed = 0 in 1200/1200, 20 repliche per ROI x modello.
- Deviazioni non dichiarate di forma (M5 descritto all'inverso, C-R3.3c su pixel e non su poligono, C-R3.4 normalizzata con λ, esclusioni di Luca in CP-4).
- Problemi di interpretazione rilevanti (gravita' alta): d* = 4 µm violato dal 22-28% delle coppie SR in A4 (RA-14);
  orientazioni SR quantizzate a {0°, ±90°} (RA-17); CP-3 dominato dal contatto/bordo delle maschere (RA-18, RA-19);
  selezione delle interne non neutra per S1.3 (RA-21). P4 parte da una premessa falsa in A1-A5 (RA-20).
- C-R3.6 non verificabile: rerun ancora in corso.

Tabella completa: `findings.csv` (24 rilievi).

## Comandi (dalla root del repo)
```
.venv/bin/python results/R3/review/claims/rv_verdicts.py     # 54 verdetti + CP-1/CP-2/CP-4 per ROI; n e n_failed dei nulli
Rscript --vanilla results/R3/review/claims/rv_dip.R          # B-R3.3c (diptest), anche su log e sqrt
.venv/bin/python results/R3/review/claims/rv_interp.py       # d* vs nn reali, P3 con 1/λ, q non direzionale, ordine nematico (submit_job, ~1 min)
.venv/bin/python results/R3/review/claims/rv_cp3_masks.py    # contatto fra maschere ed erosione 1-4 px (submit_job, ~3 min)
.venv/bin/python results/R3/review/claims/rv_interior.py     # frazione interne reale vs CSR vs P4; quantizzazione θ_N SR; eq_r lontano dal bordo
.venv/bin/python results/R3/review/claims/rv_masks_sr.py     # valid_ds4 vs tissue_ds4 (buchi), valori di θ_N SR
.venv/bin/python results/R3/review/claims/rv_sr_lattice.py   # orientazione dei territori SR (non quantizzata)
.venv/bin/python results/R3/review/claims/rv_estimand.py     # B-R3.2 stimatore della previsione vs asserzione; Spearman per ROI
.venv/bin/python results/R3/review/claims/rv_cp4_boot.py     # IC bootstrap di ΔCV (CP-4); conteggio previsioni
```
Uscite: `rv_verdicts_recomputed.csv`, `rv_cp1_roi.csv`, `rv_cp2_roi.csv`, `rv_cp3.csv`, `rv_cp4.csv`, `rv_cp4_bootstrap.csv`,
`rv_B34_spearman_roi.csv`, `rv_B32_estimand.csv`, `rv_nn_dstar_interior.csv`, `rv_cp2_alternatives_roi.csv`, `rv_P2_by_qloc.csv`,
`rv_P2_by_q_phi.csv`, `rv_P2_by_ecc.csv`, `rv_P2_roi_replication.csv`, `rv_cp3_masks.csv`, `rv_interior_fraction.csv`, `rv_eqr_far_from_frame.csv`.

## Indipendenza e limiti
- Le metriche per cellula (aree, e_T, θ_T, frac_out) sono lette da `cells_all.parquet` prodotto dalla pipeline: il ricalcolo e' indipendente
  dall'aggregazione e dalle regole di decisione (R3_analysis.R), non dalla tassellazione. La tassellazione e' coperta da C-R3.1/C-R3.2 (non rifatti).
- Permutazioni CP-2b con RNG diverso (numpy, seed 12345, 2000 permutazioni); ordine nematico con 500 permutazioni.
- Analisi di erosione/contatto (RA-18/19) su 14 casi: A4 5 ROI x {SR, StarDist}, A6 r1 e r3 x {SR, StarDist}, A1 r1, A2 r1 (non su tutti i 40).
- q_loc, ordine nematico, erosione, eq_r lontano dal bordo sono analisi del revisore, post hoc, descrittive.
- C-R3.6 (md5 del secondo run) non verificato: alle 17:44 `/mnt/micron/geo_spatialtrans/R3_rerun` aveva 577 file null (repliche 01-10) e calib vuota.
