# R3 — revisione avversariale, revisore «code»

Data: 2026-10-08 · branch `step1-cell-layer` @ 508a74a (working tree con modifiche non committate a `results/R3/R3_selftest_py.csv` e `tools/R3_figures.R`, non toccate).
Scritto solo in `results/R3/review/code/`. Nessun file tracciato modificato, nessun commit.

## Cosa è stato ricalcolato, e come
Implementazione indipendente in Python (`.venv/bin/python`: shapely 2.1.2, scipy 1.18.1, numpy), senza usare `tools/R3_*`,
partendo solo dagli ingressi grezzi: `results/R2/R2_nuclei_all.parquet` (scale==1 & keep), `results/R2/tissue_masks/<A>_<roi>_valid_ds4.png`,
`results/R2/R2_rois_checked.csv` (side_um = side_px·um_per_px), label native `/mnt/micron/geo_spatialtrans/R2/masks/<A>_<roi>_<metodo>_native.npz`.

- Maschera: unione dei rettangoli dei run di pixel validi (convenzione dei lati), `shapely.union_all`.
- Voronoi: `shapely.voronoi_polygons(..., extend_to=ROI+50 µm, ordered=True)`, intersecato col riquadro; ritaglio con `clip_by_rect` + `intersection`.
- Interna: tile contenuta nella maschera, oppure tile entro il ROI con area ritagliata/area tile > 0.999.
- Momenti del territorio: decomposizione in triangoli + autovettori (`numpy.linalg.eigh`), quindi non la formula di Green usata da R3.
- Nuclei: momenti dei pixel nativi (x = col·upp, y = row·upp) + autovettori; contenimento con punto-in-poligono sul territorio (`shapely.contains_xy`), non con il KD-tree di R3.
- CV_loc: kernel gaussiano σ = 5/√(n/area), correzione di bordo uniforme e(u) = (kernel convoluto con la maschera) via FFT, con e senza leave-one-out.
- Δθ: nuclei con ecc_N ≥ 0.8 e territori con e_T ≥ 0.5; 2000 permutazioni (numpy, seed 20261008); controprova a specchio θ_N → −θ_N.

Copertura: tutti i 30 ROI col segmentatore primario, più StarDist nei 10 ROI di A4 e A6 (40 ROI×metodo).
Ci sono anche i test di jitter (A4 p1 SR ±0.5 e ±1 µm; A1 r1 e A2 r4 ±1 µm; 3 semi ciascuno).

## Comandi (dalla root del repo, sul server)
```
.venv/bin/python results/R3/review/code/rc_voronoi.py <A> <roi> <metodo> [jitter_um] [seed]   # un ROI -> roi/*.json, cells/*.parquet
# batch: jobs.txt (generato nel job) | xargs -P 4 -L 1 ... rc_voronoi.py   (max ~8 core)
.venv/bin/python results/R3/review/code/rc_compare.py         # -> rc_roi_all.csv, rc_compare_roi.csv (ROI: ricalcolato vs R3_roi_summary/R3_cp2_roi/cells_all)
.venv/bin/python results/R3/review/code/rc_compare_cells.py   # -> rc_compare_cells.csv (per cellula, appaiato per label, 7 ROI)
.venv/bin/python results/R3/review/code/rc_pool.py            # -> rc_compare_archetype.csv (pool di archetipo vs R3_cells_summary)
.venv/bin/python results/R3/review/code/rc_quant.py           # -> rc_quantization.csv (griglia dei centroidi, orientazioni sugli assi, centroidi sui bordi maschera)
.venv/bin/python results/R3/review/code/rc_mechanism.py       # -> rc_alignment_mechanism.csv (direzione del vicino, Δθ per terzile di gap)
.venv/bin/python results/R3/review/code/rc_six.py             # -> rc_six_roi_table.csv (i 6 ROI richiesti, dichiarato / ricalcolato)
Rscript --vanilla results/R3/review/code/rc_null_seeds.R      # seed dei 1200 nulli
Rscript --vanilla results/R3/review/code/rc_a6r3.R            # diagnosi del nucleo in più in A6 r3
```

## Esito in breve
Tutte le grandezze richieste (1)–(9) sono riprodotte con una seconda implementazione, a precisione di macchina, in 40/40 ROI×metodo.
Le eccezioni sono un nucleo su un bordo di pixel in A6 r3 e i pareggi di contenimento in A4-SR, entrambe senza effetto sui verdetti.
Seed, banda, leave-one-out e ordine tile↔punto sono corretti.

L'allineamento territorio–nucleo in A1–A4 non viene da convenzioni di orientazione (controprova a specchio ≈ 45°) né dalla quantizzazione dei centroidi SR (jitter; StarDist).
È compatibile con il volume escluso di nuclei segmentati non sovrapposti (RA-code-15): l'interpretazione biologica di CP-2b va rivista.
In A4 l'orientazione nucleare SR è quantizzata sugli assi (RA-code-14) e N/C, rapporto dei raggi e CP-3 dipendono dal segmentatore (RA-code-16).
Dettaglio in `findings.csv`.
