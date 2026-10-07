# S1.1 — revisione avversariale indipendente (2026-10-07)
Oggetto: extract_regions() (R/04b1_extract_regions.R) al commit e7f2562, branch step1-cell-layer. Il repo non è stato modificato.
Tabella dei rilievi: S1.1_adversarial_review.csv (26 rilievi: 13 confermati, 2 discrepanze, 11 rilievi di metodo).
Script, in ordine:
- 01_export.R: legge gli RDS di input ed esporta parquet/json (nessuna funzione del repo)
- 02_components.py: componenti 4-connesse indipendenti (scipy.ndimage.label) → py_out/*_components.parquet
- 03_sources.py: ricostruzione delle mappe reali dalle sorgenti, md5, composizione del nullo, stride di I6 → py_out/sources_check.json
- 04_export_outputs.R: esporta r0_df dagli output salvati dal test (per il confronto dei multinsiemi)
- 05_polygons.R: unico uso di extract_regions per le aree, confrontate con le aree scipy → 05_polygons_summary.csv, 05_sample_*.csv
- 06_moore.py: Moore e marching squares sui centri, implementati da zero → py_out/moore_*.csv
- 07_perf.sh: C11 rilanciato 3 volte → 07_perf.csv; 07b_perf_worst.R: mappe frammentate → 07b_perf_worst.csv
- 08_compare.py: confronti con i CSV riportati → py_out/compare.json
- 09_adversarial.R: 16 casi avversari nuovi → 09_adversarial_results.csv

