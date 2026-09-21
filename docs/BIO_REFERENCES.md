# BIO_REFERENCES — valori di riferimento biologici

Regola: una riga entra qui solo con **fonte tracciabile**: DOI di un articolo, oppure `script` + `dataset` (in `data/real/README.md`) che ha prodotto la stima. Niente valori a memoria.

| id | Archetipo | Metrica | Valore | Intervallo / distribuzione | Unità | Fonte (DOI o script@commit + dataset) | Sessione | Note |
|---|---|---|---|---|---|---|---|---|
| B-001 | A1 | densità nucleare (consenso) | 10876 | min–max 5 ROI: 10788–11473; Cellpose 10586, StarDist 5483 | nuclei/mm² tessuto | tools/R2_consensus.py@63c4ccf + Visium_HD_Mouse_Small_Intestine | R2 | **provvisorio** (BL-024): consenso Cellpose∩StarDist + esclusivi con evidenza ematossilina; 5 ROI 1 mm²; n Cellpose = 38687 |
| B-002 | A1 | area nucleare mediana (Cellpose-RGB) | 12.9 | min–max mediane ROI: 12.6–13.7; StarDist 14.8 | µm² (mediana delle mediane per ROI) | tools/R2_metrics.py@63c4ccf + Visium_HD_Mouse_Small_Intestine | R2 | profili di sezione (FFPE ~5 µm); StarDist +10–15 %; n = 38687; A4 sotto la letteratura (BL-025) |
| B-003 | A1 | diametro equivalente mediano | 4.1 |  | µm | tools/R2_metrics.py@63c4ccf | R2 | da area Cellpose-RGB |
| B-004 | A1 | eccentricità nucleare mediana | 0.80 | StarDist 0.71 | — | tools/R2_metrics.py@63c4ccf | R2 | regionprops; A6 e A1 allungati, A4/A5 rotondi |
| B-005 | A1 | distanza al primo vicino mediana | 4.7 | StarDist 5.7 | µm | tools/R2_metrics.py@63c4ccf | R2 | centroidi Cellpose-RGB nel tessuto valido |
| B-006 | A1 | raggio di inibizione g(r) | 3.5 | g(5 µm) 1.84; g(10) 1.46; g(30) 1.21; Clark–Evans R 0.95 | µm | tools/R2_spatial.R@63c4ccf | R2 | primo r in cui g rientra nell'envelope CSR 95 %; A1 e A6 senza inibizione rilevabile (§4.3 report) |
| B-007 | A1 | frazione di area nucleare | 0.15 | StarDist 0.09 | — | tools/R2_metrics.py@63c4ccf | R2 | somma aree nucleari / area tessuto valido |
| B-008 | A2 | densità nucleare (consenso) | 6739 | min–max 5 ROI: 5511–8231; Cellpose 7901, StarDist 4256 | nuclei/mm² tessuto | tools/R2_consensus.py@63c4ccf + Visium_HD_Human_Colon_Cancer (ROI tumorali) | R2 | **provvisorio** (BL-024): consenso Cellpose∩StarDist + esclusivi con evidenza ematossilina; 5 ROI 1 mm²; n Cellpose = 32721 |
| B-009 | A2 | area nucleare mediana (Cellpose-RGB) | 31.1 | min–max mediane ROI: 28.3–32.3; StarDist 31.7 | µm² (mediana delle mediane per ROI) | tools/R2_metrics.py@63c4ccf + Visium_HD_Human_Colon_Cancer (ROI tumorali) | R2 | profili di sezione (FFPE ~5 µm); StarDist +10–15 %; n = 32721; A4 sotto la letteratura (BL-025) |
| B-010 | A2 | diametro equivalente mediano | 6.3 |  | µm | tools/R2_metrics.py@63c4ccf | R2 | da area Cellpose-RGB |
| B-011 | A2 | eccentricità nucleare mediana | 0.79 | StarDist 0.74 | — | tools/R2_metrics.py@63c4ccf | R2 | regionprops; A6 e A1 allungati, A4/A5 rotondi |
| B-012 | A2 | distanza al primo vicino mediana | 6.3 | StarDist 7.5 | µm | tools/R2_metrics.py@63c4ccf | R2 | centroidi Cellpose-RGB nel tessuto valido |
| B-013 | A2 | raggio di inibizione g(r) | 5.0 | g(5 µm) 1.13; g(10) 1.37; g(30) 1.21; Clark–Evans R 1.12 | µm | tools/R2_spatial.R@63c4ccf | R2 | primo r in cui g rientra nell'envelope CSR 95 %; A1 e A6 senza inibizione rilevabile (§4.3 report) |
| B-014 | A2 | frazione di area nucleare | 0.29 | StarDist 0.15 | — | tools/R2_metrics.py@63c4ccf | R2 | somma aree nucleari / area tessuto valido |
| B-015 | A3 | densità nucleare (consenso) | 4117 | min–max 5 ROI: 3206–4558; Cellpose 4300, StarDist 2639 | nuclei/mm² tessuto | tools/R2_consensus.py@63c4ccf + Visium_HD_Human_Colon_Cancer (ROI stromali) | R2 | **provvisorio** (BL-024): consenso Cellpose∩StarDist + esclusivi con evidenza ematossilina; 5 ROI 1 mm²; n Cellpose = 21298 |
| B-016 | A3 | area nucleare mediana (Cellpose-RGB) | 31.3 | min–max mediane ROI: 27.7–33.6; StarDist 31.9 | µm² (mediana delle mediane per ROI) | tools/R2_metrics.py@63c4ccf + Visium_HD_Human_Colon_Cancer (ROI stromali) | R2 | profili di sezione (FFPE ~5 µm); StarDist +10–15 %; n = 21298; A4 sotto la letteratura (BL-025) |
| B-017 | A3 | diametro equivalente mediano | 6.3 |  | µm | tools/R2_metrics.py@63c4ccf | R2 | da area Cellpose-RGB |
| B-018 | A3 | eccentricità nucleare mediana | 0.81 | StarDist 0.76 | — | tools/R2_metrics.py@63c4ccf | R2 | regionprops; A6 e A1 allungati, A4/A5 rotondi |
| B-019 | A3 | distanza al primo vicino mediana | 7.8 | StarDist 9.0 | µm | tools/R2_metrics.py@63c4ccf | R2 | centroidi Cellpose-RGB nel tessuto valido |
| B-020 | A3 | raggio di inibizione g(r) | 5.0 | g(5 µm) 0.93; g(10) 1.53; g(30) 1.29; Clark–Evans R 1.01 | µm | tools/R2_spatial.R@63c4ccf | R2 | primo r in cui g rientra nell'envelope CSR 95 %; A1 e A6 senza inibizione rilevabile (§4.3 report) |
| B-021 | A3 | frazione di area nucleare | 0.16 | StarDist 0.09 | — | tools/R2_metrics.py@63c4ccf | R2 | somma aree nucleari / area tessuto valido |
| B-022 | A4 | densità nucleare (consenso) | 24491 | min–max 5 ROI: 19781–25574; Cellpose 24367, StarDist 21544 | nuclei/mm² tessuto | tools/R2_consensus.py@63c4ccf + Visium_HD_Human_Lymph_Node_FFPE | R2 | **provvisorio** (BL-024): consenso Cellpose∩StarDist + esclusivi con evidenza ematossilina; 5 ROI 1 mm²; n Cellpose = 114773 |
| B-023 | A4 | area nucleare mediana (Cellpose-RGB) | 13.9 | min–max mediane ROI: 13.5–14.6; StarDist 16.0 | µm² (mediana delle mediane per ROI) | tools/R2_metrics.py@63c4ccf + Visium_HD_Human_Lymph_Node_FFPE | R2 | profili di sezione (FFPE ~5 µm); StarDist +10–15 %; n = 114773; A4 sotto la letteratura (BL-025) |
| B-024 | A4 | diametro equivalente mediano | 4.2 |  | µm | tools/R2_metrics.py@63c4ccf | R2 | da area Cellpose-RGB |
| B-025 | A4 | eccentricità nucleare mediana | 0.61 | StarDist 0.60 | — | tools/R2_metrics.py@63c4ccf | R2 | regionprops; A6 e A1 allungati, A4/A5 rotondi |
| B-026 | A4 | distanza al primo vicino mediana | 4.9 | StarDist 5.0 | µm | tools/R2_metrics.py@63c4ccf | R2 | centroidi Cellpose-RGB nel tessuto valido |
| B-027 | A4 | raggio di inibizione g(r) | 5.0 | g(5 µm) 1.30; g(10) 1.07; g(30) 1.04; Clark–Evans R 1.50 | µm | tools/R2_spatial.R@63c4ccf | R2 | primo r in cui g rientra nell'envelope CSR 95 %; A1 e A6 senza inibizione rilevabile (§4.3 report) |
| B-028 | A4 | frazione di area nucleare | 0.34 | StarDist 0.36 | — | tools/R2_metrics.py@63c4ccf | R2 | somma aree nucleari / area tessuto valido |
| B-029 | A5 | densità nucleare (consenso) | 1127 | min–max 5 ROI: 895–2337; Cellpose 925, StarDist 1172 | nuclei/mm² tessuto | tools/R2_consensus.py@63c4ccf + Visium_HD_6p5mm_Mouse_Brain (isocorteccia) | R2 | **provvisorio** (BL-024); esclusivi Cellpose in A5 per lo piu' senza evidenza nucleare (BL-032): consenso Cellpose∩StarDist + esclusivi con evidenza ematossilina; 5 ROI 1 mm²; n Cellpose = 4142 |
| B-030 | A5 | area nucleare mediana (Cellpose-RGB) | 100.5 | min–max mediane ROI: 74.3–109.8; StarDist 54.7 | µm² (mediana delle mediane per ROI) | tools/R2_metrics.py@63c4ccf + Visium_HD_6p5mm_Mouse_Brain (isocorteccia) | R2 | profili di sezione (FFPE ~5 µm); StarDist +10–15 %; n = 4142; A4 sotto la letteratura (BL-025) |
| B-031 | A5 | diametro equivalente mediano | 11.3 |  | µm | tools/R2_metrics.py@63c4ccf | R2 | da area Cellpose-RGB |
| B-032 | A5 | eccentricità nucleare mediana | 0.55 | StarDist 0.62 | — | tools/R2_metrics.py@63c4ccf | R2 | regionprops; A6 e A1 allungati, A4/A5 rotondi |
| B-033 | A5 | distanza al primo vicino mediana | 18.0 | StarDist 16.4 | µm | tools/R2_metrics.py@63c4ccf | R2 | centroidi Cellpose-RGB nel tessuto valido |
| B-034 | A5 | raggio di inibizione g(r) | 11.5 | g(5 µm) 0.13; g(10) 0.66; g(30) 1.15; Clark–Evans R 1.06 | µm | tools/R2_spatial.R@63c4ccf | R2 | primo r in cui g rientra nell'envelope CSR 95 %; A1 e A6 senza inibizione rilevabile (§4.3 report) |
| B-035 | A5 | frazione di area nucleare | 0.09 | StarDist 0.08 | — | tools/R2_metrics.py@63c4ccf | R2 | somma aree nucleari / area tessuto valido |
| B-036 | A6 | densità nucleare (consenso) | 1044 | min–max 5 ROI: 922–1293; Cellpose 630, StarDist 974 | nuclei/mm² tessuto | tools/R2_consensus.py@63c4ccf + Visium_HD_6p5mm_Human_Heart | R2 | **provvisorio** (BL-024); Cellpose non affidabile in A6 (BL-033): usare StarDist/SR/consenso: consenso Cellpose∩StarDist + esclusivi con evidenza ematossilina; 5 ROI 1 mm²; n Cellpose = 2691 |
| B-037 | A6 | area nucleare mediana (Cellpose-RGB) | 18.2 | min–max mediane ROI: 14.2–22.4; StarDist 23.2 | µm² (mediana delle mediane per ROI) | tools/R2_metrics.py@63c4ccf + Visium_HD_6p5mm_Human_Heart | R2 | profili di sezione (FFPE ~5 µm); StarDist +10–15 %; n = 2691; A4 sotto la letteratura (BL-025) |
| B-038 | A6 | diametro equivalente mediano | 4.8 |  | µm | tools/R2_metrics.py@63c4ccf | R2 | da area Cellpose-RGB |
| B-039 | A6 | eccentricità nucleare mediana | 0.86 | StarDist 0.81 | — | tools/R2_metrics.py@63c4ccf | R2 | regionprops; A6 e A1 allungati, A4/A5 rotondi |
| B-040 | A6 | distanza al primo vicino mediana | 16.2 | StarDist 15.9 | µm | tools/R2_metrics.py@63c4ccf | R2 | centroidi Cellpose-RGB nel tessuto valido |
| B-041 | A6 | raggio di inibizione g(r) | 0.5 | g(5 µm) 1.11; g(10) 1.59; g(30) 1.60; Clark–Evans R 0.77 | µm | tools/R2_spatial.R@63c4ccf | R2 | primo r in cui g rientra nell'envelope CSR 95 %; A1 e A6 senza inibizione rilevabile (§4.3 report) |
| B-042 | A6 | frazione di area nucleare | 0.02 | StarDist 0.02 | — | tools/R2_metrics.py@63c4ccf | R2 | somma aree nucleari / area tessuto valido |

## Stato
- **R2 (2026-09-20/21)**: prima versione della sezione geometria (B-001…B-042): 6 archetipi × 7 metriche, da 30 ROI 1 mm² (5 per archetipo) segmentati con Cellpose-SAM 4.2 (RGB), StarDist 2D_versatile_he e poligoni Space Ranger 4. Tutte le righe di densità sono **provvisorie**: i due segmentatori disaccordano del 12–55 % (direzione tessuto-dipendente) e il valore riportato è una stima di consenso senza conteggio manuale (BL-024). Dettagli, figure e controprove in `reports/R2.html`.
- **S0 (2026-09-18)**: schema creato, **nessun valore inserito**. Le ipotesi numeriche presenti nel design doc (densità del catalogo `cell_types`, `nucleus_to_eq_ratio`, intervallo C8 [5, 25] µm) e in `full_test.R` (sparsità, UMI) **non** sono riferimenti: restano ipotesi finché una riga di questa tabella non le copre.

- **R1 (2026-09-19)**: dataset reali disponibili (vedi `data/real/README.md`), **ancora nessun valore biologico**. Unico dato accertato: scala 0.2737–0.2740 µm/px per le 5 immagini (da `scalefactors_json.json`, verificata con la griglia dei bin, `results/R1/R1_scale.csv`).

## Metriche previste per archetipo (colonne della tabella, da riempire in R2–R5)
Geometria (R2–R3): densità nucleare (cell/mm²), area nucleare (µm²), distanza al primo vicino (µm), g(r) / Ripley K, area ed eq_radius del territorio Voronoi (µm², µm), rapporto N/C, eccentricità del territorio.
Espressione (R4–R5): library size a 2 µm e 8 µm (UMI/bin), frazione di bin vuoti, sparsità, UMI nucleo vs citoplasma, `diffusion_sigma_um`.

## Convenzioni
- `id`: `B-###` progressivo. `Fonte`: DOI **oppure** `script@commit` + nome del dataset come in `data/real/README.md`.
- Una riga stimata da dati reali riporta sempre n (cellule/bin) e l'intervallo (IQR o IC 95 %).
