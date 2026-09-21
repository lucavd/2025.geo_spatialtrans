# R2 — Pre-registrazione (scritta PRIMA di ogni run; approvata da Luca il 2026-09-20: «procedi»)

Obiettivo unico: segmentazione nucleare su H&E per archetipo A1–A6 su regioni annotate (5 ROI 1×1 mm per archetipo),
tre segmentazioni (Cellpose primario; StarDist 2D_versatile_he; Space Ranger nucleus_segmentations per A2/A4/A5/A6),
prima versione con fonte di docs/BIO_REFERENCES.md (geometria). Nessun codice del cell layer, nessun Voronoi.

Fatto emerso in apertura: Space Ranger 4.x usa un'implementazione custom di StarDist per la segmentazione nucleare
(10x docs, Algorithms overview: Image processing). StarDist-nostro vs SR non e' un confronto indipendente; Cellpose lo e'.

## Check C
| id | Attesa | Soglia | Fonte |
|---|---|---|---|
| C-R2.1 | torch vede la GPU; Cellpose e StarDist caricano e segmentano un ritaglio di prova | PASS/FAIL | ambiente |
| C-R2.2 | Riproducibilita': due run stesso seed stesso ROI -> stesso n oggetti e stessa somma pixel | identita' | — |
| C-R2.3 | Scala: area um2 = area px * (um/px)^2 con um/px da scalefactors_json (R1) | esatta | R1 |
| C-R2.4 | Ogni ROI interamente dentro immagine e maschera in_tissue; nessun ROI su bin fuori immagine (BL-021) | 100 % | R1 |
| C-R2.5 | Orientamento bin<->pixel gestito (A6 ruotato 180°, BL-020): poligoni SR sovrapposti ai nuclei H&E, IoU>0 su >=95 % | >= 95 % | R1 |
| C-R2.6 | Cuore: maschera bolle di montaggio applicata, area mascherata riportata (BL-017) | riportata | R1 |
| C-R2.7 | Nuclei con area < 4 um2 o > 400 um2 flaggati; frazione flag < 5 % | < 5 % | convenzione, da confermare |

## Check B
| id | Attesa | Soglia | Fonte |
|---|---|---|---|
| B-R2.1 | Ordinamento densita' nucleare: A4 > A1 >= A2 > A5 > {A3, A6} | ordine dei mediani | roadmap §3 (ipotesi qualitativa) |
| B-R2.2 | A4: 1 900–6 300 nuclei/mm2 (TLS su H&E: min 0.0019, media 0.0040 ± 0.0010, max 0.0063 /um2) | mediana nell'intervallo | DOI 10.1371/journal.pone.0256907 (surrogato: TLS polmonari) |
| B-R2.3 | A4: area nucleare linfociti mediana 20–45 um2 (letteratura: 35–40 um2 linfociti reattivi) | mediana nell'intervallo | PMC2983036 |
| B-R2.4 | A5: 1 900–3 500 nuclei/mm2 da densita' 3D (neuroni 6.68e4–12.3e4/mm3, glia 3.6e4, endotelio 7.0e4) con Abercrombie N_A = N_V (t + D), t = 5 um, D = 6–10 um; aree bimodali (dip test p < 0.05) | mediana nell'intervallo; p < 0.05 | DOI 10.3389/fnana.2018.00083; DOI 10.1002/ar.1090940210; t: protocollo 10x FFPE (da citare) |
| B-R2.5 | A6: nuclei miociti 28 000 ± 7 200/mm3 -> 250–500/mm2; totale 600–2 000/mm2; tra i piu' bassi con A3 | mediana nell'intervallo | DOI 10.1016/j.cell.2015.05.026 |
| B-R2.6 | A6: eccentricita' nucleare mediana > A1 e A4 di almeno 0.1 | Δ mediane > 0.1 | roadmap §3 |
| B-R2.7 | A1, A2, A3: nessuna fonte 2D; preset design (5500, 1500, 8000) = ipotesi; se |stima − preset| > 30 % -> BL-007 diventa revisione design | registrare | design v1.1 |

## Controprova
| id | Test | Soglia | Fonte |
|---|---|---|---|
| CP-R2.1 | Cellpose vs StarDist stessi ROI: F1@IoU0.5, Dice unione, Δ densita' | F1 >= 0.70 PASS / 0.50–0.70 WARN / < 0.50 FAIL; |Δ| <= 10 % PASS / <= 25 % WARN | convenzione (benchmark H&E F1@0.5 0.67–0.82, arXiv 2411.00078) |
| CP-R2.2 | Cellpose vs SR e StarDist vs SR (A2, A4, A5, A6): stessi indici; attesa StarDist↔SR > Cellpose↔SR | come sopra | 10x docs |
| CP-R2.3 | BL-015: nativo vs 2x vs 4x. Ipotesi: 2x Δ densita' < 10 %; 4x cala > 10 % (A4 piu' sensibile) | come dichiarato | ipotesi |
| CP-R2.4 | Modello nullo: g(r) reale vs CSR; g(r) < 1 per r < ~5 um in tutti gli archetipi, -> 1 oltre 30 um | envelope 95 % | spatstat |
| CP-R2.5 | A3 vs A2: se densita' non distinguibili (IC 95 % sovrapposti) il surrogato BL-018 e' insufficiente | IC 95 % disgiunti | BL-018 |

## Scelte operative approvate
1. 5 ROI 1x1 mm per archetipo scelti da Claude (guida: cluster SR + ritagli R1), miniature consegnate a Luca PRIMA della segmentazione; A3 = stroma nel CRC; A4 = 2 follicoli + 3 paracorticale.
2. Cellpose (Cellpose-SAM; fallback 'nuclei' su ematossilina) + StarDist 2D_versatile_he, GPU, .venv del progetto.
3. Metriche in Python (results/R2/), g(r) e test spaziali in R (spatstat), report reports/R2.qmd.
