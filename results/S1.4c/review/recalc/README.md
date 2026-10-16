# S1.4c — revisione avversariale, ruolo RICALCOLO (recalc)

Data: 2026-10-15. Repo a `d710cfb` (branch `step1-cell-layer`), nessun file del repo modificato fuori da questa cartella,
nessun commit. Scratch: `/mnt/micron/geo_spatialtrans/S1.4c_review/recalc/`. Codice proprio, in Python dove possibile
(lingua diversa dal banco); nessuna funzione di `tools/S1.4c_analysis.R` o `tools/S1.4c_decision.R` usata.

## Script (eseguiti dalla root del repo)
| Script | Cosa fa |
|---|---|
| `r01_dump.R` | SOLO I/O: obs (900), yard dei loo S1.4, b4/phi dei gen, righe V-c di S1.4b, celle S1.4 + ev/shape_label del prep, territori SR (A4/A6), finestra del ROI (`r3_roi_window`, input) |
| `r02_verdicts.py` | riferimenti reali (real_shape ∩ ev), metri inter-ROI, verdetti di fold/archetipo, CP-2, punteggi, bootstrap (semi 777 e 31415), decisione, CP-1, B-4 |
| `r03_obs.py` | ri-deriva le distanze per seme B1,B2,B7,B8,S1,S2,S3 da `sim_out` con Voronoi proprio (shapely) + regola degli orfani D-S1.3.2 ricodificata; valida l'operatore sui generatori reali; ricalcola b7_real/b8_real |
| `r04_masks.py` | maschere native → n_bnd, n_c2, S-1, coppie a contatto con ciclo per etichetta e dilatazione binaria 5×5; appaiamento A4/A6 rifatto con shapely; CP-0 con Delaunay scipy (score implementato e score pre-registrato) |
| `r05_raster.py` | rasterizzazione propria (shapely.contains_xy sui centri dei pixel) di 6 file gen_in (A1 r1, A4 f1; fold 1, s01; K-c/K-1/K-2) |
| `r06_c0.R` | C-0: K-c con `derive_cell_geometry()` (R/04b4 invariato), seme seed+100000·s, s=0..10, fold A1:1, A2:2, A4:3, A6:5; metriche ricalcolate (ks.test, Spearman) contro V-c di S1.4b |
| `r07_margins.py` | celle in cui K-2 e K-1 differiscono, verdetti a un fold dal cambio |
| `r08_phi_b4.py` | phi di K-2 (q99 tipo 7) e B-4 per fold |

Uscite copiate qui: `recalc_*.csv/json`, `r0*_*.csv`, `*.log`.

## Esiti (tutti ricalcolati; confronto con `results/S1.4c/*.csv`)
| Quantità | Riportato | Ricalcolato | Accordo |
|---|---|---|---|
| Riferimenti reali per ROI (s1, s2_median, s3_median, b7, b8, n_ev; 30 ROI) | real_refs.csv | max scarto 4.4e-16 | sì |
| Verdetti di fold (6 arch × 3 var × 5 fold × 7 metriche = 630) | fold_verdicts.csv | 0 discordanti; max scarto d 5.6e-16, metro 1.0e-15 | sì |
| Verdetti d'archetipo (126) | verdicts.csv | 0 discordanti (anche n_pass/n_warn/n_fail) | sì |
| CP-2: archetipi che superano | S1 4/6, S2 4/6, S3 6/6 → tutte declassate | idem; n_roi_above identici in 18/18 | sì |
| Pannello | B1,B2,B7,B8 | B1,B2,B7,B8 | sì |
| Punteggi (pannello) | K-c 3, K-1 7, K-2 8 | K-c 3, K-1 7, K-2 8 (su 7 metriche: 3, 10, 12) | sì |
| P(bootstrap) V > K-c | 1.0, 1.0 | 1.0, 1.0 (due semi indipendenti); P(K-2 > K-1) 0.925 (riportato 0.935) | sì |
| Scelta | K-2 | K-2 | sì |
| CP-1 | 6 / 2 = 3 ≤ 3 → smentita | idem (sul pannello: 4/2 = 2 > 1 → non smentita) | sì |
| B-4 (18 celle) | K-2 PASS 6/6 | max scarto 6.1e-16; K-2 PASS 6/6 | sì |
| Distanze per seme ri-derivate da sim_out (900 × 7) | obs/*.rds | max scarto ≤ 9.4e-16 per tutte; n_eval, n_excluded identici | sì |
| Operatore Voronoi proprio sui generatori reali | a_t, d_gen di S1.4 | max rel a_t 8.1e-14, max abs d 6.0e-14 (30 ROI) | sì |
| b7_real, b8_real ricalcolati da S1.4 | prep | max scarto 5.0e-16 | sì |
| S-1 reale dalle maschere native (30 ROI) | real_refs.csv | max scarto 4.4e-16; n_bnd, n_c2 per nucleo identici; coppie di contatto: differenza simmetrica 0 | sì |
| Appaiamento A4/A6 rifatto | prep shape_label | matched e etichette identici 10/10 ROI | sì |
| CP-0 AUC per ROI (30) e riunita (6) | cp0_roi.csv, cp0.csv | max scarto 4.4e-16; n coppie identico | sì |
| Raster indipendente (6 file) | sim_out | n_px per nucleo identico 6/6 (0 nuclei diversi), centroidi identici, n_bnd/n_c2 identici, 0 conflitti | sì |
| C-0 (4 fold × 11 semi = 44 righe) | 330/330 identici | max scarto 1.1e-16 (3/44 bit a bit: KS/Spearman propri) | sì |
| phi (30 fold) | gen | max scarto 5.0e-16 | sì |

## Rilievi (non cambiano verdetti né decisione, salvo dove detto)
- **R-1 scelta K-2 contro K-1 su una sola cella.** K-2 e K-1 differiscono nel pannello solo in A2×B2 (K-2 PASS 4/5, K-1 WARN 3/5)
  e in A6×B7 (K-1 WARN, K-2 FAIL, senza effetto sul punteggio). P(K-2 > K-1) nel bootstrap 0.925 (mio) / 0.935 (riportato);
  P(K-2 ≥ K-1) = 1.0. La regola 4 non chiede margine fra K-1 e K-2: la scelta è conforme, ma poggia su 1 punto.
- **R-2 B-4 di K-2 sul filo in A2 e oltre soglia in 5 fold.** Media dei fold A2 0.0470 (soglia 0.05); fold oltre 0.05: A2 f1 0.0505,
  A2 f4 0.0591, A3 f4 0.0532, A3 f5 0.0804, A5 f1 0.1263. Il PASS 6/6 dipende dalla regola «media sui fold» (dev notes).
- **R-3 CP-1 smentita all'uguaglianza.** 3 ≥ 6/2 = 3; sul pannello effettivo (dopo CP-2) 1 < 4/2 = 2 → non smentita. La
  pre-registrazione fissa le 7 metriche: «smentita» è l'esito corretto, ma 2 dei 3 punti di guadagno in A5+A6 vengono da
  metriche S poi declassate.
- **R-4 score di CP-0.** La pre-registrazione scrive «r_i + r_j − d_ij»; il codice (dichiarato nelle dev notes prima del run) usa
  (r_i + r_j)/d_ij. Le due non sono trasformazioni monotone l'una dell'altra. Con lo score pre-registrato: AUC A1 0.959, A2 0.942,
  A3 0.965, A4 0.975, A5 0.995, A6 0.997 (implementato: 0.960, 0.947, 0.967, 0.974, 0.994, 0.997). Premessa non smentita in entrambi i casi.
- **R-5 sensibilità alla regola degli orfani (D-S1.3.2).** Senza riassegnazione dei frammenti (celle multiparte intere) B1 cambia in
  57/900 osservazioni (max 4.5e-3), 1 verdetto di fold e 1 d'archetipo cambiano, punteggi invariati (3, 7, 8). Con la regola
  ricodificata: 0 scarti.
- **R-6 fragilità locale.** A6×B2 è PASS in tutte e tre le varianti grazie al fold 3 con d/mediana 0.991–0.998: se cambiasse,
  ogni variante perderebbe 1 punto (differenze invariate). A4×B2 K-2 ha il fold 5 a 1.005 (WARN): con un PASS, K-2 = 9.
- **R-7 S-1 e S-2 ridondanti.** La mediana reale di S-2 è 0 in 24/30 ROI (non nulla: A2 r1–r5, A3 r5): la KS di S-2 è dominata dalla massa in 0 e quasi
  coincide con |ΔS-1| (CP-2 A6: d_seg S1 e S2 uguali in 4/5 ROI). Nel punteggio a 7 metriche contano quasi due volte la stessa
  informazione; nessun effetto sulla decisione (entrambe declassate).
- **R-8 robustezza di CP-2.** S1: superano A1 (5/5), A2 (5/5), A4 (3/5), A6 (4/5); S2: A1, A2, A4 (5/5), A6 (4/5); S3: 6/6. La
  declassazione regge anche se l'archetipo più vicino alla soglia (A4 in S1) cambiasse.
- **R-9 contatto simulato molto sopra il reale con K-1/K-2.** Raster proprio, A4 f1 s01: S-1 sim K-2 0.712, K-1 0.706, K-c 0.017
  contro reale 0.373; A1 r1: 0.526 / 0.524 / 0.016 contro 0.431. Le metriche S (declassate) avrebbero comunque penalizzato tutte
  le varianti; K-2 vince sulle metriche B, non sulla forma di contatto.

## Limiti di questa revisione
C-0 ricalcolato su 4 fold su 30 (44 righe su 330); raster su 6 file su 900 (ma tutte le 900 distanze sono ri-derivate da sim_out);
bootstrap con generatore diverso (numpy, semi 777 e 31415), quindi confronto qualitativo; la finestra del ROI è letta con
`r3_roi_window()` (input, non ricalcolata); S-3 (solidità di skimage) non ricalcolata con metodo diverso, solo riletta.
