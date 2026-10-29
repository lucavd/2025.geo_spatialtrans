# CHANGELOG — geo_spatialtrans

Una voce per sessione (id sessione nel titolo). Lo stato di ogni sessione è in `docs/ROADMAP.md`; le idee rinviate in `docs/BACKLOG.md`.

## S1.4f, 2026-10-29: nucleo orientato in parte lungo il territorio (branch `step1-cell-layer`), **chiusa**

File dello step `docs/steps/S1.4f.md`, report `reports/S1.4f.html`; decisione tecnica DEC-028 (Cauchy avvolta; regola O: legge della compressione, per spareggio; regola N2: resta il disco). decisione di Luca DEC-029 (si chiude il lavoro sul disegno del nucleo, disco nei preset, eccesso di contatti come limite dichiarato). Letta da Luca il 2026-10-29. BL-171, BL-174, BL-166, BL-142 chiuse; nuova BL-175; aggiornate BL-063, BL-123.

### Codice
- `R/04b4_derive_cell_geometry.R`: orientamento parziale lungo il territorio (colonne opzionali del catalogo `nucleus_orient_conc`, `nucleus_orient_law`; leggi costante, allungamento, compressione; 2δ dalla Cauchy avvolta con la stessa u dell'orientamento a caso). Senza le colonne nuove l'output è identico bit per bit (G-3: K2 ed ET = S1.4e 300/300 e 240/240).
- `tools/S1.4f_bench_lib.R`, `tools/S1.4f_bench.R` (taratura, banchi A e B, G-3, md5), `tools/S1.4f_measure.py` (allineamento per classe di distanza, coppie, P-1, P-2, CP-1, G-4), `tools/S1.4f_calib.py` (taratura lascia-fuori-un-ROI, regola O), `tools/S1.4f_analysis.py` (verdetti, regola N2, controprove, S-1), `tools/S1.4f_refs.py`, `tools/S1.4f_explore*.py`, `tools/S1.4f_posthoc.py` (PH-1, PH-2), `tools/S1.4f_figures.py`, `tools/S1.4f_report_parts.py`, `tools/run_S1.4f.sh`.
- Test `R/testing/test_S1.4f.R` (28 asserzioni) e `tools/S1.4f_test.py` (19), runner `R/testing/run_S1.4f.sh` (18 mutanti, 6 pre-registrati e 12 della revisione).

### Risultati
- Allineamento dei nuclei vicini confermato 5/5; contatti tolti circa la metà dell'ellisse orientata esattamente (tumore del colon 0.639, linfonodo 0.781, ancora smentiti; stroma 0.340 indeterminato); allineamento a 5-8 µm smentito per tutte le leggi (PH-1: i territori sono giusti, il nucleo vero segue il territorio solo vicino al vicino). Regola N2: resta il disco. Controprove superate; sensibilità ai doppioni al bordo senza effetto.

## S1.4e, 2026-10-28: disegno del nucleo e contatti a parità di distanza (branch `step1-cell-layer`), **chiusa**

File dello step `docs/steps/S1.4e.md`, report `reports/S1.4e.html`; decisione tecnica DEC-026 (regola N: resta il disco; coppie di CP-1); decisione di Luca DEC-027 (disco nei preset, orientamento parziale come prossimo tentativo; nessun secondo giudizio di Federica). Letta da Luca il 2026-10-28. BL-166 chiusa in parte; nuove BL-171...BL-174 (BL-172 chiusa per decisione di Luca); aggiornate BL-063, BL-123, BL-129, BL-142.

### Codice
- Nessuna modifica al pacchetto: le ellissi orientate a caso o lungo il territorio esistevano già come opzioni di `derive_cell_geometry()`.
- `tools/S1.4e_bench_lib.R` (bracci K2, ER, ET; ridisegno dei nuclei di un livello di `place_cells()`; invarianti; centri del banco A), `tools/S1.4e_bench.R` (banchi A e B, controlli G-3, md5), `tools/S1.4e_measure.py` (contatto di coppia, P-1, P-2, CP-2, allineamento ed eccentricità, CP-1, G-4), `tools/S1.4e_analysis.py` (verdetti, regola N, controprove), `tools/S1.4e_refs.py`, `tools/S1.4e_explore*.py`, `tools/S1.4e_posthoc_raters.py`, `tools/S1.4e_figures.py`, `tools/S1.4e_report_parts.py`, `tools/run_S1.4e.sh`.
- Test `R/testing/test_S1.4e.R` (20 asserzioni) e `tools/S1.4e_pairs_test.py` (11), runner `R/testing/run_S1.4e.sh` (13 mutanti).

### Risultati
- Ellisse lungo il territorio: nuclei a contatto da -0.080 a -0.098 nei quattro tessuti densi; ellisse a caso invariata; regola N 4 contro 3 (serviva +2): resta il disco. Allineamento dell'ellisse orientata oltre il vero in 5 tessuti su 5.
- Sui punti di Luca, sotto 6 µm: intestino riprodotto (0.33 contro 0.35); tumore del colon e stroma ancora in eccesso (0.72 e 0.74 contro 0.30 e 0.47); nel tumore Luca e Federica divergono sulle stesse coppie (0.17 contro 0.44).

## S1.5e, 2026-10-26: vuoti che tolgono cellule e vuoti allungati nella semina (branch `step1-cell-layer`), **chiusa**

File dello step `docs/steps/S1.5e.md`, report `reports/S1.5e.html`; decisione tecnica DEC-024 (esito delle regole A, B, Q; opzioni nuove); decisione di Luca DEC-025 (i controlli falliti come scritti non bloccano; le quote si tolgono al consolidamento, BL-170). Letta da Luca il 2026-10-27 («vanno bene i tuoi consigli»).

### Codice
- `R/04b2_seed_centroids.R`: colonne facoltative della tabella `heterogeneity`: `void_mode` («move» come S1.5d, «remove»: i vuoti tolgono cellule, numero per regione riferito all'area senza vuoti), `void_aspect` (campo dei vuoti anisotropo, orientamento per regione), `dens_cap` (tetto della densità all'80 % della saturazione, fattore riscalato a media 1 fuori dai vuoti). Senza le colonne la semina è identica a S1.5d (48/48 configurazioni, revisione).
- `R/04b0_tissue_presets.R`: passaggio delle colonne nuove; preset invariati.
- Test `R/testing/test_S1.5e.R` (16 asserzioni), runner `R/testing/run_S1.5e.sh` (13 mutanti come copie del codice).
- Strumenti `tools/S1.5e_*`: esplorazione (densità di riferimento, affollamento, contatti per distanza, forma dei vuoti), calibrazioni (X, XA, XAD, XADs), banco (2 100 unità), finestre, metriche, analisi, controlli G-4, post hoc, figure, parti del report e registro delle verifiche.

### Risultati
- Regola A: restano i vuoti che spostano le cellule (19 contro 15, 15, 16 su 40). Regola B non adottata (58 contro 19, taglie 19 contro 21). Regola Q: senza quote (60 contro 58, 48 contro 46); applicazione: decisione di Luca.
- I contatti in eccesso vengono dal disegno del nucleo (a parità di distanza 0.96 contro 0.47 nel tumore del colon), non dalla semina. Il conteggio manuale misura la densità del tessuto pieno. Togliere cellule porta la densità dell'intestino al limite basso del vero ma la peggiora nel linfonodo e nel cervello, dove la calibrazione usa i vuoti come sostituto di distanze irregolari; i vuoti allungati migliorano il cervello e sono troppo pochi e allungati nell'intestino.
- Cancelli: G-1 PASS (16/16 dopo la revisione), G-2 9/9 (13/13), **G-3 fallito come scritto** (6 cellule non piazzate nella stessa unità di Q di S1.5d), G-4 PASS, G-5 16 942/16 945. CP-1, CP-2 fallite nel cuore (rumore della calibrazione), CP-3 vuota per costruzione (errore di pre-registrazione). Revisione: XAD rifatto come pre-registrato (XADs), nessun taglio degli intervalli (1 verdetto su 1 074 cambia), 4 mutanti nuovi rilevati. Rinvii BL-166...BL-170.

## S1.5d, 2026-10-25/26: spazi vuoti e densità indipendenti dal tipo nella semina, modello di taglia dagli scarti veri (branch `step1-cell-layer`), **chiusa**

File dello step `docs/steps/S1.5d.md`, report `reports/S1.5d.html`; decisione tecnica DEC-022 (esito delle regole A, C, B); decisione di Luca DEC-023 (G-2 e G-3 falliti come scritti non bloccano, tre errori di pre-registrazione, regola C sulla semina uniforme come scritta). Letta da Luca il 2026-10-26 («chiudi lo step»); prossimo step S1.5e (semina eterogenea, parte 3). BIO_REFERENCES B-232...B-264. BACKLOG: chiusa BL-159; aggiornate BL-157, BL-158, BL-160; nuove BL-162...BL-165.

### Codice
- `R/04b2_seed_centroids.R`: `seed_centroids(heterogeneity = ..., quotas = ...)`: campo dei vuoti (soglia τ = Φ⁻¹(1 - φ_v)) e campo della densità exp(σ_d W), indipendenti dal tipo; a_t pesati con l'eterogeneità; griglia di lambda_max raffinata per ℓ_d; campi registrati nelle informazioni della semina.
- `R/04b4_derive_cell_geometry.R`: residui empirici del modello di taglia per gruppo e decimo di log A_T (`resid_table`, `resid_mode = "decile"|"pooled"`), stessa sequenza casuale del modello gaussiano; colonna facoltativa `nucleus_area_resid_shift` (braccio di sensibilità U3c).
- `R/04b0_tissue_presets.R`: tabella `.TP_HET` (soli vuoti, B-259...B-264), compartimenti ricalibrati (B-232...B-248), `.TP_SIZE_RESID_DEFAULT = "gaussian"`, tabelle `inst/extdata/size_resid.csv`.
- `R/04b_cell_placement.R`: `place_cells(override_heterogeneity, quotas, size_resid, resid_mode)`; documentazione: senza `override_heterogeneity` i compartimenti usano i vuoti del preset (RC-11).
- Test `R/testing/test_S1.5d.R` (32 asserzioni), runner `R/testing/run_S1.5d.sh` (20 mutanti come copie del codice).
- Strumenti `tools/S1.5d_*`: esplorazione, riferimenti veri, modello di taglia, calibrazioni, banco (3 600 unità), metriche, analisi, controllo della calibrazione, post hoc, figure, parti del report e registro delle verifiche.

### Risultati
- Con i soli vuoti: spazio vuoto 6/6 tessuti (S1.5c 3/6), territori grandi 4/6 (0/6), variabilità dei territori dell'intestino +0.02 (S1.5c -0.41); contatti sopra il vero in tumore del colon (0.686) e stroma (0.360). Regola A: soli vuoti, 11 contro 10.
- Taglia: residui empirici per decimo di territorio riportano le cellule quasi senza citoplasma a 0.0221 (vuoti 0.0435), ma la regola C è in parità (34 contro 34): resta il gaussiano. La controprova senza decimi (0.0375) mostra che conta la dispersione che cresce col territorio.
- Regola B: non adottata (58 contro 19, taglie 19 contro 21). G-2 9/11 come scritto (20/20 dopo la revisione); G-3 fallito (1 076 cellule non piazzate in bracci non vincenti); G-1 con due errori di tolleranza dichiarati prima. Revisione: 11 393 confronti identici, nessun errore che cambi i numeri.

## S1.5c, 2026-10-23: semina a compartimenti per tipo (branch `step1-cell-layer`), **chiusa**

File dello step `docs/steps/S1.5c.md`, report `reports/S1.5c.html`; decisione tecnica DEC-020 (regole delle tre scelte); decisione di Luca DEC-021 (G-4 fallito come scritto è un errore della tolleranza pre-registrata, non blocca lo step). Letta da Luca il 2026-10-23 («chiudi lo step»).

### Codice
- `R/04b2_seed_centroids.R`: `seed_centroids(spatial = "compartments", compartments, labels = "field"|"blind")`: campi gaussiani latenti per tipo (caratteristiche di Fourier casuali), probabilità locali (softmax), densità locale dalla densità del tipo puro, spostamenti a_t sulla griglia della regione, diradamento, quote per tipo, distanza minima della coppia; contatori diagnostici; ripiego per regioni troppo sottili. Il modo `"uniform"` è invariato (300/300 unità identiche a S1.5b).
- `R/04b0_tissue_presets.R`: tabella `.TP_COMPART` (23 tipi, B-209...B-231), `tissue_preset_spec(spatial, min_dist_mode = "tissue"|"nucleus")`, `.TP_SPATIAL_DEFAULT = "uniform"`.
- `R/04b_cell_placement.R`: `place_cells(spatial, min_dist_mode, labels, override_compartments)`.
- Test `R/testing/test_S1.5c.R` (28 asserzioni), runner `R/testing/run_S1.5c.sh` (12 mutanti come copie del codice, senza interruttori in produzione).
- Strumenti `tools/S1.5c_*`: esplorazione, calibrazione, aree nucleari per tipo, metriche comuni a reale e simulato, riferimenti reali, banco (1 800 unità), contatti per finestra, metriche per unità, analisi, controllo della calibrazione e indagine Monte Carlo, post hoc, figure.

### Risultati
- Conferme spaziali 52 su 92 contro 19 della semina uniforme (raggruppamento dei tipi 23/34 contro 4/34); taglie invariate; controprova con i tipi mescolati superata.
- Territori ancora troppo uniformi in 4 tessuti (coda dei territori grandi accanto agli spazi vuoti, PH-1); contatti del linfonodo invariati; cellule senza citoplasma del tumore del colon 3.4 %: la regola 3 lascia la semina uniforme come predefinita.
- G-4 fallito come scritto (75/79), spiegato dall'errore Monte Carlo del riferimento; decisione a Luca. Revisione: ricalcolo identico, nessun errore che cambi numeri; 4 test aggiunti. Rinvii BL-157...BL-161.

## S1.5b — 2026-10-22/23 — I tipi cellulari nei sei preset (branch `step1-cell-layer`) — **chiusa**

File dello step `docs/steps/S1.5b.md`, report `reports/S1.5b.html`; decisioni DEC-016 (gruppi di taglia, densità del tipo puro; tecnica), DEC-017 (composizione del cervello con i nuclei persi tipizzati dall'espressione; Luca), DEC-018 (stroma misto provvisorio; Luca), DEC-019 (braccio ricentrato, semina cieca al tipo, tipi predefiniti; tecnica). BIO_REFERENCES B-186…B-208. BACKLOG: chiuse BL-055, BL-138, BL-150, BL-099 (parte dei tipi); nuove BL-153…BL-156.

### Codice
- `R/04b0_tissue_presets.R`: tabella dei tipi `.TP_TYPES` (23 tipi, righe B-186…B-208), `tissue_preset_spec(types = "multi", size_arm = "recentred", composition = "tissue"|"follicle"|"paracortex")` predefinito; `types = "single"` = preset di S1.5 invariato.
- `R/04b2_seed_centroids.R`: `type_order = "blind"` (una semina per regione, tipi permutati dopo).
- `R/04b_cell_placement.R`: `place_cells(types, size_arm, composition)`, semina cieca con i preset per tipo, `size_group` ed `epithelium_label` in `cell_df`.
- `tools/S1.5_single_preset.R`: gli script e i test di S1.5 restano a un tipo (md5 identici).
- Test `R/testing/test_S1.5b.R` (22 asserzioni), runner `R/testing/run_S1.5b.sh` (8 mutanti; propaga i fallimenti di S1.4d, regressioni, OKF).
- Strumenti `tools/S1.5b_*.py|R`: calibrazione, riferimenti reali, tipizzazione dei nuclei persi del cervello, attesi, banco (900 unità × 2 corse), analisi, controllo della calibrazione, figure, post hoc.

### Risultati (banco v2, dopo la correzione della semina)
- Taglia del nucleo per tipo confermata in 14 tipi su 17 (nullo senza tipi: 9); cellule quasi senza citoplasma nel tumore del colon 2.2 % (≤ 2.9 % reale); neuroni/glia 3.55; valle fra le due taglie del cervello non riprodotta (BL-154); rapporto nucleo/territorio per tipo smentito in 11 tipi su 17 (BL-153).

## R6 — 2026-10-22/23 — Classificazione per tipo dei nuclei reali dall'espressione a 2 µm (branch `step1-cell-layer`) — **chiusa**

File dello step `docs/steps/R6.md`, report `reports/R6.html`; decisione tecnica DEC-012 (impronta, clustering, marcatori); decisioni di Luca DEC-013 (riferimenti dal criterio corretto), DEC-014 (epitelio del colon per istologia), DEC-015 (linfonodo provvisorio, una taglia per i linfociti), 2026-10-23.

### Codice (nessun file degli step passati modificato)
- `tools/R6_common.py` (impronte bin→nucleo: espansa, solo nucleo, spostata; aggregazione dei conteggi), `tools/R6_aggregate.py` (30 ROI × 3 impronte, appaiamento con le cellule di Space Ranger), `tools/R6_dexp.py` (espansione misurata da Space Ranger), `tools/R6_export_terr.R` (territori reali in parquet), `tools/R6_cluster.py` (Leiden + criterio pre-registrato), `tools/R6_catalogue.py`, `tools/R6_relabel.py` (PH-1), `tools/R6_analysis.py` (riferimenti per tipo, P-1…P-6, CP-1, CP-2, sensibilità; bootstrap a blocchi sulle ROI, B = 2000), `tools/R6_posthoc.py` (PH-2…PH-6), `tools/R6_figures.py`, `tools/R6_report_parts.py`; cancelli `tools/R6_test.py` (G-1), `tools/R6_md5.sh` (G-2), `tools/R6_repro.sh` (G-3), `tools/R6_check.py` (G-4).
- Correzioni dopo la revisione avversariale: filtro Low UMI di Space Ranger (pandas 3), soglia del 5 % di P-2 su tutti i nuclei, bootstrap a due stadi per P-5.

### Risultati
- Cancelli tutti superati; revisione: 316 numeri ricalcolati, 3 discrepanze corrette, di cui due causate da un bug vero in `tools/R6_analysis.py` (filtro «Low UMI» di Space Ranger inefficace in P-6, trovato dal revisore numerico) e una da una soglia applicata diversamente dalla pre-registrazione; il revisore del codice non ha trovato altri bug che cambino i numeri; analisi rieseguita dopo le correzioni.
- Criterio di etichettatura pre-registrato difettoso (cellule di Paneth → «endoteliali», tre cluster di cardiomiociti «non risolti»); correzione PH-1 con accordo con Space Ranger 0.91 (cervello) e 0.79 (cuore); kappa 0.83 e 0.29 contro 0.13 e 0.10 con impronte spostate.
- Neuroni 111 µm² contro glia 34 (P-1); i tipi spiegano la taglia in intestino, tumore e stroma del colon (P-2 confermata), non nel linfonodo (B e T uguali) né nel cuore; il nucleo del cardiomiocita non è più allungato (P-3 smentita); le quote del design smentite in 5 tessuti su 6 (P-5).
- BIO_REFERENCES B-161…B-185 (riferimenti per tipo, provvisori fino alle decisioni); BACKLOG BL-147…BL-152; BL-099, BL-138, BL-064, BL-055 aggiornate.

## S1.5 — 2026-10-21 — Il simulatore del livello cellulare tutto insieme: orchestratore, controlli, sei preset, scala (branch `step1-cell-layer`) — **chiusa**

File dello step `docs/steps/S1.5.md`, report `reports/S1.5.html`; decisioni DEC-008 (perimetro), DEC-009 (distanza minima dai tessuti veri), DEC-010 (nuclei in piu' parti: nessuna modifica), DEC-011 (caso limite «scala sbagliata → C8» fallito come scritto: errore di pre-registrazione).

### Codice (nessun file degli step passati modificato nel comportamento predefinito)
- `R/04b_cell_placement.R`: `place_cells()` (catena completa, tre semi indipendenti per semina, taglia e nuclei pallidi; `stride` dai metadati o errore, BL-042; fusione delle componenti piccole, BL-043; conversione del modello di taglia alla densita' manuale `size_conversion = "preserve_area"`; `keep_tiles`; `allow_estimate`), `merge_small_components()`, `plot_cell_layer()` (4 modi).
- `R/04b0_tissue_presets.R`: `tissue_preset_spec()` per i 6 preset con la tabella delle fonti (B-043...B-054, B-091...B-111, B-113, B-121...B-125, B-139...B-143, B-150...B-154).
- `R/04b5_cell_layer_validation.R`: `validate_cell_layer()` (C1-C9; C2 senza `st_union`, coppie nei due ordini, tolleranze sulla griglia, NA = ERROR; C9 lega righe e geometrie) e `overlay_prec()` (BL-070, BL-079, BL-081, BL-093...BL-097, BL-116).
- `R/04b6_estimate_pixel_size.R`: `estimate_pixel_size()` dai metadati (`scalefactors_json` controllato sulla griglia dei bin) e dall'immagine (laplaciano della gaussiana, rifiuto, controllo di coerenza); la seconda fuori dall'uso raccomandato.
- Test `R/testing/test_S1.5.R` (83 asserzioni), runner `R/testing/run_S1.5.sh` (test e mutanti M1-M9); banco `tools/S1.5_bench.R`, `tools/S1.5_real_refs.R`, `tools/S1.5_real_multipart.py`, `tools/S1.5_faults.R`, `tools/S1.5_cp1_resample.py`, `tools/S1.5_windows.py`, `tools/S1.5_analysis.R`, `tools/S1.5_posthoc.R`, `tools/S1.5_regress_cmp.R`, `tools/S1.5_fig_crops.py`, `tools/S1.5_figures.R`; controprova dopo la revisione `tools/S1.5_review_mutants.sh|py`, `tools/S1.5_review_mutants_analyse.py`, `tools/S1.5_review_faults.R`.
- Correzioni dalla revisione del codice (commit `ab03c27`): FND-01...FND-04, FND-06...FND-10, FND-12...FND-15.
- Design v1.6.

### Risultati
- Cancelli: G-1 81/83 (E10, E10b: caso limite «scala sbagliata → C8» non realizzabile, errore di pre-registrazione, DEC-011); G-2 9/9; G-3 600/600 senza ERROR, rerun e catena manuale 6/6, ricontrollo dopo le correzioni 600/600 identico; G-4 guasti 630/630, nulli 1 800/1 800; G-5 274/277 (3 rilievi di testo e sensibilita').
- Regola della scelta 2: vince la conversione che conserva il legame nucleo/territorio (P-1 6 contro 4).
- Previsioni: P-1 6/6, P-3 6/6, P-7 6/6 confermate; P-2 smentita in intestino e tumore del colon; P-4 smentita nel tumore del colon; P-5 smentita in tumore del colon e linfonodo, indeterminata in stroma e cuore; P-6 smentita in 5 tessuti (cuore confermata, fragile); P-8 smentita nel linfonodo (310 s, 8.50 GB). Post hoc: con i generatori reali contatti e rapporto nucleo/territorio tornano (causa: semina omogenea).
- CP-1 smentita: la stima della scala dall'immagine esce dall'uso raccomandato.
- BIO_REFERENCES B-155...B-160 (densita' di riferimento del modello di taglia). BACKLOG: chiuse BL-016, BL-029, BL-042, BL-043, BL-050, BL-070, BL-076, BL-079, BL-081, BL-093, BL-094, BL-095, BL-097, BL-106, BL-116, BL-139; nuove BL-143...BL-146; rinvii agli step dei tipi (BL-055, BL-064, BL-099, BL-138) e della semina (BL-051, BL-068).

## R3b — 2026-10-20/21 — Maschera nucleare unica per linfonodo e cuore (branch `step1-cell-layer`) — **chiusa**

File dello step `docs/steps/R3b.md`, report `reports/R3b.html`; decisioni DEC-006 (maschera CO, taglia dai nuclei di forma StarDist) e DEC-007 (G-4 fallito come scritto: errore di scrittura).

### Codice (nessun file degli step passati modificato)
- `tools/R3b_consensus.py`: maschera CO (appaiamento reciproco di massima sovrapposizione, StarDist appaiati sopra gli Space Ranger rimasti), tabella dei nuclei R2, filtro cp2, cancello G-4; test `tools/R3b_consensus_test.py` (34/34) e runner `R/testing/run_R3b.sh` (mutanti 5/5, mutanti del revisore 12/12).
- `tools/R3b_bench.R` (generatori R3, territori S1.4, M2 per fold, phi, nucleo ritagliato, controprova con i punti di Luca), `tools/R3b_nuc.py`, `tools/R3b_p1p2.py`, `tools/R3b_windows.py`, `tools/R3b_export.R`, `tools/R3b_analysis.py`, `tools/run_R3b.sh`; post hoc `tools/R3b_posthoc*.py|R`, `tools/R3b_refs.R`, `tools/R3b_fig_mask.py`, `tools/R3b_report_parts.py`.

### Risultati
- Regola di decisione: CO in linfonodo e cuore. Conteggio CO / punti di Luca 1.04 e 1.09 (StarDist 0.77 e 0.84); contatti simulati 0.69 e 0.05; con i punti di Luca come centri 0.72 e 0.015; troncati 1-3 %.
- Post hoc: contatti nella maschera CO del linfonodo 44/59 contro 32 di Luca (64 % via poligoni Space Ranger); poligoni Space Ranger rimasti precisione 0.64 e 0.70; taglia dai soli nuclei di forma StarDist (sigma 0.210 contro 0.288).
- BIO_REFERENCES B-133...B-154 sostituiscono B-073...B-078, B-085...B-090, B-100...B-102, B-106...B-108, B-112, B-114 e la parte q99 di B-124, B-126.
- BACKLOG: chiuse BL-100, BL-115, BL-121, BL-130; nuove BL-141, BL-142; BL-064 rinviata a S1.5.

## S1.4d — 2026-10-20 — Il nucleo ritagliato nel territorio e riportato all'area giusta (K-2) nel pacchetto (branch `step1-cell-layer`) — **chiusa**

Primo step con il protocollo di sessione v2 (DEC-002); file dello step `docs/steps/S1.4d.md`, report `reports/S1.4d.html`.

### Passaggio alla versione 2
- `tools/okf.py` adattato da MMM con i suoi test (`tools/test_okf.py`, 8/8; DEC-001); bundle OKF `docs/decisions/`, `docs/steps/`, `docs/index.md`, `docs/log.md`; decisione sul protocollo v2 registrata (DEC-002).

### Codice
- `R/04b4_derive_cell_geometry.R`: `nucleus_model = c("clip", "scale")`, «clip» predefinito (K-2: ritaglio nel territorio rientrato di 1e-5 µm, ingrandimento fino all'area del modello di taglia, tetto `nucleus_nc_max` = φ per tipo), «scale» identico a `db4afea`; solutore con falsa posizione di Illinois, ricaduta sulla bisezione e arresto legato alla griglia (euristico, BL-140); nuclei in `sf` via WKB (DEC-003).
- Test `R/testing/test_S1.4d.R` (36), runner `R/testing/run_S1.4d.sh`, mutanti `tools/S1.4d_mutants.R`; banco `tools/S1.4d_bench.R`, finestre `tools/S1.4d_windows.py`, analisi `tools/S1.4d_analysis.R`, controllo cella per cella `tools/S1.4d_g4a_cells.R`, post hoc `tools/S1.4d_posthoc_nc.R`, `tools/run_S1.4d.sh`.
- Script di S1.4, S1.4b, S1.4c e delle loro revisioni: 19 chiamate con `nucleus_model = "scale"` (riproducibilità); `test_S1.4.R` in modo «scale».
- Design v1.5: §5.4 riscritta (modello «clip»), operatore di osservazione dei check B (BL-112), frazione di contatti sulle finestre a peso uguale (BL-131), `nucleus_nc_max` in §4.4.

### Risultati
- Cancelli: test e mutanti superati; 0 nuclei fuori e 0 sovrapposti su 6 823 226; equivalenza con il prototipo cella per cella (2 411 650 celle) e del modo «scale» con `db4afea` (330/330); revisione avversariale 38/38 numeri concordanti. Controllo di precisione delle aree fallito come scritto per un nucleo di 0.18 µm² in un nullo: errore della soglia (DEC-004, Luca).
- Frazione di nuclei a contatto contro il conteggio manuale: confermata in 5 tessuti, linfonodo indeterminato (0.744 contro ≤ 0.750); con centri StarDist 0.41 (BL-130 spiegata). Controprova: il modello vecchio non produce contatti nei tessuti densi.
- Tumore del colon: 5.8 % di cellule quasi senza citoplasma contro ≤ 2.9 % reale, dalla coda del modello di taglia; φ resta 1, correzione con i tipi cellulari (DEC-005, Luca; BL-138).
- Tempo sulla sezione S 224.6 s (prototipo 622.7 s), 6.36 GiB.

### Rinvii
- BL-138 (coda del modello di taglia, S1.5), BL-139 (nuclei in più parti, S1.5), BL-140 (arresto del solutore, da assegnare); BL-130 → R3b. Risolte BL-112, BL-118, BL-122, BL-127, BL-131; BL-098 implementata.

## R2c — 2026-10-16/20 — Verità manuale sul contatto fra nuclei (BL-120) (branch `step1-cell-layer`) — **chiusa con riserve**

### Decisioni (ROADMAP §6)
- Disegno (Luca, 2026-10-16): compito sì/no su coppie proposte; 15 nuclei focali per finestra; Federica su metà delle finestre.
- Deviazioni D-1…D-4 approvate da Luca prima delle annotazioni (candidati senza Delaunay; d_max combinato; fascia verificata a
  campione, proposta di Luca; trappole oltre d_geo). D-5, D-6 trovati dai test prima delle risposte.
- **Riferimento S-1 da S1.4d in poi = valore manuale per archetipo con IC bootstrap a due stadi (B-127…B-132)**; segmentatori solo
  descrittivi; S-2/S-3 fuori (Luca, 2026-10-20).

### Codice (`tools/`, nessun file del pacchetto)
- `R2c_candidates.py` (d_max, candidati, fascia, focali, item), `R2c_annotator.py` + `R2c_annotator_core.js` (annotatore sì/no,
  zip con PNG senza perdita), `R2c_check_annotator.py` (C-2), `R2c_refs.py` (segmentatori sui focali, K-2 nelle finestre),
  `R2c_analysis.py` (+ `R2c_analysis_test.py`, 24 test, mutanti 3/3), `R2c_dryrun.py`, `R2c_posthoc.py`, `R2c_rereview_items.py`,
  `R2c_cp6.py`, `R2c_figures.py`, `R2c_fig_eligibility.py`; `R/testing/test_R2c.sh` (9 asserzioni); `reports/R2c.qmd`.

### Risultati (`results/R2c/`; zip degli annotatori in `/mnt/micron/geo_spatialtrans/R2c/`)
- Risposte: Luca 1 607/1 607 (+ 94 del riesame), Federica 834/834.
- Check C: C-1 PASS (dopo D-1/D-2; regola pre-registrata FAIL A1/A2, `prereg_rule/`), C-2 PASS (ritagli pixel-identici, marcatori
  1e-13 µm, cecità, ordine riproducibile), C-3 48/48 identici a cp2, C-4 3/3.
- S-1 manuale (Luca): A1 0.50, A2 0.37, A3 0.18, A4 0.55, A5 0.08, A6 0.02; B-dec «nessuno» in A1–A4, A5 StarDist (per 0.0005),
  A6 indeciso; K-2 PASS A1, A2, A5, A6, FAIL A3 (fragile), WARN A4 (0.77). P-1, P-2 smentite; P-3 non valutabile.
- Controprove: CP-1 κ 0.765, CP-2 0.810 / 1.000, CP-3, CP-4 AUC 0.83–0.95, D-3 0/90, CP-5 3/4 PASS; CP-6 (post hoc, riesame):
  abitudine a «separati» non sostenuta (2/64 contro 1/15, p 0.91).
- Revisione avversariale: 2 revisori; 167/170 quantità identiche, 0 bug attivi, 1 latente (RV-12), verdetti P-1…P-3 mancanti
  (RV-10, aggiunti), scelte non dichiarate (RV-11, dichiarate); sensibilità post hoc (IC onesti, confronto appaiato,
  leave-one-window-out, K-2 a peso uguale, accordo a 3 categorie).

## S1.4c — 2026-10-15/16 — BL-098 riformulata: forma del nucleo, scala contro ritaglio (branch `step1-cell-layer`) — **chiusa con riserve**

### Decisioni (ROADMAP §6)
- Perimetro (Luca): banco + decisione; implementazione in `04b4` in S1.4d. Varianti: K-c, K-1, K-2 nel Voronoi.
- **BL-098 → K-2** (ritaglio nel territorio rientrato di 1e-5 µm + gonfiaggio a pari area, tetto φ = min(1, q99 N/C)), per regola
  pre-registrata; Luca: «Nessun veto: K-2 per S1.4d, fragilità tra le riserve».

### Codice (prototipo, non pacchetto)
- `tools/S1.4c_shape.R` (K-1, K-2, `obs_territories()`, B-8), `tools/S1.4c_raster.py` (operatore raster unico: contatto ≤ 2 px,
  solidità, `frac_out`; modi test/sim/real/cp2), `tools/S1.4c_bench.R` (prep, gen, obs, mut, perf), `tools/S1.4c_decision.R`
  (verdetti, bootstrap appaiato, regola, AUC, CP-2), `tools/S1.4c_analysis.R`, `tools/S1.4c_figures.R`, `tools/S1.4c_fig_maps.py`,
  `tools/S1.4c_bio.R`, `tools/run_S1.4c.sh`, test `R/testing/test_S1.4c.R` (17), `R/testing/test_S1.4c_decision.R` (15), test Python (24); `reports/S1.4c.qmd`.
- Difetti di sviluppo D-1…D-5 (ritaglio sulla griglia, buffer negativo con precisione, sentinella, tolleranza, MU2); difetti
  dell'analisi dopo il run D-6 (eleggibilità: C-perf incluso per errore → scelta K-1 invece di K-2), D-7 (overflow dell'AUC), dichiarati.

### Risultati (`results/S1.4c/`; grezzi in `/mnt/micron/geo_spatialtrans/S1.4c/`)
- Check C: C-0 330/330, C3 0, C-NO 0, vuoti 0, C-area 5.8e-7, C-OP (maschere = R3 a scarto 0; territori osservati = S1.4 a 1e-9),
  test tutti PASS, mutanti 5/5; **C-perf WARN** (622.7 s, 5.6 GB).
- Check B (pannello {B-1, B-2, B-7, B-8}): PASS K-c 3, K-1 7, K-2 8; K-c non produce contatti (B-7 ≈ 0.005, B-8 = 0 ovunque);
  B-4 K-2 1.5–4.7 % (PASS 6/6), K-1 8.5–50 %, K-c 16–48 % in A1–A4.
- CP-0: AUC 0.95–0.97 in A1–A4 (premessa non smentita); CP-1 smentita al limite (3 ≥ 3); CP-2: S-1/S-2/S-3 dipendono dal segmentatore → descrittive.
- Previsioni P1, P2, P4, P6 confermate; P3, P5 smentite. BIO_REFERENCES B-121…B-126.

### Revisione avversariale
- 2 revisori («ricalcolo», «codice»): 0 discrepanze numeriche che cambino verdetti, 0 bug nella pipeline; decisione dipendente da D-6
  e con 1 punto su K-1; asimmetrie A4/A6 quantificate. BACKLOG BL-118…BL-127; BL-098 decisa, BL-111 superata, BL-112 risolta nel banco.

## S1.4b — 2026-10-15 — BL-098: relazione nucleo–territorio, banco delle opzioni (branch `step1-cell-layer`) — **chiusa con riserve**

### Decisioni (ROADMAP §6)
- Perimetro (Luca): banco + decisione; implementazione in S1.4c. Decisione per regola pre-registrata con veto motivato di Luca.
- **BL-098 → (c)** contenimento nel Voronoi con troncamento dichiarato (D-S1.4.2 resta), per regola: parità 4–4 → (c);
  Luca: «Nessun veto: resta (c), BL-098 riaperta sulla forma del nucleo».

### Codice (prototipo, non pacchetto)
- `tools/S1.4b_laguerre.R` (diagramma di potenza da qhull/`geometry`, celle per Sutherland–Hodgman, griglia F1 1e-6 µm;
  nuclei disgiunti in una passata), `tools/S1.4b_mutants.R` (MU1–MU4), `tools/S1.4b_bench.R` (prep, loo, mut, perf),
  `tools/S1.4b_arbiter.py` (C-L2, MU5), `tools/S1.4b_cp1.py` (CP-1 sulle maschere), `tools/S1.4b_analysis.R`,
  `tools/S1.4b_posthoc.R` (PH-C), `tools/S1.4b_posthoc_arbiter.py` (PH-A, PH-B), `tools/S1.4b_figures.R`,
  `tools/S1.4b_recheck.R`, `tools/run_S1.4b.sh`, `tools/run_S1.4b_recheck.sh`, `R/testing/test_S1.4b.R` (C-L1 19/19), `reports/S1.4b.qmd`.
- Libreria: + `geometry` 0.5.2 (PPM) in `tools/setup_r_library.R` (`renv.lock` invariato, BL-078).
- Difetti D-1…D-6 (dev notes): griglia mista, precisione geos→sf, mutanti, ricerca dei vicini O(n²) (ricalcolo 60/60 identico),
  «doppi» dell'arbitro (rieseguito identico), unità A1 × MU4 fermata.

### Risultati (`results/S1.4b/`; grezzi in `/mnt/micron/geo_spatialtrans/S1.4b/`)
- Check C: C-0 330/330, C-L1 19/19 (mutanti 4–6 FAIL ciascuno), C-L3 3.8e-15, C3 0 su 2 652 815 nuclei, C-NO 0, C-mut 5/5,
  C-perf 54.7 s / 3.6 GB; **C-L2 22/30 FAIL** (sotto-criterio delle aree mal calibrato; 0 errori puntuali su 5.3 M punti; PH-A 8/8, PH-B).
- Check B: PASS {B-1, B-2, B-3, B-7} 4 (V-a/V-b) vs 4 (V-c); B-4 FAIL ovunque tranne V-c A6; V-a/V-b troncano più di V-c
  (A1–A4 36–56 % vs 16–48 %). B-8: cerchi reali di pari area sovrapposti 33–53 % in A1–A4 (≈ B-7 reale); PH-C: ellissi dei momenti 24–47 %.
- CP-1: Laguerre con pesi reali batte 20/20 permutazioni in 30/30 ROI ((a) non smentita), guadagno modesto (m2 −4…−23 % A1–A4).
- Previsioni P1–P4 smentite, P5–P6 confermate.
- BIO_REFERENCES B-115…B-120 (contatto fra nuclei reali: cerchi ed ellissi sovrapposti; m1 Voronoi / Laguerre).

### Revisione avversariale
- 2 revisori («recalc», «code»): 0 discrepanze numeriche; decisione non robusta (semi RC-2, operatore di V-c CODE-07, C-L2 CODE-05),
  non sovrapposizione ~2× il minimo (CODE-06), bug n = 1 (CODE-01). BACKLOG BL-109…BL-117.

## S1.4 — 2026-10-15 — `derive_cell_geometry()` (branch `step1-cell-layer`) — **chiusa con riserve**

### Decisioni (ROADMAP §6)
- Perimetro (Luca): solo geometria; BL-087 in una sessione dedicata.
- D-S1.4.1 (Luca): modello di taglia a due parametri (β, σ) che contiene il design; D-S1.4.2: clamp in tutti i modi (sul raggio circoscritto).

### Codice
- `R/04b4_derive_cell_geometry.R`: `derive_cell_geometry()`, `fit_nucleus_size_model()`; interni `.dcg_*`; D-1 (sfc misto con MULTIPOLYGON) corretto prima degli stadi ufficiali; guardia «centro sul bordo» dopo la revisione (CODE-05).
- `R/testing/test_S1.4.R` (C-1 63/63), `tools/S1.4_mutants.R` (M1–M5 pre-registrati, X1–X9 della revisione; 14/14 rilevati), `tools/S1.4_bench.R` (stadi real, loo, rec, null), `tools/S1.4_perf.R`, `tools/S1.4_analysis.R`, `tools/S1.4_posthoc.R` (PH-1…PH-3), `tools/S1.4_figures.R`, `reports/S1.4.qmd`.
- Design doc v1.4 (§4.2, §4.4, §5.4, §9).

### Risultati (`results/S1.4/`; grezzi in `/mnt/micron/geo_spatialtrans/S1.4/`)
- Check C PASS tranne C-rec 31/36 (troncamento del clamp con β = 0; PH-2 36/36 senza clamp). C3: 0 violazioni su 347 971 nuclei reali e 964 660 dei nulli. C-perf 44 s, 3.6 GB.
- Modello calibrato per archetipo: β 0.11 (A1), 0.32 (A2), 0.15 (A3), 0.44 (A4), 0.62 (A5), 0.13 (A6); il design è smentito (CP-1 30/30 fold).
- Check B: B-4 FAIL A1–A5; B-1/B-2 per lo più FAIL/WARN. Causa (PH-1): i nuclei reali escono dal Voronoi dei centroidi nel 33–52 % in A1–A4.
- BIO_REFERENCES B-091…B-114 (modello di taglia, intervalli C8). BL-066, BL-006, BL-071 chiuse; nuove BL-098…BL-108.

### Revisione avversariale
- 2 revisori («recalc», «code»): 0 bug, 0 verdetti cambiati; correzioni: test e mutanti aggiuntivi, guardia del bordo, `stopifnot` nel banco, sensibilità al seme nel report (deviazione RC-8).

## S1.3d — 2026-10-14/15 — BL-084: banco di prova degli overlay (branch `step1-cell-layer`) — **chiusa**

### Decisioni (ROADMAP §6)
- Perimetro (Luca, 2026-10-14): completo F1–F6, organizzato per operazione (O1 partizione, O2 cella × bin, O3 aggregazione 2→8 µm).
- Criterio di sostituzione (Luca, 2026-10-14): garanzia documentata, poi tempo (F1 resta salvo ≥ 2× più veloce con memoria ≤ 1.5×).
- Per regola pre-registrata: **F1 (V2, snap-rounding GEOS 1e-6 µm) per O1, O2, O3**, un solo helper (BL-081); CP-3 su 300 nulli nuovi PASS.

### Codice
- `tools/S1.3d_{common,candidates,f4,compare,cases,run,perf,b1,validity,analysis,figures}.R`, `tools/S1.3d_{arbiter,arbiter_batch,arbiter_merge,arbiter_check,coverage}.py`,
  `tools/S1.3d_perf.sh`, `tools/run_S1.3d.sh`, `R/testing/test_S1.3d.R` (C-1 79/79), `reports/S1.3d.qmd`. Libreria: + `exactextractr` 0.10.1, `terra` 1.9-50 (PPM).
- Arbitro esatto di S1.3b riusato invariato (md5 `0c6a43fe…`), esteso con cella ∩ regione indicizzata e Sutherland–Hodgman esatto per O2.
- Difetti in sviluppo (6, dev_notes) e durante il run: D-1 clausola «0 invalide» estesa a F1/F3/F4; D-3 chiavi `paste()` «1e+05» (falsi errori O2 su S), validità F2 su S; corretti, decisione invariata.

### Risultati (`results/S1.3d/`; grezzi in `/mnt/micron/geo_spatialtrans/S1.3d/`)
- Check C tutti PASS (C-0 1 800/1 800; C-2 14 011 uguaglianze esatte; mutanti 4/4 + M-arb; C-4 558/558; C-5 F4 = `04b3` 1 986/1 986).
- F1: 0 errori in K, U, N, S, CP-1, H; O1 su S 82 s (metodi sf ~1 270 s). V0: O1 errato in 70/1 800 nulli, O2/O3 mai (7.85 M pezzi su S).
- F5 (polyclip 1.10.7/Clipper 6.4.0): errori topologici silenziosi in O2 su input validi (causa: Clipper 6.4.0 `AddLocalMinPoly`, corretto in 6.4.2); polyclip tronca.
- F4 (taglio unico): copertura valida 1 800/1 800, ma 14 celle invalide (strozzature della regione) e 1 660 s su S. F6: corretto, 3.2 GB, non garantito.
- B-1: 62/62 verdetti B di S1.3 invariati per ogni candidato. Previsioni smentite: F5 «0 errori»; V0 copertura invalida ≥ 80 % (61.6 %).
- Revisione avversariale: «ricalcolo» (20/20 righe, 5 467 righe V0/F1 identiche) e «esatto» (verità 100 %, 565/565 violazioni reali); 24 rilievi, nessuno cambia la decisione.
  BACKLOG BL-091…BL-097; BL-080, BL-084 chiuse.
- Stato: proposto «chiusa con riserve»; Luca (2026-10-15): «Le riserve sono pretenziose» → **chiusa**; i punti aperti restano in BACKLOG.

## S1.3c — 2026-10-12/14 — BL-077: spazio extracellulare intra-regione dai conteggi a 2 µm (branch `step1-cell-layer`) — **non chiusa**

### Decisioni (ROADMAP §6)
- UMI per bin sui geni della matrice filtrata (Luca, dopo C-3); stato «non chiusa» (BL-077 aperta, decisione 1 in vigore per difetto).

### Codice
- `tools/S1.3c_ecm.py` (stadi umi/roi/ctrl: UMI per bin, soglie di sfondo, finestre, nullo per permutazione nel territorio, modello A e BA, componenti, C-2, CP-1, CP-2, controlli sintetici C-5), `tools/S1.3c_analysis.py` (regola, tabelle; `--test` 8/8), `tools/S1.3c_posthoc.py` (PH-1…PH-4), `tools/S1.3c_figures.py`, `R/testing/test_S1.3c.py` (31/31), `reports/S1.3c.qmd`.
- Difetti: in sviluppo `m.empty` (pandas) → `is_empty`; dalla revisione M-01 (conteggio CP-1 con `numpy.bool_`), RC-2 (mappatura sulle maschere ds4), M-02 (spostamento ottimo di C-2 non applicato come pre-registrato); aggiunto R sotto il nullo (M-05). Catena rieseguita, verdetti invariati.

### Risultati (`results/S1.3c/`; grezzi in `/mnt/micron/geo_spatialtrans/S1.3c/`)
- E < 0.05 in 30/30 ROI (max 0.0294 primaria, 0.0333 sensibilità) → C formale; P2 FAIL, P1/P4 non valutabili, P3/P5 non informativi.
- CP-2 FAIL 6/6: i lumi veri non scendono al livello del vetro lontano (mediane 4.4–10.9× a 16–32 µm dal bordo); coda di diffusione anche sul vetro esterno (revisione). PH-3: recupero di E 10–97 % di un vuoto perfetto; minimo rilevabile 10–20 % (A4 tipo A > 30 %).
- CP-1 pre-registrato: A6 4/5, A2 3/5, A4 2/5, A1 1/5, A3 0/5, A5 0/5 (BA non definita). R: A5 0.59 e A6 0.63 contro il nullo 0.90 e 0.98; A1–A3 R ≈ nullo.
- C-2: 16 PASS / 4 WARN / 10 FAIL; PH-2: A2 non riproducibile fra metà (rumore), A5 0–4 µm riproducibile.
- Maschera R2: esclude neuropilo in A5 (luminosità > 0.97, ~12 % del ROI) e miocardio pallido in A6 r1–r2.
- Revisione avversariale: «Ricalcolo» (141/141 identici) e «Metodo»; 27 rilievi, azioni in `results/S1.3c/S1.3c_review_actions.csv`. BACKLOG BL-085…BL-090; BL-077 aggiornata.

## S1.3b — 2026-10-10/12 — BL-069: operazioni di sovrapposizione GEOS sui territori (branch `step1-cell-layer`)

### Decisioni (ROADMAP §6)
- Perimetro (Luca, 2026-10-10): S1.3b = solo BL-069; poi BL-077 (S1.3c), poi S1.4.
- Modello di precisione per unioni/intersezioni/differenze su territori = **V2**, snap-rounding GEOS a 1e-6 µm (`geos::geos_set_precision`, `*_prec`), per regola pre-registrata (unica eleggibile; V1/V3 esclusi da C-4). Applicazione **al momento dell'overlay** (Luca); `tessellate_voronoi()` invariato.

### Codice
- `tools/S1.3b_overlay.R` (varianti V0, V0g, V1, V2, V3; tau; misure; annodamento V3 con fusione + giunzioni a T; perturbazione), `tools/S1.3b_run.R` (stadi md5, null, grid, cp1, cp2, mutants), `tools/S1.3b_analysis.R`, `tools/S1.3b_b2.R` (verdetti S1.3 ricalcolati in worktree), `tools/S1.3b_posthoc.R`, `tools/S1.3b_figures.R`, `R/testing/test_S1.3b.R` (C-unit 57/57). `geos` 0.2.5 + `libgeos` 3.11.1 in `tools/setup_r_library.R`.
- Difetti trovati in sviluppo: nomi di lista in `s3b_rebuild()`; predicato non preparato in C-6 (> 15 min → 3 s); indici double in `factor()` ("1e+05", falso FAIL di C-6). Dopo il run: ordine di preferenza della regola, area di differenza vuota, CP-2 con eccezioni per variante. Dopo la revisione: **placebo CP-1b nullo** (RC-01) corretto in ulp e rieseguito.

### Risultati (`results/S1.3b/`; grezzi in `/mnt/micron/geo_spatialtrans/S1.3b/`)
- C-0 1 800/1 800 nulli di S1.3 riprodotti, stessi 15 falliti. V1, V2, V3: C-1, C-2, C-3, C-5, C-6 0 violazioni; C-4: V1 e V3 2 celle invalide (stesse), V2 0; mutanti 3/3.
- Grezzo: `st_union` perde celle intere (15 nulli); `st_intersection` restituisce una cella intera in ≥ 11 coppie (10 a sovrapposizione, 9 invisibili al check d'unione; 1 a lacuna trovata dalla revisione); intersezione esatta = schegge 2e-15–2.4e-14 µm² (aritmetica razionale): sbagliano le operazioni, non la geometria. GEOS 3.11.1/3.12.1/3.14.1 identici; griglia 1e-6 elimina 11/11.
- Condizione: vertice condiviso di due celle ritagliate, sul bordo della regione, diverso di 1–4 ulp (CP-1a, fusione-senza-T 15/15, test inverso 24/40 vs 0/80); non sufficiente (83 % dei nulli che passano) né specifica (**CP-1b FAIL** 18/20 dopo la correzione).
- CP-2: V0 0/180 errori silenziosi (non confermata, potenza 0.395; eccezioni solo su input invalidi); V1–V3 0 FAIL. C-6 cella × bin 2 µm: 0 errori nemmeno nel grezzo (39 nulli con la revisione). B-2: 62/62 previsioni B di S1.3 invariate; B-1 FAIL V1/V2 (errore di pre-registrazione dichiarato prima del run).
- Revisione avversariale: 2 revisori (ricalcolo indipendente 44 grandezze, aritmetica esatta, terzo GEOS), 24 rilievi (1 bug, 4 discrepanze, 8 metodo, 11 confermati); azioni in `results/S1.3b/S1.3b_review_actions.csv`.
- BACKLOG BL-078…BL-083; BL-069 decisa.

## S1.3 — 2026-10-08/09 — `tessellate_voronoi()`: territori di Voronoi per regione (branch `step1-cell-layer`)

### Decisioni (ROADMAP §6)
- D-S1.3.1 Voronoi per regione; D-S1.3.2 frammenti orfani al vicino con il confine condiviso piu' lungo; D-S1.3.3 Chaikin contenuto, `n_iter = max(1, round(3·cs))` (Luca, «approvo»). Controprova della roadmap per A6 sostituita (BL-064).
- **Motore GEOS** (`sf::st_voronoi`, `point_order = TRUE`) invece di deldir (decisione 6 del design): confronto pre-registrato (addendum C-10), regola senza eleggibili per cancelli non discriminanti, decisione di Luca sull'evidenza discriminante. Design v1.3.
- Dopo la chiusura: registrata **BL-077** (spazio extracellulare intra-regione assente per costruzione, decisione 1 del design); ordine di S1.4: BL-069 → BL-077 → `derive_cell_geometry()`.

### Codice
- `R/04b3_tessellate_voronoi.R` (nuovo): `tessellate_voronoi(centroids, regions, corner_smoothing = 0, keep_tiles = FALSE, verbose = TRUE)` → `cell_territories`, `territory_df`, `region_check`, `tiles`, `info` (frammenti, agganci, multiparte, riparazioni). Un solo motore (`.tv_engine()`).
- Test `R/testing/test_S1.3.R` (294 asserzioni, con C-4b e C-8c dopo la revisione), mutanti `tools/S1.3_mutants.R` (7), variante deldir `tools/S1.3_variants.R`, arbitro per semipiani `tools/S1.3_arbiter.R`, banco ROI `tools/S1.3_roi.R`, controprove `tools/S1.3_cp.R`, prestazioni `tools/S1.3_perf.R`, analisi `tools/S1.3_analysis.R`, figure `tools/S1.3_figures.R`, post hoc `tools/S1.3_c2_coverage.R`, `tools/S1.3_c2_pairs_check.R`, controlli della revisione `tools/S1.3_review_checks.R`, previsioni `tools/S1.3_predictions.R`, `R/testing/run_S1.3.sh`.
- Correzioni dopo la revisione (`26baaf0`): confine dei frammenti > 1e-6 µm e decisioni simultanee (contatti puntiformi, catene), donatore escluso, smussatura per parte con riparazioni contate, arbitro con certificato, C-4b/C-8c/M7.
- Difetti trovati in sviluppo: lunghezza condivisa nulla e unioni multiparte con deldir (tolleranza 1e-7 µm + aggancio); definizioni del test (esagonale, avversari vuoti); riquadri ricalcolati per ogni orfano (263 → 34 s su A1 r1); regex degli avversari nello stadio C-10a; analisi che contava come FAIL i multiparte ammessi dalla deviazione 3.

### Risultati (`results/S1.3/`; uscite grezze in `/mnt/micron/geo_spatialtrans/S1.3/`)
- Check C (GEOS) 294/294, mutanti 7/7, C-perf 35.9 s / 1 076 MB; C-5 11/40 e C-6 186/1 200 FAIL (errore di pre-registrazione, BL-074; attribuzione verificata: regola R3 39/50 identiche, 11 solo lati); C-2 su 15/1 800 nulli FAIL da artefatto di `st_union` (copertura per punti 0/2.5 M/0, BL-069/070).
- Confronto dei motori: arbitro (certificato, lati > 1e-9 µm) 16/16 GEOS vs 0/16 deldir; check C deldir 290/294 (geometria non annodata); eventi di ripiego sugli stessi ingressi 5 vs 705; C-perf deldir 4 691 s / 8 958 MB (130.7×); esponenti 0.93 vs 1.70. Regola dell'addendum: nessun eleggibile → «non chiusa»; chiusura in deroga, decisione di Luca.
- Check B 59/62 previsioni confermate (G1 23/24, G2 24/24, G3 12/14); CV delle aree fuori tolleranza in A1–A5 per G1–G3 (semina omogenea); G2 passa A6; G3 (regola di default) fallisce anche A6; sensibilita' al riferimento 48/48 identica; nessun verdetto cambiato dopo le correzioni.
- CP-1 confermata (testo del design: 93.7 % delle regioni con 3–1 000 cellule scoperte > 1 %); CP-2 senza potenza, CP-2b post hoc (Δ ≥ 0.05 da k = 2); CP-3 nessuna sovrapposizione senza contenimento ma 9 521 territori fuori dalla regione (29/30 ROI); D-1 1.071.
- Revisione avversariale: 2 revisori, 36 rilievi (codice 15: 2 bug, 5 discrepanze, 3 metodo, 5 confermati; affermazioni 21: 9 discrepanze, 9 metodo, 3 confermati); reimplementazione Python 62 121 celle a 5.3e-13; ricalcolo 188/188 verdetti; correzioni in `26baaf0`, azioni in `results/S1.3/S1.3_review_actions.csv`.
- BACKLOG BL-069…BL-076; chiuse BL-047, BL-056 (misurata), BL-061.

## R3 — 2026-10-08 — Voronoi sui nuclei reali (branch `step1-cell-layer`)

### Decisioni in apertura (ROADMAP §6)
- D-R3.1 segmentatore di R2b, sensibilità StarDist in A4/A6; D-R3.2 BL-035 come calibrazione sulle finestre R2b; D-R3.3 R3 prima di S1.3 (Luca, «approvo»).

### Codice (strumenti di riferimento, non pacchetto)
- `tools/R3_voronoi_metrics.R`: momenti esatti dei poligoni (origine locale), eccentricità/orientazione (definizione regionprops), tassellazione `deldir` ritagliata con `st_crop` + `st_intersection`, celle interne, intensità locale, riassunto — riusabile dai check B di S1.3–S1.4.
- `tools/R3_run.R` (stadi real, null, calib), `tools/R3_containment.py` (contenimento, momenti nucleari dai pixel, hed), `tools/R3_analysis.R`, `tools/R3_figures.R`, `tools/R3_posthoc.R` (P1–P3), `tools/R3_posthoc_P4.R` + `tools/R3_export_tissue.py` (P4), `R/testing/test_R3.R`, `tools/R3_mutants.sh`, `R/testing/run_R3.sh`.
- Difetto trovato da C-R3.8 e corretto prima delle analisi: momenti nucleari con E[x²] − c² (ecc ~1e-4 su nuclei isotropi) → due passate (`508a74a`).

### Risultati (`results/R3/`; per cellula, nulli e rerun in `/mnt/micron/geo_spatialtrans/R3*`)
- Check C 28/28 sintetici, 11/11 reali, mutanti 5/5, C-R3.6 1 408 file md5 identici.
- Check B 8/14 previsioni confermate; controprove: CP-1 processo omogeneo inadeguato in A1–A3, RSA(d*) troppo regolare in A4; CP-2 A6 isotropo (previsione smentita), allineamento nei tessuti densi; CP-3 dipendente dal segmentatore; CP-4 senza potenza.
- Revisione avversariale: 2 revisori (code: reimplementazione Python, 40/40 ROI identici; claims: 55/55 verdetti), 44 rilievi, azioni in `results/R3/R3_review_actions.csv`. BIO_REFERENCES B-055…B-090; BACKLOG BL-060…BL-068 (chiuse BL-028, BL-054).

## S1.2 — 2026-10-08 — `seed_centroids()`: centroidi per regione (branch `step1-cell-layer`)

### Decisioni in apertura (ROADMAP §6)
- BL-042: stride dai metadati + errore (Luca); BL-043: componenti sotto soglia alla regione adiacente (Claude su delega; confermata da B1b); D3: densità di miscela armonica 1/Σ(f_i/ρ_i) (Luca). Design v1.2 parziale (§5.2, §9).

### Codice
- `R/04b2_seed_centroids.R` (nuovo): `seed_centroids(regions, cell_types, region_composition, random_seed = 42, dmin_factor = 2/3, max_attempts_factor = 1000, max_consecutive_rejections = 10000, verbose = TRUE)` → `centroids`, `region_df` (+ target armonico, lambda, n_target, n_cells, n_failed, densità ottenuta), `type_df`, `info`. RSA a n fissato con griglia dei vicini globale (d_ij = (d_i+d_j)/2, anche fra regioni), arrotondamento Bernoulli per regione e Madow per tipo, `min_dist_um` opzionale, stato RNG globale ripristinato.
- Deviazioni approvate dal design §5.2: RSA invece di Bridson (massimale, ≈ 5× la densità target), arrotondamento stocastico, distanza globale, `min_dist_um`, preset rinviati a S1.5, `n_failed` invece del ciclo infinito.
- `R/testing/test_S1.2.R` (77 asserzioni, 8 mutanti in `tools/S1.2_mutants.R`), `R/testing/run_S1.2.sh`, `tools/S1.2_perf.R`, `tools/S1.2_checkB.R` (check B e controprove su 30 ROI R2, ~20 000 pattern simulati), `tools/S1.2_analysis.R`, `tools/S1.2_fixcheck.R`.
- Correzioni dopo la revisione avversariale: factor in ingresso → errore (RA-code-07/08); ordine dei tipi con `radix` (RA-code-09); arresto per saturazione dopo 10 000 rifiuti consecutivi (RA-code-11: striscia ruotata 1 200 s → 17–34 s); validazione di geometrie vuote e area incoerente (RA-code-16); warning di griglia grossolana (RA-code-14); documentazione di `max_attempts_factor` (RA-code-21); C-sat riportato all'input pre-registrato (RA-code-12); C-d2 + M7 (RA-code-06); D di nuovo integrato da r = 0 come pre-registrato (RA-checkB-07).

### Risultati (`results/S1.2/`; simulazioni e log in `/mnt/micron/geo_spatialtrans/S1.2/`)
- Check C 77/77 PASS, 8/8 mutanti rilevati; C-perf 449 536 punti in 4.6 s, 346 MB.
- Check B: B1a FAIL in A1, A2, A5 per tolleranza operativa sbagliata (scarto ≤ 2 cellule); B1b PASS A1–A5, FAIL A6 (0/5: semantica di esclusione di S1.1, conferma BL-043); B2 FAIL 6/6 (aggregazione 5–30 µm).
- Controprove: CP-1 4/6 PASS al limite (A1 0.266 vs 0.25; 3/6 in 3 varianti di sensibilità su 4; previsione smentita in A3); CP-2 regola ⅔·eq_r adeguata solo in A1, A2 (A5 smentito), d* 2.5–4 µm (A5 7); CP-3 6/6 (A6 al limite).
- Revisione avversariale: 2 revisori, 44 rilievi (23 code: 4 bug, 1 discrepanza; 21 checkB: 13 confermati, 8 di metodo); verifica delle correzioni 60/60 unità di simulazione identiche bit per bit. BACKLOG BL-050…BL-059; BIO_REFERENCES B-049…B-054 (d*, provvisori).

## S1.1 — 2026-10-07 — `extract_regions()`: dalla mappa dei cluster alle regioni connesse (branch `step1-cell-layer`)

### Codice
- `R/04b1_extract_regions.R` (nuovo): `extract_regions(clust, pixel_size_um, min_region_area_um2 = 100, cluster_col, simplify_tol_um = 0.5, stride = NULL)` → `region_df`, `region_polygons` (sfc POLYGON in µm), `excluded_df`, `info`. Componenti 4-connesse con union-find vettorizzato sulle run di riga (R base); poligono = unione GEOS dei rettangoli delle run (contorno sui lati dei pixel, area esatta, buchi inclusi); filtro sull'area esatta in pixel; stride rilevato (warning se > 1 e non dichiarato); `pixel_size_um` obbligatorio (C-A).
- Deviazioni dichiarate dal design §5.1 (approvate da Luca): contorno sui lati invece di Moore sui centri; union-find invece di BFS.
- `R/testing/test_S1.1.R` (716 asserzioni), `R/testing/run_S1.1.sh` (comando unico, ~3 min), `tools/S1.1_inputs.R` (57 input + 2 di prestazione), `tools/S1.1_moore_variant.R` (CP-1), `tools/S1.1_perf.R` (C11), `tools/S1.1_cp1_posthoc.R` (esplorativo), `tools/S1.1_figures.R`.
- Difetti trovati e corretti in sessione: CoverageUnion GEOS su rettangoli con giunzioni a T (poligoni non validi); `anyDuplicated()` su data.frame (24 s, ~5 GB su 16 M pixel); dopo la revisione: stride su input sparsi (RA-17, warning), etichette fuse da `as.numeric()` (RA-18), `n_holes` per non-POLYGON (RA-25), warning su tolleranza > mezzo pixel (RA-20/FX-24).

### Risultati (tracciati in `results/S1.1/`; input/output derivati in `/mnt/micron/geo_spatialtrans/S1.1/`)
- Check C: 711 PASS, 1 WARN (C6 `adv_pinch_hole`, semplificazione), 4 FAIL (CP-1) su 716; C11 8.3 s / 1.3 GB (primo giro 30.6 s / 5.2 GB, WARN, difetto corretto).
- Check B (descrittivo): sintetici 5–30 regioni, ~0 % escluso; ROI reali 263–2 104 regioni/mm², 2.7–27 % di area < 100 µm² (A6 ~22 %).
- Controprove: CP-1 FAIL (Moore sui centri −0.7…−1.4 % su I1–I4: previsione smentita; −74…+18 % sui ROI, esplorativo); CP-2, CP-3 (frazione esclusa reale 5–21 % vs nullo 51–75 %), CP-4, CP-4b PASS.
- Revisione avversariale: 26 rilievi (13 confermati, 2 discrepanze accolte, 11 di metodo) + verifica delle correzioni (114/114 esecuzioni identiche). BACKLOG BL-042…BL-049.

## R2b — 2026-09-21 → 2026-10-07 — Verità manuale per la densità nucleare (branch `step1-cell-layer`)

### Codice (tutto in `tools/`; dati derivati fuori git su `/mnt/micron/geo_spatialtrans/R2b/`)
- `R2b_windows.py` (24 finestre seedate a lato adattato, C-R2b.1/4), `R2b_annotator.py` + `R2b_annotator_core.js` (annotatore HTML autonomo alla cieca, v1.1 solo conteggio; core testato con node in `R/testing/test_R2b_core.js`), `R2b_match.py` (appaiamento punto↔oggetto 1:1, inter-rater, sommari; `--selftest` C-R2b.3), `R2b_import_model.py` (VLM → formato annotatore), `R2b_run_vlm.py` (chiamate API OpenAI-compatibili, risposte grezze conservate), `R2b_null_and_figures.py` (CP-R2b.4 nullo dei punti casuali, sovrapposizioni), `R2b_od_strata.py` (CP-R2b.5), `R2b_sheet.py`, `docs_patch_R2b.py`. `R/testing/test_R2b.sh` → `test_R2b.py` (PASS/FAIL dalle tabelle).
- Tracciati in `results/R2b/`: finestre, bersagli, pre-registrazione (con 6 aggiunte datate), annotazioni (Luca, Federica, 3 VLM) e risposte grezze dei modelli, tutte le tabelle, figure, `R2b_test_results.csv`, `R2b_adversarial_review.csv`, checklist della patologa.

### Esiti
- Check C 5 PASS + 1 WARN (C-R2b.1b: maschera bolle R2 parziale, BL-040). Check B: B-R2b.1, 2, 3a, 3b **FAIL**, B-R2b.4 PASS, B-R2b.5 non eseguito. Controprove: CP-R2b.1a WARN (Luca–Federica Δ 2.3 % sul totale, ma 6/12 finestre e 4/6 archetipi entro 10 %; A4 13 %), CP-R2b.1b FAIL (F1 3 µm 0.76), CP-R2b.2 FAIL (13–29 % dei nuclei umani persi da Cellpose **e** StarDist), CP-R2b.3 FAIL (Cellpose in A5/A6: recall basso, non precisione bassa → BL-032/033 riformulati), CP-R2b.4 PASS, **CP-R2b.5 PASS 6/6** (i nuclei persi sono più pallidi ma sopra il fondo: classe ambigua reale, tesi di Luca).
- VLM: GPT-6 Astra localizza (F1 5 µm 0.59–0.92) ma sovraconta 10–47 %; deepseek-flash non è un annotatore (thinking: non ripetibile, ×1.2–3.1; no-thinking: conteggi ricorrenti, posizioni al nullo).
- `docs/BIO_REFERENCES.md` v2: B-043…B-048 (densità come intervallo [evidenti, totale] da due annotatori); B-001/008/015/022/029/036 superate.
- Deviazioni documentate: finestre a lato adattato (non 100 µm); annotatore ridotto a solo conteggio su richiesta di Luca; area nucleare rinviata; doppia corsa DeepSeek-thinking (usata come ripetibilità); A3_w2 fuori archetipo per giudizio di Claude (la patologa non l'ha valutata).
- Chiusura: `5164bd3`.

## R2 — 2026-09-20/21 — Segmentazione nucleare su H&E e BIO_REFERENCES v1 (geometria) (branch `step1-cell-layer`)

### Ambiente
- Due venv Python (uv): `.venv` += torch 2.14 (cu130) + Cellpose 4.2.1.1 (`tools/requirements-r2-cellpose.txt` → `requirements-cellpose.lock`); `.venv-stardist` = TensorFlow 2.21 + StarDist 0.9.2 (`requirements-r2-stardist.txt` → `requirements-stardist.lock`). `tools/py_stardist.sh` mette le librerie CUDA del venv nel loader path (senza, TensorFlow non vede la GPU). `tools/setup_python_env.sh` esteso (idempotente).
- Libreria R: aggiunti `diptest`, `arrow`, `spatstat.random` (`tools/setup_r_library.R`).

### Codice (tutto in `tools/`, dati derivati fuori git su `/mnt/micron/geo_spatialtrans/R2/`)
- `R2_overview.py` (panoramiche H&E 1/16 + cluster graphclust 8 µm, griglia mm), `R2_rois.py` (30 ROI 1×1 mm da `results/R2/R2_rois.csv`, check C-R2.4, ritagli nativi, miniature), `R2_common.py` (maschera tessuto, maschera bolle BL-017, canale ematossilina, tabella nuclei, cache maschere), `R2_segment_cellpose.py` (Cellpose-SAM `cpsam_v2`, varianti rgb/hed, scale 1/2/4, riprendibile), `R2_segment_stardist.py` (2D_versatile_he, riprendibile), `R2_spaceranger_rois.py` (rasterizzazione poligoni `nucleus_segmentations.geojson` nei ROI + C-R2.5), `R2_metrics.py` (sommari, appaiamento F1@IoU0.5, effetto scala, riproducibilità), `R2_exclusive_check.py` (arbitro: OD ematossilina degli oggetti esclusivi), `R2_consensus.py` (densità di consenso), `R2_export_masks.py`, `R2_spatial.R` (g(r) con envelope CSR, dip test, Clark–Evans), `R2_figures.py` (9 figure), `R2_variant_eval.py` (rivalutazione varianti Cellpose). `R/testing/test_R2.sh` (un comando, PASS/FAIL).
- Risultati tracciati: `results/R2/*.csv`, `results/R2/spatial/`, `results/R2/figures/`, `results/R2/rois_thumbs/`, `results/R2/variant_eval/`, `results/R2/R2_preregistration.md`. Fuori git: ritagli, label image, `R2_nuclei_all.parquet` (in `results/R2/nuclei/` ignorato se > 50 MB — vedi `.gitignore`).

### Esiti
- Check C: 4 PASS, 3 WARN (C-R2.5 criterio pre-registrato mal posto, orientamento confermato dal criterio indipendente; C-R2.6 maschera bolle parziale; C-R2.7 14 % di frammenti in A5_r1). Riproducibilità bit-identica.
- Check B: B-R2.1, B-R2.2, B-R2.3, B-R2.4 **FAIL** (attese non modificate: A4 ≈ 24 500/mm² vs 1 900–6 300; area linfociti 14 µm² vs 20–45; A5 ≈ 1 300 vs 1 900–3 500; ordinamento con A3 > A5); B-R2.5 PASS; B-R2.6 WARN.
- Controprova: i due segmentatori disaccordano del 12–55 % sulla densità con direzione tessuto-dipendente (CP-R2.1 WARN); risoluzione nativa necessaria (CP-R2.3, BL-015 → BL-029); inibizione a corto raggio in A2–A5 ma non in A1/A6 (CP-R2.4); A3 distinguibile da A2 (CP-R2.5 PASS).
- Deviazione metodologica documentata: variante Cellpose-ematossilina prima liquidata come «tassellatura», poi rivalutata su richiesta di Luca e riclassificata come segmentazione nucleo+alone (BL-028).
- `docs/BIO_REFERENCES.md` v1: B-001…B-042 (6 archetipi × 7 metriche), tutte **provvisorie** (BL-024). BACKLOG BL-024…BL-031; chiusi BL-015, BL-020.

## R1 — 2026-09-18/19 — Inventario e download dei dataset Visium HD di riferimento (branch `step1-cell-layer`)

### Dati
- Inventario di 46 pagine dataset Visium HD 10x (582 file con md5 pubblicato), letto con Claude in Chrome il 2026-09-18: `data/real/inventory/datasets_10x_pages_2026-09-18.csv`, `files_10x_pages_2026-09-18.csv`.
- Scaricati 5 dataset (38 file, 79.84 GB) su `/mnt/micron/geo_spatialtrans/data_real` (symlink `data/real/datasets`), md5 = pubblicato 38/38 (download + ricalcolo indipendente): A1 `Visium_HD_Mouse_Small_Intestine` (SR 3.0.0), A2/A3 `Visium_HD_Human_Colon_Cancer` (4.0.1), A4 `Visium_HD_Human_Lymph_Node_FFPE` (4.0.1), A5 `Visium_HD_6p5mm_Mouse_Brain` (4.1.0), A6 `Visium_HD_6p5mm_Human_Heart` (4.1.0). Niente FASTQ/cloupe/molecule_info.
- `data/real/README.md` (mappa dataset → archetipi, scala, verifica istologica), `data/real/checksums.md5`, `data/real/inventory/R1_download_manifest.csv`, `R1_archetype_map.csv`, `R1_download_status.tsv`.

### Codice / ambiente
- `tools/R1_download.sh` (download riprendibile con verifica md5), `tools/R1_verify.py` (metadati BigTIFF, scala, controprova griglia, ritagli), `R/testing/test_R1.sh` (un comando, PASS/FAIL; `R1_FAST=1` salta il md5).
- Ambiente Python di progetto: `tools/setup_python_env.sh` (uv, Python 3.12, `.venv/` ignorato) + `tools/requirements-r1.txt` / `tools/requirements.lock` (tifffile, imagecodecs, zarr, pyarrow, h5py, scikit-image).
- Tabelle di verifica tracciate in `results/R1/*.csv`; ritagli in `results/R1/crops/`.

### Esiti
- Check C 30/30 PASS; controprova 1 (scala dalla griglia dei bin) 5/5 entro 0.05 %; controprova 2 (istologia) confermata da Federica Pezzuto (anatomopatologa) per A1, A2, A4, A5, A6 con due rilievi (A4 ritaglio non rappresentativo; A6 bolle di montaggio).
- Deviazioni registrate: budget disco ridefinito (79.8 GB su micron); previsione B-R1.2 smentita (esistono dataset A4 e A6); soglia B-R1.3 senza fonte (il paper di Cellpose non fissa un minimo di µm/px).
- Reperto: il tag TIFF di risoluzione è placeholder (96 dpi) in tutte le immagini, mentre l'OME-XML `PhysicalSizeX` è corretto (entro 0.08 % da `scalefactors_json`) → `estimate_pixel_size()` può usare l'OME-XML, mai il tag, e deve verificare contro la griglia dei bin.
- Revisione avversariale (2026-09-20, sotto-agente indipendente): 30/30 affermazioni confermate, 0 discrepanze, 8 rilievi non bloccanti (5 accolti, 4 in BACKLOG BL-020…BL-023).

### Documenti
- `reports/R1.qmd` + `R1.html`; `docs/sessions/2026-09-18_R1.md`; ROADMAP §1, §3, §6, §7; BACKLOG BL-014…BL-023; BIO_REFERENCES (stato); `results/R1/R1_cp2_histology.csv` (verifica istologica), `results/R1/R1_adversarial_review.csv`.

## S0 — 2026-09-18 — Infrastruttura server, baseline riproducibile, design v1.1 (branch `step1-cell-layer`)

### Server / ambiente
- Checkout `~/2025.geo_spatialtrans` portato da `dev@53103b0` (17 commit indietro) a `dev@7836829`; creato branch `step1-cell-layer`.
- Clone stantio `~/Projects/2025.geo_spatialtrans` (main, dic. 2024, remote git@ non funzionante) rinominato in `~/Projects/_archive/2025.geo_spatialtrans_2024-12_stale` (nessuna cancellazione).
- Clone locale sul Mac di Luca rimosso (sviluppo solo su server).
- Nuova libreria R di progetto `renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu` costruita da binari Posit PPM con `tools/setup_r_library.R` in tre passaggi (127 → 164 → 186 pacchetti: richiesti, namespace risolti dal sistema, chiusura ricorsiva Depends/Imports/LinkingTo); 0 compilazioni da sorgente. Include sf, deldir, polylabelr, imager, testthat, quarto e le dipendenze della pipeline attuale. Autosufficiente: nessun namespace dalla site-library di sistema.
- `.Rprofile`: imposta `RENV_PATHS_LIBRARY` su `renv/library` prima di `renv/activate.R` (il progetto ha un `DESCRIPTION`, quindi renv avrebbe usato `~/.cache/R/renv/library/...`); così `quarto render` e le sessioni R interattive vedono la stessa libreria degli script `--vanilla`.
- `renv.lock` rigenerato per R 4.6.1 (era R 4.5.0), 186 pacchetti.

### Codice
- `R/testing/full_test.R` e `full_test_visHD.R`: header con `.libPaths()` esplicito sulla libreria di progetto; nuove opzioni `--seed=<int>` e `--tag=<str>`; gli output non condividono più lo stesso path (`R/testing/full_test_result_<tag>.rds`, `results/full_simulation_data_<tag>.rds`). Nessun cambiamento alla configurazione o alla pipeline.
- Nuovi `R/testing/test_S0_library.R`, `R/testing/test_S0_baseline.R`, `R/testing/run_S0.sh` (un comando, PASS/FAIL per asserzione; `R_LIBS_SITE`/`R_LIBS_USER` inesistenti per isolamento esplicito; scrive `results/S0_checksums.md5`).
- `R/testing/full_tissue_complex.png` rigenerato dal run baseline seed 42 (effetto collaterale di `full_test.R`, BL-013).
- Tabelle di verifica tracciate in `results/S0_*.csv` (i `.rds` restano fuori git).
- RDS legacy spostati in `R/testing/legacy/` (vedi nota di provenienza `R/testing/legacy/README.md`): non erano riproducibili da `full_test.R` (prodotti da `full_test_visHD.R` e da un run a seed casuale 1200).

### Documenti
- `docs/superpowers/specs/2026-04-29-step1-cell-layer-design.md` → **v1.1**: C-A (scala px→µm obbligatoria, `estimate_pixel_size()`, `pixel_size_source`, check C8), C-B (diffusione mRNA come requisito di Step 2), pannello A1–A6 → 6 preset, refuso C6 corretto, storico §9.
- Nuovi `docs/ROADMAP.md` (v4 + stato sessioni), `CHANGELOG.md`, `docs/BACKLOG.md`, `docs/BIO_REFERENCES.md` (solo schema, zero valori), `docs/sessions/2026-09-18_S0.md`, `reports/_template.qmd`, `reports/S0.qmd` + `S0.html`.
