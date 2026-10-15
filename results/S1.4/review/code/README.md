# S1.4 — revisione avversariale: CODICE e METODO (revisore `code`)

Commit esaminato: `ec2a874` (branch `step1-cell-layer`). Data: 2026-10-15. Nessun file versionato è stato modificato.
Script e output: `/mnt/micron/geo_spatialtrans/S1.4/review_code/` (adversarial.R, angle_check.R, loo_recompute.R, round2.R, round3.R, review_mutants.R, test_S1.4_review.R, verd.py; CSV e log in quella cartella).

## Sintesi

Revisione avversariale del codice e del metodo di S1.4 (commit ec2a874). Nessun bug che cambi numeri o verdetti. Dove possibile ho ricalcolato con stimatori indipendenti (lm, stats::ks.test, geos_within+STRtree, momenti per anello): le stime LOO, il metro inter-ROI, le KS simulate (40 fold), i semi, l'appaiamento BL-060 e i verdetti B/CP coincidono a 1e-16. Il clamp garantisce C3 per costruzione (cerchio; ellisse: limite inferiore ammissibile, ricerca monotona per convessità, bisezione massimale). Su 64 casi avversari (38 694 nuclei) ci sono 0 violazioni. La convenzione degli angoli pacchetto/R3 è identica: Δ dell'asse 0.0°; con il segno invertito l'allineamento nucleo-territorio scompare. Rilievi principali: (1) i test non hanno potere su 5 difetti plausibili su 9, fra cui il margine 1/0.85 dell'ellisse e la bisezione (D-S1.4.2), l'ordine del flusso casuale e i buchi nell'asse; (2) la soglia di CP-2 non è calibrata: il FAIL di A5 ha probabilità 0.47 anche con un pipeline corretto; (3) 'M1/M3 basta' in CP-1 dipende in parte dal rumore di un solo seme; (4) C5 non è garantito con centro = generatore sul bordo (area 0 che st_within accetta), caso assente nei dati; (5) C-area dipende dalla scala delle coordinate (fino a 1.4e-7 a 1e4 µm); (6) rho/A_ref incoerenti nello stadio null (solo check strutturali). Le scelte non fissate dalla pre-registrazione (definizione di 'batte', M0 con rho stimato, B-4 su tutte le cellule) sono dichiarate nelle dev notes. B-4 sull'insieme di valutazione dà gli stessi verdetti.

Conteggio: bug 0, deviazioni 0, rilievi 9, conferme 13. Nessun rilievo cambia un verdetto pre-registrato.

## Rilievi e conferme

### CODE-01 — conferma

**Cosa.** Modello di taglia conforme a D-S1.4.1 (log A_N = log(rho^2 A_ref) + beta(log A_T - log A_ref) + sigma z; modo assoluto log(pi r^2) + sigma z; A_ref = 1e6/density).

**Evidenza.** R/04b4 righe 176-183 (.dcg_target_area). fit_nucleus_size_model righe 92-106: 'design' -> rho = media geometrica di sqrt(A_N/A_T) (verificato algebricamente). Ricalcolo indipendente con lm() sui 40 fold LOO: |Δrho|, |Δbeta|, |Δsigma| <= 4.4e-16 (loo_recompute.csv).

### CODE-02 — conferma

**Cosa.** Fattore del poligono corretto: area del 32-gono inscritto nell'ellisse (k a s, k b s) = pi a b s^2 = s^2 A_target, con k = sqrt(2pi/(n sin(2pi/n))) = 1.00322.

**Evidenza.** Riga 185; derivazione: (n/2) k^2 a b sin(2pi/n) = pi a b. Area GEOS vs nucleus_area: max rel 2.0e-13 sui reali (verdicts), <= 1.2e-8 nei miei 64 casi avversari (vedi CODE-09).

### CODE-03 — conferma

**Cosa.** Clamp del cerchio garantisce C3 per costruzione; per l'ellisse il limite inferiore (cerchio circoscritto, s_circ = 0.85 d/(k a)) è sempre ammissibile e la ricerca è monotona.

**Evidenza.** Riga 231: k·a·s <= 0.85 d => poligono nel disco di raggio 0.85 d < d attorno al centro, disco contenuto nel territorio (d = distanza dal bordo, buchi e parti inclusi). Monotonia: l'ellisse-poligono è convessa e contiene il centro, quindi s'P ⊂ sP per s' < s: fits(s) non crescente. Empirico (ellipse_clamp_optimality.csv, 12 config su A4 p1/A1 r1/A2 r1, aspetto 1.6/3, orientazione territorio/casuale, 15 213-30 140 ellissi clampate per config): s(1+1e-6) non entra nel 100 %, il dilatato (1-1e-9)/0.85 entra nel 100 % (le 7+2 'violazioni' a 1/0.85 esatto hanno distanza 0: contatto, arrotondamento); la bisezione guadagna mediana +2..+13 % su s_circ, mai ferma al limite inferiore.

### CODE-04 — conferma

**Cosa.** Funzione non rotta da territori avversari: 0 violazioni C3 (sf::st_within) e 0 errori su 64 casi / 38 694 nuclei: schegge h = 1e-2..1e-4 µm (anche a 1e4 µm), celle da 0.09 e 1e-4 µm², buchi (anche anelli orari), MULTIPOLYGON con buco e generatore nella parte piccola, pettine e stella non convessi, Voronoi 3000 celle a 1e4 µm con aspetto 1/1.6/20, orientazione casuale/territorio/fissa, centro POI/generatore, sigma 0.8, beta 1.5.

**Evidenza.** /mnt/micron/geo_spatialtrans/S1.4/review_code/adversarial.R -> adversarial_geometry.csv; frazione d'area del nucleo fuori dal territorio (st_difference) = 0 in tutti i casi.

### CODE-05 — rilievo

**Cosa.** C5 non garantito in modalità 'generator' se il generatore cade sul bordo del territorio (d = 0): s = 0, area esattamente 0, e sia geos_within sia sf::st_within dichiarano 'dentro' il poligono degenere (C3 non lo vede). Nessuna guardia nel codice.

**Evidenza.** Caso 'generatore sul bordo' in adversarial_geometry.csv: min_area = 0, within_sf 2/2. Sui 40 set reali min d_gen = 0.031 µm (A6_r3_spaceranger) > 0 (real_geometry_census.csv): non si verifica nei dati di S1.4; i nuclei del B-5 'gen' hanno area > 0.

### CODE-06 — conferma

**Cosa.** Flusso casuale come dichiarato: set.seed(random_seed); z <- rnorm(n); u <- runif(n) in ordine di cell_id, sempre (anche sigma = 0); angolo casuale pi·u - pi/2 nella convenzione (-pi/2, pi/2).

**Evidenza.** R/04b4 riga 47 e 189; cen ordinato per cell_id (riga 117); T14 'righe permutate' PASS.

### CODE-07 — conferma

**Cosa.** CONVENZIONE DEGLI ANGOLI identica fra pacchetto e R3 (y verso il basso, da +x verso +y, 0.5·atan2(2mxy, mxx-myy), x = colonna·upp, y = riga·upp): i Δθ di B-5 sono corretti.

**Evidenza.** Codice: R3_containment.py pixel_moments (x = c*upp, y = r*upp, th = 0.5*arctan2(2mxy, mxx-myy)); R3_voronoi_metrics.R r3_shape idem; .dcg_territory_axis idem. Dati (angle_convention.csv, 40 set): asse del territorio del pacchetto vs theta_T di R3 sulle interne valutate: mediana e q95 della differenza assiale = 0.0°, con segno invertito 40-47°; centroide pixel del nucleo vs generatore 0.0 µm (assi scambiati: 354-519 µm). Controprova: Δθ nucleo-territorio reale mediano 29-41° (A1-A5) con la convenzione usata, 44-46° (≈ uniforme) con il segno invertito di theta_N: l'allineamento esiste solo nella convenzione coerente.

### CODE-08 — conferma

**Cosa.** .dcg_territory_axis corretto con buchi e multiparte (segno dei buchi, anelli orari, momenti combinati delle parti).

**Evidenza.** Confronto con un calcolo indipendente (r3_poly_moments per anello, buchi sottratti) sui 291 territori reali con buchi (236) o multiparte: max differenza 5.1e-10 gradi (holes_axis_check.csv). Il theta_T di R3 differisce su questi territori (fino a 89.7°) perché R3 usa le tessere deldir non ritagliate: 288/291 non sono interne, le 3 interne con buco differiscono di 10-19°; B-5 usa l'asse del pacchetto, quindi nessun impatto.

### CODE-09 — rilievo

**Cosa.** C-area (soglia 1e-9 relativa) dipende dalla scala delle coordinate: con coordinate ~1e4 µm e nuclei minuscoli o aspetto 20 l'errore relativo area GEOS / nucleus_area arriva a 1.4e-7 (cancellazione nella shoelace in coordinate assolute).

**Evidenza.** adversarial_geometry.csv: scheggia h=1e-3 @1e4 ellisse ortogonale 1.4e-7; minuscole 1e-4 µm² @1e4 ellisse20 1.2e-8; a coordinate <= 1000 µm (tutti i ROI reali, max_coord 1000.1) l'errore resta <= 2e-13. Non tocca i verdetti di S1.4; da considerare per sezioni grandi (S1.3d, 4000 µm, dove C-area non è misurata).

### CODE-10 — conferma

**Cosa.** LOO senza fuga d'informazione: A_ref, (rho, beta, sigma) di M0-M3 e aspetto dell'ellisse B-5 stimati solo sui 4 ROI di training; la simulazione usa i territori veri del ROI escluso (come pre-registrato).

**Evidenza.** loo_recompute.csv (40/40 fold): A_ref salvato = Σ A_T / n dei 4 training (rel. 2.4e-16), diverso dall'A_ref a 5 ROI (es. A5 fold 1: 1126.24 vs 961.83); stime M2 ricalcolate con lm() identiche (<= 4.4e-16); e_med(ecc_N) training identico (0); roi_test = ids[k] 40/40; n_ev identico 40/40.

### CODE-11 — conferma

**Cosa.** Insiemi di valutazione, metro inter-ROI e semi corretti; le KS e |ΔSpearman| simulate si riproducono esattamente con uno stimatore indipendente (stats::ks.test).

**Evidenza.** loo_recompute.csv: metro (4 ROI × 3 statistiche) max |Δ| 1.1e-16; ri-simulazione M0-M3 con parametri e seme salvati: max |Δ| di ks_nc, ks_an, d_sp, clamp 2.2e-16; seme = 20261015 + 1000·i + 10·k 40/40.

### CODE-12 — conferma

**Cosa.** Appaiamento BL-060 corretto: nessun centroide StarDist su un confine condiviso, nessun duplicato, nessun centroide fuori dai territori; ordine idx dei parquet StarDist gen/nuc identico (non verificato dal banco, ma vero nei dati).

**Evidenza.** bl060_recheck.csv (10 set A4/A6): appaiamento rifatto con geos_within + STRtree: 0 disaccordi, 0 discrepanze di area e di theta_N; frazione di interne appaiate A4 0.708-0.759, A6 0.472-0.770 (come dev notes). tools/S1.4_bench.R riga 46-52: manca stopifnot(identical(sg$idx, sn$idx)) per il ramo spaceranger (rilievo minore di robustezza).

### CODE-13 — rilievo

**Cosa.** Selezione BL-060 non neutra in A4: le interne NON appaiate hanno territori più piccoli delle appaiate (mediana A_T 26.7-35.3 vs 31.4-36.5 µm²), quindi stima e valutazione di A4 escludono preferenzialmente le cellule piccole (dove il clamp morde di più).

**Evidenza.** bl060_recheck.csv colonne at_med_matched / at_med_int_unmatched. Minaccia già dichiarata (CP-3); il dato la quantifica.

### CODE-14 — conferma

**Cosa.** Regole di verdetto applicate come nelle dev notes; ricalcolo indipendente dei verdetti B-1/B-2/B-3 (fold e archetipo), B-4, CP-1 dalla tabella dei fold: tutti identici a S1.4_verdicts.csv.

**Evidenza.** verd.py su S1.4_loo_folds.csv: es. A3 M2 B-1 PPPPP PASS, B-2 WWWWF WARN, B-3 WPWWW WARN; CP-1 M2>M0 5/5 in 6/6; M1 non peggiore 0,2,0,1,4,1; M3 0,0,2,1,0,3 (A1..A6).

### CODE-15 — conferma

**Cosa.** Le scelte non definite nella pre-registrazione sono dichiarate nelle dev notes prima dei run (nessuna deviazione non dichiarata trovata): 'batte' = S = KS(N/C)+KS(A_N) più basso; M0 con rho stimato (media geometrica); B-4 su tutte le cellule; regole WARN/FAIL d'archetipo; CP-1 'smentita' equivalente alla pre-registrazione (>= 5/6 ⇔ < 4/5 in al più 1).

**Evidenza.** S1.4_dev_notes.md sez. Banco; tools/S1.4_analysis.R righe 22, 36, 101, 139-146. B-4 sull'insieme di valutazione invece che su tutte: A1 0.176, A2 0.336, A3 0.173, A4 0.438, A5 0.057, A6 0.0035 -> stessi verdetti (FAIL ×5, PASS A6).

### CODE-16 — rilievo

**Cosa.** M0 'design' ha rho calibrato sui dati (non il rho di catalogo del design): è il design nella forma, con taglia media corretta. È dichiarato, ma l'etichetta 'M0 design' in P2/CP-1 va letta così; rende P2 e CP-1 conservativi (un M0 col rho di catalogo sarebbe peggiore).

**Evidenza.** fit_nucleus_size_model(model='design'): b = 1, a = mean(y - x), sigma = 0 (righe 101-104); dev notes 'M0 = design (β = 1, σ = 0, ρ = media geometrica ...)'.

### CODE-17 — rilievo

**Cosa.** CP-2: la soglia pre-registrata |beta_perm| <= 0.02 non è calibrata sull'errore standard: in A5 (n = 1993) e A6 (n = 2221) la SD di beta sotto permutazione è 0.027 e 0.017, quindi un pipeline corretto fallisce CP-2 con probabilità 0.47 (A5) e 0.26 (A6); P(6/6 PASS | nessuna fuga) ≈ 0.39. Il FAIL di A5 (0.0224, SE 0.0218) non è evidenza di difetto.

**Evidenza.** cp2_permutations.csv: 500 permutazioni per archetipo; con il seme ufficiale si riproduce beta_perm A5 = 0.0224, A6 = -0.0164. Il verdetto pre-registrato resta FAIL (non si cambia), ma va interpretato come non informativo in A5/A6. Nota: S1.4_verdicts.csv non ha una riga complessiva CP-2 (6/6 richiesto: esito 5/6).

### CODE-18 — rilievo

**Cosa.** CP-1, conclusioni 'M1 basta' / 'M3 basta': il confronto si fa con un solo seme e diversi scarti sono dentro il rumore Monte Carlo (SD di S fra 10 semi = 0.002-0.047). 'M3 basta' in A6 poggia su 3 fold con |z| < 0.3; in A5 'M1 basta' 3 dei 4 fold hanno |z| < 0.8; viceversa in A1 M3 è 'peggiore' in 5/5 fold ma solo di ~1 SD di seme (z 0.89-1.17). La conclusione principale (M2 batte M0) è robusta: z >= 6.6 in 29/30 fold (A6 fold 1: z = 1.33).

**Evidenza.** cp1_seed_noise.csv: z = (S_Mk - S_M2)/(sqrt(2)·SD_seme(S_M2)). Le regole sono quelle pre-registrate; si suggerisce di riportare la media sui 10 semi accanto.

### CODE-19 — rilievo

**Cosa.** Test e mutanti: i 5 mutanti pre-registrati sono rilevati (riprodotto: FAIL 5,2,2,5,9), ma 5 difetti plausibili su 9 aggiuntivi passano 57/57: X1 ellisse senza bisezione (solo cerchio circoscritto), X2 ellisse senza il margine 1/0.85 (nucleo a contatto del bordo, viola D-S1.4.2), X4 ordine del flusso u prima di z (viola la dichiarazione del flusso), X5 buchi non sottratti nell'asse del territorio, X8 angolo casuale in [0, pi). Rilevati: X3 segno dell'asse, X6 clamp del cerchio a 1.0 d, X7 area non scalata, X9 POI grossolano.

**Evidenza.** /mnt/micron/geo_spatialtrans/S1.4/review_code/review_mutants.R + test_S1.4_review.R (copia del test con RES nella mia cartella), log in tests/log_*.txt. Es. T5 (ellisse ortogonale): s = 0.1502 corretto, 0.1502 con X1, 0.1767 con X2, e T5 passa in tutti e tre. Il codice attuale è corretto su questi punti (CODE-03, CODE-06, CODE-08), quindi nessun numero cambia; il difetto è di potere dei test.

### CODE-20 — rilievo

**Cosa.** Stadio null: il rho stimato con l'A_ref dei 5 ROI reali è usato con A_ref = 1e6/densità del nullo; con beta != 1 l'intercetta si sposta di (1-beta)·log(A_ref_nullo/A_ref_reale). Non dichiarato nelle dev notes.

**Evidenza.** tools/S1.4_bench.R riga 176: ctab(fs$fit, 1e6 / dens). Effetto stimato <= ~6 % sull'area (A5: beta 0.62, A_ref per ROI 894-1126 vs 962). Riguarda solo C3/C5/C-area/C7 sui nulli (strutturali): nessun verdetto cambia.

### CODE-21 — rilievo

**Cosa.** Numero di core: tools/S1.4_bench.R usa di default 40 core (S14_CORES) e tools/S1.4_posthoc.R mc.cores = 30, sopra il limite di 20 indicato per il server.

**Evidenza.** tools/S1.4_bench.R riga 16; tools/S1.4_posthoc.R righe 41, 59. Nessun effetto sui numeri (mclapply con semi espliciti per job).

### CODE-22 — conferma

**Cosa.** C-rec FAIL (31/36) è dovuto al troncamento del clamp nelle configurazioni beta = 0 (clamp 5-10 %, precondizione 'clamp ≈ 0' della pre-registrazione non soddisfatta), non allo stimatore.

**Evidenza.** S1.4_rec.csv: i 5 FAIL sono tutti beta = 0 (A3 σ0.4, A5 σ0.2/0.4, A6 σ0.2/0.4) con clamp 0.054-0.096; S1.4_ph2.csv (senza clamp) 36/36, max |Δβ| 0.029 (letto dal CSV, non ricalcolato).

## Ricalcoli

| Quantità | Riportato | Ricalcolato | Coincide |
|---|---|---|---|
| A_ref, rho, beta, sigma (M0-M3), 40 fold LOO | S1.4/loo/*.rds fits | lm() indipendente sui 4 ROI di training: max |Δ| 4.4e-16 | sì |
| Metro inter-ROI (KS N/C, KS A_N, |ΔSpearman|), 40 fold × 4 ROI | loo rds yard | stats::ks.test + cor: max |Δ| 1.1e-16 | sì |
| KS/ΔSpearman/clamp simulati M0-M3, 40 fold | loo rds sims | ri-simulazione con seme salvato: max |Δ| 2.2e-16 | sì |
| Semi dei fold | 20261015 + 1000·i + 10·k | 40/40 conformi | sì |
| Verdetti B-1/B-2/B-3/B-4 e CP-1 per archetipo (M2, primario) | S1.4_verdicts.csv | ricalcolati da S1.4_loo_folds.csv: identici (es. A3 B-1 PASS, A6 B-1 WARN, B-4 A6 0.0042 PASS) | sì |
| CP-2 beta_perm A5 / A6 | 0.0224 / -0.0164 | 0.0224 / -0.0164 (seme 20261015+i); P(|β|>0.02) sotto nullo 0.47 / 0.26 | sì |
| Appaiamento BL-060, territori appaiati (A4 f1, A6 r1) | 16706 / 325 (frazione interne 0.708 / 0.472) | 16706 / 325, 0 disaccordi in 10 set | sì |
| Asse del territorio pacchetto vs theta_T R3 (interne valutate, 40 set) | convenzione dichiarata identica | Δ mediano 0.0°, q95 0.0°; segno invertito 40-47° | sì |
| Mutanti C-mut (FAIL per mutante) | 5,2,2,5,9 | 5,2,2,5,9 (copia del test) | sì |
| B-4 frazione clampata su insieme di valutazione (non dichiarata, controprova) | tutte le cellule: A1 0.163 A2 0.328 A3 0.174 A4 0.477 A5 0.065 A6 0.004 | valutazione: 0.176 0.336 0.173 0.438 0.057 0.0035 (stessi verdetti) | sì |

## Limiti di questa revisione

- Non ho ricalcolato in modo indipendente gli stadi rec, null e perf, né PH-1/PH-3. Di questi ho letto i CSV; di PH-2 ho letto solo l'esito (36/36).
- Ho verificato l'ottimalità della bisezione su 3 ROI reali (A4 p1, A1 r1, A2 r1) × 4 configurazioni, non su tutti i 40 set.
- I mutanti aggiuntivi (X1-X9) li ho scelti io: sono un campione di difetti plausibili, non un'enumerazione completa.
- Le figure (tools/S1.4_figures.R) non le ho riviste.
