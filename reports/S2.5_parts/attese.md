Commit dei cancelli: `3badae6` («S2.5: perimetro, esplorazioni, funzione, cancelli»). Ok di Luca sulle previsioni P-1...P-9, sulle soglie (intervalli corretti di R7; «meno di un blocco vuoto per area» dove il vero è sempre zero), sulla regola R-2 e sulle controprove: 2026-11-09, «Ok, come proposte (Recommended)». Commit: vedi `git log` («S2.5: ok di Luca sulle previsioni»).

Semi: pannello, seme dell'utente per (area, k) estratto con Mersenne-Twister dal seme 20261109 (k = 0, 1); l'orchestratore ne estrae i semi dei quattro passi; esplorazione seme 20261110; bootstrap B = 2000, semi da 20261121. Unità indipendente: l'area (5 per tessuto), con le correzioni di R7 per le poche aree e le zone comuni. Processi: al più 6, dal file `WORKERS` (server condiviso, limite termico). Uscite in `/mnt/micron/geo_spatialtrans/S2.5/` (un file per esecuzione, rilanciabile); tabelle in `results/S2.5/`.

### Le funzioni

File `R/05b5_simulate_visium_hd.R`; batteria `R/testing/test_S2.5.R`.

1. **`simulate_visium_hd(cell_layer, seed, archetype, layers, frame_open, profiles, depth_factor, chunk_size, chunk_transcripts, component_seeds, hd, workers)`**: `build_hd_grid()` (o la griglia data, controllata contro il livello), `simulate_stripes()`, `simulate_cell_profiles()` con il fattore di profondità per M, `place_and_move_transcripts()` con la seconda coda, `capture_and_bin()`. Strati spegnibili: strisce, fondo, campo, seconda coda. Uscita come la sezione 1 del documento di progetto: `counts_2um`, `counts_8um`, `bins`, `bins_8um`, `truth$overlap` e `overlap_8um`, `origin` e `origin_8um`, `background` e `background_8um`, `arrived_2um` (trascritti arrivati per quadratino prima della cattura), `cell_counts` (generati, M volte la scala di R4), `cells` (tipo, sottotipo, generati, persi, scartati, catturati, atteso catturato), `genes`, `layers` (strisce, fondo, seconda coda, spostamento, campo, profondità), `metadata` (scala, archetipo, semi, strati, blocchi, M, tempi).
2. **`visium_hd_chunk_size(n_transcripts, chunk_transcripts = 4e6, max_cells = 2000)`**: cellule per blocco = ⌊chunk_transcripts / media dei trascritti generati per cellula⌋, fra 1 e 2 000 (regola R-1).
3. **`visium_hd_archetype(cell_layer)`**: archetipo dal preset del livello.
4. Semi: quattro estrazioni Mersenne-Twister dal seme dell'utente (strisce, profili, posizione, cattura); stato del generatore dell'utente ripristinato; risultato indipendente da `workers`.

### Il pannello

Le 30 aree dello Step 1 (livelli U6 di S1.6), griglie di S2.3 (controllate), `frame_open = TRUE`, tutti gli strati accesi, 2 semi (variante `main`). Varianti: `depth_off` (fattore di profondità 1, semi di main k = 0; CP-1), `s24_identity` (A2_r1, A6_r1 con i semi dei passi e i blocchi di 2 000 cellule del pannello di S2.4; G-4), `main_r2` (dopo la regola R-2, stessi semi di main).

### Definizioni delle misure

Stesse funzioni sui due lati (`tools/S2_5_common.py`; lato vero `tools/S2.5_real.py`, lato simulato `tools/S2.5_stats.py`).

- Quadratini: definizione **reg** (vedi E-1). Blocchi di 8 µm: almeno 8 dei 16 quadratini in reg. Misure: UMI medio, frazione di vuoti, geni rilevati medi a 2 e 8 µm; rapporto nucleo/corona (statistica di R4 PH-1, densità media nel nucleo / nella corona di 8 µm, diviso la stessa con la maschera spostata di 20 µm; funzione del pannello di S2.3, lato vero di S2.3); orientamento dei blocchi di 8 µm (E-4).
- Valore del tessuto: media delle 5 aree (simulato: prima la media dei 2 semi nell'area).
- Intervallo del vero (= intervallo di rilevanza, [DEC-061](../docs/decisions/DEC-061.md)): bootstrap percentile della media sulle 5 aree, semi-ampiezze per il fattore del tessuto `kappa_boot_bin` di R7, per il fattore dell'asimmetria se il CV fra aree supera 0.45, limite inferiore portato a 0 ([DEC-060](../docs/decisions/DEC-060.md) punti 1, 2, 6, 7).
- Intervallo del simulato: stessa costruzione con il fattore del lato simulato (riga `sim_5`, 1.632; [DEC-060](../docs/decisions/DEC-060.md) punto 3), senza fattore dell'asimmetria.
- Esito: intervallo del simulato tutto dentro l'intervallo del vero, confermata; tutto fuori, smentita; a cavallo, indeterminata.
- Intervalli del vero di ampiezza nulla (vuoti a 8 µm in intestino, tumore e cervello: 0 in tutte le aree): l'intervallo di rilevanza diventa [0, 1/n₈], con n₈ i blocchi in tessuto medi per area (meno di un blocco vuoto per area).
- Sensibilità (senza esito): le stesse misure con la definizione sr contro gli intervalli della tabella adottata di R7; il rapporto nucleo/corona anche contro la banda di ±10 % approvata in S2.3.

### Regole di decisione dei dati (scritte prima del calcolo)

**R-1, blocchi di cellule** (BL-227). `chunk_transcripts` = 4 × 10⁶ trascritti generati per blocco, dalla memoria per fase di S2.4 (E-5); la funzione di posizione e spostamento lavora a blocchi di `visium_hd_chunk_size()` cellule. I blocchi cambiano le estrazioni, non il modello. Si giudica con P-8.

**R-2, fattore di profondità sulla catena completa** (BL-205, D2-9). S2.2 ha tarato il fattore senza le perdite dello spostamento e il fondo dei valori veri di R5, e ha rinviato la verifica a questo step. Dopo il pannello `main`, per ogni tessuto: f' = f · (U_vero − b) / (U_sim − b), con U l'UMI medio per quadratino di 2 µm in tessuto (definizione reg; U_sim del pannello main) e b il fondo della cattura del tessuto (UMI per quadratino; la media di c sul tessuto è 1). Un passo solo, in tutti e sei i tessuti; si rifà il pannello (`main_r2`) e si rigiudicano P-1...P-9. Le previsioni del pannello main dicono come si comporta il simulatore di oggi; quelle di `main_r2` sono il verdetto sul preset che resta. Se dopo R-2 l'UMI non è confermato in un tessuto, nessuna seconda iterazione: dubbio.

### Cancelli (al massimo 5)

Soglia per i confronti multipli: |z| < Φ⁻¹(1 − 0.0005/m), m = numero di confronti.

| id | controllo | tolleranza | comando |
|---|---|---|---|
| G-1 | **batteria di test** (35 asserzioni): (a) equivalenza con la catena a mano (due tessuti, `frame_open`, griglia costruita); (b) identità delle uscite: `counts_8um`, `origin_8um`, `overlap_8um`, `background_8um` uguali a `sum_to_8um()` di S2.1; UMI per quadratino = origine + fondo; generati = catturati + scartati + persi; arrivati = generati − persi; atteso catturato = generati / M; metadati; (c) ogni strato spento uguale alla catena a mano con lo strato spento; `depth_factor`; (d) regola dei blocchi; (e) semi: stato dell'utente ripristinato (due generatori), stesso seme stesso risultato, semi dei passi = 4 estrazioni Mersenne-Twister, 1 e 2 processi identici; (f) errori espliciti e archetipo dal preset. Almeno 10 mutanti (profondità senza M, semi in progressione, blocchi senza M, blocco di 8 µm sfasato nell'origine, seconda coda ignorata, fondo ignorato, `frame_open` non passato, griglia non controllata, generatore non ripristinato, atteso catturato senza M, frazioni a 8 µm su 4 µm²), tutti rilevati | identità esatte; mutanti 100 % | `R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla R/testing/test_S2.5.R`; `.venv/bin/python tools/S2.5_mutants.py` |
| G-2 | **invarianti del pannello** (main, depth_off, main_r2): conservazione per cellula, gene, quadratino; somma a 8 µm; verità a 8 µm; nessun NA; catturati entro la soglia di Bonferroni dall'atteso Σ arrivati · c/C; fondo dall'atteso b · Σ c | identità; z sotto la soglia | `Rscript --vanilla tools/S2.5_panel.R check <dir>` |
| G-3 | **riproducibilità**: lato vero, pannello main k = 0, misure e tabelle degli esiti rifatti in un'altra cartella con un numero di processi diverso | md5 identici (JSON senza i tempi) | `tools/S2.5_repro.sh` |
| G-4 | **integrità e continuità**: (a) md5 di ogni file letto prima e dopo (matrici e posizioni delle 5 sezioni, tabelle di R2, R4, R7, livelli di S1.6, griglie e strati esportati di S2.3, lato vero di S2.3, profili di S2.2, pannello di S2.4); (b) batterie di S2.1, S2.2, S2.3, S2.4 invariate; (c) A2_r1 e A6_r1 con i semi e i blocchi di S2.4: totali per quadratino identici al pannello di S2.4; (d) lato vero sr = `R4_roi.csv` (6 × 30) | identità; (d) entro 1e-9 | `tools/S2.5_md5.sh before / after`; batterie; `tools/S2.5_identity.sh` |
| G-5 | **revisione avversariale**: sotto-agenti con codice proprio ricalcolano (a) le misure vere (reg) di un'area per tessuto dalle matrici raw; (b) le misure simulate di 6 esecuzioni dalle matrici salvate; (c) intervalli corretti ed esiti dalle tabelle per area; (d) la regola R-2; (e) l'orientamento nel cuore (vero e simulato); (f) leggono la funzione contro la sezione 1 del documento di progetto e questo file | quantità deterministiche entro 1e-9 relativo; ogni discrepanza è un dubbio | sotto-agenti, report |

### Previsioni

Intervallo di rilevanza: l'intervallo corretto del valore vero del tessuto ([DEC-061](../docs/decisions/DEC-061.md)), definizione reg (decisione chiesta a Luca, sotto). Esiti per tessuto, sul pannello main e dopo R-2.

| id | quantità | tessuti | atteso (pannello main) | atteso dopo R-2 | 2 errori standard tipici |
|---|---|---|---|---|---|
| P-1 | UMI medio per quadratino di 2 µm | A1-A6 | smentita per difetto nell'intestino (circa 25 contro 26.2-32.3) e nel cervello (circa 39 contro 40.4-47.8); confermata negli altri quattro (E-3) | confermata 6/6 | simulato ±1-3 % (E-3); vero ±10 % (intestino) ... ±50 % (stroma) |
| P-2 | frazione di quadratini di 2 µm vuoti | A1-A6 | smentita per difetto nell'intestino e nel tumore (circa un decimo del vero, BL-226); confermata in stroma, linfonodo, cervello, cuore (intervalli larghi) | invariata (la profondità sposta poco i vuoti) | vero ±30-100 % |
| P-3 | geni rilevati per quadratino di 2 µm | A1-A6 | come P-1 (i geni seguono gli UMI: E-2 0.96 nell'intestino, 0.86 nel cervello) | confermata 6/6 | come P-1 |
| P-4 | UMI medio per blocco di 8 µm | A1-A6 | come P-1 | confermata 6/6 | come P-1 |
| P-5 | frazione di blocchi di 8 µm vuoti | A1-A6 | confermata 6/6 (simulato 0 in E-2; veri ≤ 0.008) | confermata 6/6 | — |
| P-6 | geni rilevati per blocco di 8 µm | A1-A6 | smentita nel cervello (circa 500 contro 514-615); indeterminata o confermata nell'intestino (333 contro 324-385); confermata negli altri | confermata 6/6 | come P-1 |
| P-7 | rapporto nucleo/corona | A1-A6 | smentita per eccesso nel tumore del colon (0.91 contro 0.82-0.88 in E-2); indeterminata nell'intestino (0.95 contro 0.92-0.95); confermata negli altri (cuore con intervallo largo, 1.52-2.02) | invariata | vero ±1.5 % (intestino) ... ±14 % (stroma, cuore) |
| P-8 | tempo e memoria per area, un processo | A1-A6 | confermata 6/6: ≤ 300 s e ≤ 8 GB (E-2: 20-166 s, 2.9-6.7 GB) | confermata 6/6 | da misurare; server condiviso |
| P-9 | orientamento del campo dei blocchi di 8 µm | A1-A6 | smentita per difetto nel cuore (CP-2); indeterminata o smentita per eccesso in intestino, tumore e linfonodo (contrasto dei bordi); confermata in stroma e cervello | invariata | vero: cuore 0.09-0.41, intestino 0.10-0.24, altri da 0 a 0.18 |

P-8: intervallo di rilevanza 0-300 s e 0-8 GB (bilancio del documento di progetto, sezione 12), media delle 5 aree ± 2 errori standard come in S2.4.

### Controprova

- **CP-1, il confronto vede un simulatore sbagliato?** Pannello con il fattore di profondità tolto (1, un seme): UMI per quadratino di 2 µm smentita nel linfonodo (f = 0.79: circa 14.0 contro 8.8-13.6), nel cervello (1.25: circa 31 contro 40.4-47.8) e nel cuore (2.62: circa 1.2 contro 2.1-5.0). Se uno dei tre non risulta smentito, lì il confronto non ha potere. Intestino, tumore (f vicino a 1) e stroma (intervallo del vero largo) non dicono nulla.
- **CP-2, il cuore come controllo negativo** (documento di progetto, sezione 15). Il simulatore ha cellule isotrope; i conteggi veri del cuore hanno un orientamento che viene dall'interno del tessuto (E-4). P-9 deve risultare smentita per difetto nel cuore. Se non lo è, la misura non vede un limite noto del simulatore.

### Misure (senza intervallo)

- Blocchi scelti e fattore M per area; tempi per passo.
- BL-212: quadratini in tessuto senza territorio, per area (attesa: meno di un millesimo).
- Sensibilità: esiti con la definizione sr e la tabella adottata di R7; rapporto nucleo/corona contro ±10 %.
- Fattori di profondità prima e dopo R-2.
