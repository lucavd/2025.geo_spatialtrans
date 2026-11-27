# Revisione indipendente dei numeri di G1 (geo_spatialtrans, polmone 10x)

Revisore: sotto-agente senza storia, codice proprio in `/mnt/micron/geo_spatialtrans/G1_review/rev_indip/` (`rev_D.py`, `rev_M.py`, `rev_cmp.py`, `rev_src.py`, `rev_A.py`). Repo `~/2025.geo_spatialtrans` (ramo `riesame`, HEAD 17010b3) solo letto; nessuno script del progetto eseguito.

## In sintesi

1. **Tutti i numeri dichiarati per G1 tornano.** D, intervalli ed errori standard dei 13 bracci, nei 3 gruppi, coincidono con `results/sonda/S6/D.json` con scarto 0 (39 righe). Lo stesso schema di semi (default_rng(20261123), B = 2000, metodo (ii), f = 1,5838) dà gli stessi estremi.
2. **Il metro regge anche dal basso.** Ho ricalcolato M-1...M-6 per tutte le 10 aree (5 di tumore, 5 di stroma) dalle uscite della catena per quadratino, cella e blocco: 80 esecuzioni (V, T omologo di H1, B con seme 0 e 1, Y, Yalt, YH1, Yshuf). Ho ricalcolato da zero anche le etichette trasferite dal riferimento HLCA: 153 303 su 153 303 identiche. Misure e conteggi coincidono con le tabelle per area con scarto 0.
3. **I bracci Yalt e YH1 non copiano l'area misurata.** In 30 unità su 30 la sorgente dei tipi è quella dichiarata (Yalt: area k mod 5 + 1 di H2; YH1: area omologa di H1 in `S3b_sonda`; Yshuf: l'area stessa, permutata). Il tipo simulato è uguale al tipo trasferito della cellula sorgente in 100 % delle cellule di Yalt e YH1. I nuclei sono gli stessi di Y.
4. **G1, regola di RIESAME_C par. 6: caso 4 con δ = 1 h, caso 1 con δ = 0,5 h**, come dichiarato. Con δ = 1 la stima puntuale di Yalt e di YH1 sta nella fascia «riposizionare». L'intervallo di Yalt tocca invece tutte e tre le fasce; quello di YH1 sta fra «riposizionare» e «chiudere».
5. **A5 e A4 tornano**: A5 sul guadagno −0,600 ± 0,066 (non regge), tau-b di CP-gain 0,6; A4 con la catena U e = +3,108 ± 0,396, con la catena C e = +1,107 ± 0,485 (non regge con nessuna delle due).
6. Nessuna discrepanza. Ho trovato alcune anomalie nei documenti e nel registro (sezione 6); nessuna cambia un numero di G1.

## 1. D dai valori per area (compito 1)

File: `results/sonda/S3b/area_values.csv` (V, T, B, P), `results/sonda/S4/area_values.csv` (Y, G, U, A, AE), `results/sonda/S6/area_values.csv` (Yalt, YH1, Yshuf), `results/sonda/S3b/h.csv`.
Metodo (scritto da zero da REGOLA_S3c): per cella (tessuto × misura) d_k = braccio_k − V_k sulle 5 aree appaiate; q = (media(d)² − var(d, ddof 1)/5)/h²; D = radice con segno della media di q sulle 12 celle (M-4 tenuta in entrambi i tessuti: media dei singoletti di V 0,812 nel tumore, 0,775 nello stroma, sotto 0,99). B = media dei due semi per area. Bootstrap per coppie stratificato per tessuto, stessi indici per V e per tutti i bracci; intervallo percentile spostato della distorsione mediana, poi semi-ampiezze × f.

Tumore + stroma:

| braccio | D dichiarato | D ricalcolato | intervallo dichiarato | intervallo ricalcolato | se ricalcolato | scarto massimo |
|---|---|---|---|---|---|---|
| T | 0.336461 | 0.336461 | (-0.2627; 0.8210) | (-0.2627; 0.8210) | 0.2768 | 0.0e+00 |
| B | 3.150890 | 3.150890 | (2.9084; 3.4400) | (2.9084; 3.4400) | 0.1411 | 0.0e+00 |
| Y | 1.771358 | 1.771358 | (1.5141; 2.1998) | (1.5141; 2.1998) | 0.1807 | 0.0e+00 |
| Yalt | 1.734626 | 1.734626 | (1.1903; 2.4890) | (1.1903; 2.4890) | 0.3481 | 0.0e+00 |
| YH1 | 1.786390 | 1.786390 | (1.5036; 2.2594) | (1.5036; 2.2594) | 0.1968 | 0.0e+00 |
| Yshuf | 2.636707 | 2.636707 | (2.3238; 3.0759) | (2.3238; 3.0759) | 0.1884 | 0.0e+00 |
| P | 4.018906 | 4.018906 | (2.9830; 5.4718) | (2.9830; 5.4718) | 0.6192 | 0.0e+00 |
| B0 | 3.096484 | 3.096484 | (2.8449; 3.4035) | (2.8449; 3.4035) | 0.1440 | 0.0e+00 |
| B1 | 3.205471 | 3.205471 | (2.9606; 3.5111) | (2.9606; 3.5111) | 0.1443 | 0.0e+00 |
| G | 3.014042 | 3.014042 | (2.7656; 3.3508) | (2.7656; 3.3508) | 0.1548 | 0.0e+00 |
| U | 2.636419 | 2.636419 | (1.5654; 3.5818) | (1.5654; 3.5818) | 0.5254 | 0.0e+00 |
| A | 0.961373 | 0.961373 | (0.6909; 1.3150) | (0.6909; 1.3150) | 0.1617 | 0.0e+00 |
| AE | 0.508625 | 0.508625 | (0.3592; 0.6517) | (0.3592; 0.6517) | 0.0774 | 0.0e+00 |

Soglia = tetto + δ = 1.336461 (δ = 1) e 0.836461 (δ = 0,5); pavimento − δ = 2.150890 (δ = 1) e 2.650890 (δ = 0,5); pavimento − tetto = 2.814429. Gruppi tumore e stroma: anche lì scarto 0 (tabella `rev_D_compare.csv`).

Confronti appaiati (regola di S3c, δ = 1): B − T +2,814 (w 0,568) «diverso»; Yalt − Y −0,037 (w 0,674) «uguale»; YH1 − Y +0,015 (w 0,143) «uguale»; Yshuf − Y +0,865 (w 0,413) «indeterminato»; Yalt − T +1,398 (w 0,777) e YH1 − T +1,450 (w 0,669) «indeterminati». Tutti uguali a D.json (`rev_D_comparisons.csv`).

Robustezza degli estremi (non richiesta; serve a leggere il caso). Bootstrap esatto: tutti i 5⁵ ricampionamenti per tessuto, combinati fra i due tessuti. Errore Monte Carlo: deviazione standard degli estremi su 20 semi diversi.

| braccio | intervallo (B = 2000) | intervallo esatto | sd Monte Carlo (inferiore / superiore) |
|---|---|---|---|
| T | (-0.263; 0.821) | (-0.258; 0.827) | 0.089 / 0.017 |
| B | (2.908; 3.440) | (2.904; 3.465) | 0.007 / 0.013 |
| Y | (1.514; 2.200) | (1.507; 2.183) | 0.005 / 0.023 |
| Yalt | (1.190; 2.489) | (1.173; 2.522) | 0.016 / 0.025 |
| YH1 | (1.504; 2.259) | (1.498; 2.246) | 0.008 / 0.016 |
| Yshuf | (2.324; 3.076) | (2.325; 3.057) | 0.009 / 0.012 |

Con l'intervallo esatto il caso non cambia. Yalt superiore 2,522 e YH1 superiore 2,246 restano sopra 2,151 e i due inferiori restano sotto 2,151 (caso 4); con δ = 0,5 Yalt superiore 2,522 ≤ 2,651 (caso 1).

## 2. Misure dal basso (compito 2)

Copertura: **10 aree su 10 (5 per tessuto)** per V, T (H1 omologo, `S3b_sonda/runs/T`), B (seme 0 e seme 1), Y, Yalt, YH1, Yshuf: 80 esecuzioni. Il braccio P non è stato ricalcolato.
Ingressi per esecuzione: `bins.parquet`, `X.parquet`, `genes.tsv`, `cells.parquet`, `exp_labels_expanded.npy` (espansione di bin2cell), `rctd_blocks.csv` (uscita di RCTD); riferimento `S3_sonda/ref` (`ref_L_counts.mtx`, `ref_L_cells.csv`, `ref_L_genes.tsv`, `markers_L.csv`).
Metodo, codice mio:
- etichette trasferite: profilo multinomiale per tipo (somme del riferimento + 0,5), massima verosimiglianza sui conteggi del nucleo allargato, nessun tipo sotto 20 UMI;
- M-1: UMI dei quadratini in_reg con etichetta espansa, sugli UMI in_reg;
- M-2: UMI medi per cellula, sulle cellule con almeno metà dei quadratini del nucleo in in_reg;
- M-3: quota media di UMI di marcatori di altri tipi, sulle cellule tipizzate con almeno 5 UMI di marcatori;
- M-4: quota di singoletti fra i blocchi non «reject»;
- M-5: concordanza fra il primo tipo RCTD e il tipo del nucleo, sui blocchi con un nucleo di almeno 4 quadratini;
- M-6: coppie di blocchi vicini (destra, sotto) con lo stesso primo tipo.

Esito: le etichette sono identiche in 153 303 cellule su 153 303. M-1...M-6, n_cells, n_typed, n_nuc_blocks e n_pairs coincidono con le tabelle di S3b, S4 e S6 con scarto 0 in 80 esecuzioni su 80 (`rev_M_compare.csv`). Le due esecuzioni di bin2cell e di RCTD non sono state rifatte: le loro uscite sono la catena stessa.

## 3. Fasce di G1 (compito 3)

Regola del contratto (`docs/PAPER.md`): costruzione se D ≤ soglia; riposizionare se soglia < D ≤ pavimento − δ; chiudere se D > pavimento − δ. Precondizione: pavimento − tetto = 2,814 > δ (sia con 1 sia con 0,5), quindi la prima riga («chiudere» perché pavimento e tetto sono vicini) non scatta. B − T è «diverso» anche con l'intervallo (2,814 − 0,568 > 1). Fascia dell'intervallo = fasce toccate dagli estremi.

| braccio | δ | D (intervallo) | fascia della stima | fascia dell'intervallo |
|---|---|---|---|---|
| Yalt | 1.0 | 1.735 (1.190; 2.489) | riposizionare | a cavallo: costruzione / riposizionare / chiudere |
| Yalt | 0.5 | 1.735 (1.190; 2.489) | riposizionare | riposizionare |
| YH1 | 1.0 | 1.786 (1.504; 2.259) | riposizionare | a cavallo: riposizionare / chiudere |
| YH1 | 0.5 | 1.786 (1.504; 2.259) | riposizionare | riposizionare |
| Y | 1.0 | 1.771 (1.514; 2.200) | riposizionare | a cavallo: riposizionare / chiudere |
| Y | 0.5 | 1.771 (1.514; 2.200) | riposizionare | riposizionare |
| B | 1.0 | 3.151 (2.908; 3.440) | chiudere | chiudere |
| B | 0.5 | 3.151 (2.908; 3.440) | chiudere | chiudere |

RIESAME_C par. 6, con δ = 1 (2,151). Estremo superiore di Yalt 2,489 > 2,151: non è il caso 1. Estremo superiore di YH1 2,259 > 2,151: non è il caso 2. Estremi inferiori 1,190 e 1,504, entrambi < 2,151: non è il caso 3. **Caso 4 (indeterminato; a G1 con i numeri, senza proroga automatica).** Con δ = 0,5 (2,651): estremo superiore di Yalt 2,489 ≤ 2,651, quindi caso 1. La regola dice che la sensibilità con δ = 0,5 non decide. Coincide con `rule_outcome` di D.json.

## 4. A5 e A4 (compito 4)

A5 (`results/sonda/S8/per_unit.csv`, 150 righe = 30 unità × 5 metodi, nessun duplicato). Guadagno = score − score_shuf; B = media dei due semi per area; 100 divisioni (tumore 3 + 2, stroma 2 + 3); tau-b di Kendall; bootstrap 2000 con default_rng(20261201) (una divisione a caso, poi ricampionamento delle aree dentro ciascuna metà e per tessuto).

| quantità | dichiarato | ricalcolato |
|---|---|---|
| T5 (guadagno) | 1,000 ± 0,060 | 1.000000 ± 0.060439 |
| S5c (guadagno) | 0,400 ± 0,031 | 0.400000 ± 0.030899 |
| A5 (guadagno) | −0,600 ± 0,066 | -0.600000 ± 0.065774, esito «non regge» |
| CP-gain tau-b | 0,6 | 0.600000 (UMI del pannello da `S8_sonda/X1/*/cells.parquet`) |
| A5 col coseno (S8) | −0,598 ± 0,067 | -0.598000 ± 0.066585 |

Le 2000 repliche coincidono con `results/sonda/S9/A5_gain_boot.csv` entro 1,1e-16.

A4 (tabelle `results/sonda/S9/{U,C}/area_values_*.csv` e `h_*.csv`, stimatore appaiato mio, stessi indici; e = D(B) − (D(T) + 1), se = f × sd bootstrap di D(B) − D(T)):

| catena | tetto | D(B) | e ± se | esito |
|---|---|---|---|---|
| U (LuCA con RCTD) | 0.328025 | 4.436230 | +3.108204 ± 0.396161 | non regge |
| C (CellTypist) | 0.186117 | 2.293422 | +1.107305 ± 0.485223 | non regge |

Coincidono con `A4_U.json` e `A4_C.json` (scarto 0 alla precisione stampata). Controllo in più: ho ricostruito h di `results/sonda/S3b/h.csv` (bootstrap della media di V con default_rng([20261115, i, j]), i = 1 per il tumore e 2 per lo stroma, j = indice della misura da 0) con scarto massimo 5,7e-14.

## 5. Discrepanze

Nessuna.

## 6. Anomalie viste

1. `results/misure.tsv`: nessuna riga duplicata (12 sessioni, una riga ciascuna). A1 = 3,15089 (se 0,141124) compare tre volte, in S3c, S4 e S6. Non è un duplicato: A1 non cambia per costruzione.
   - La riga di S8 porta A5 = −0,598, cioè il punteggio col coseno, superato in S9.
   - La riga di S9 registra solo A4 (catena U). A5 col guadagno (−0,600) e A2 (9 voci) stanno in `results/misura.json` ma non nel registro.
   - Il file ha terminatori di riga CRLF.
2. `docs/STATO.md` (riga 71) scrive che l'intervallo di D(Yalt) è «a cavallo di 2,151». È vero, ma incompleto: l'estremo inferiore 1,190 sta anche sotto la soglia 1,336, quindi con δ = 1 l'intervallo tocca tutte e tre le fasce di G1. Il contratto non dice come leggere G1 con gli intervalli; RIESAME_C par. 6 legge solo il confine pavimento − δ.
3. T5 sul guadagno vale 1 in tutte le 100 divisioni: la classifica di X1 per guadagno è la stessa in ogni metà, quindi il tetto di A5 è saturo. Il suo se (0,060) viene solo dal ricampionamento dentro le metà.
4. Durante la revisione (2026-11-27, 17:45-17:49 ora italiana) un altro processo, non il mio, ha riscritto `/mnt/micron/geo_spatialtrans/S3_sonda/smoke/` (D.json, per_measure.csv, chain.log). Il suo D.json coincide con `results/sonda/S3c/D.json` (scarto 0). `results/misura.json` indica come fonte di A1 il D.json del check smoke e non `results/sonda/S3c/D.json`. I file che ho letto non sono cambiati (punto 7).
5. Nel repo c'è la cartella non tracciata `results/F3.2/sheets/` (`git status`).
6. Con δ = 0,5 il caso di RIESAME_C passa da 4 a 1: la sensibilità cambia il caso. La regola dice che la sensibilità non decide; lo segnalo perché un lettore del solo D.json vede «case 1».
7. md5: 634 file letti (tabelle, 80 cartelle della catena con 7 file ciascuna, riferimento HLCA, `cells_src.parquet` di S6, `cells.parquet` di X1 in S8), identici prima e dopo (`md5_before.txt` = `md5_after.txt`). Eccezione: le 4 tabelle di S9 U e C le ho registrate solo dopo la lettura. Sono tracciate da git e pulite rispetto a HEAD 17010b3, che era HEAD anche all'inizio.

## Deviazioni dal compito

- L'espansione di bin2cell e le chiamate di RCTD sono uscite della catena usate così come sono, non rieseguite.
- A4 è ricalcolata dalle tabelle per area di S9, non dalle uscite per cella delle catene U e C.
- A5 è ricalcolata dalla tabella per unità di S8, senza rifare i punteggi per cella.
- L'ordine delle chiamate al generatore casuale nei bootstrap (tumore poi stroma; divisione, poi metà 1, poi metà 2) l'ho preso dal codice del progetto: serviva a riprodurre gli intervalli con lo stesso schema di semi. Il calcolo è mio. Il bootstrap esatto e i 20 semi in più non dipendono da quell'ordine.
- Per D ho usato h di `h.csv`, `h_U.csv` e `h_C.csv`; ho ricostruito in modo indipendente solo h di S3b.
