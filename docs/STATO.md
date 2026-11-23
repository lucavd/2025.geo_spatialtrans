---
type: State
title: "Stato di geo_spatialtrans"
description: "Memoria fra le sessioni: fase, dati, numeri, bilancio dell'errore, minacce, strade fallite, prossima mossa."
status: draft
generated: { by: geo-agent/progetto-master, at: 2026-11-23T22:00:00+01:00 }
---
# Fase

Sonda, budget 5, firmata a G0 il 2026-11-20. Registro: S1, S2 (`infra`), S3b, S3c (`misura`). Usate 4 sessioni; budget prorogato da 5 a 7 il 2026-11-23: ne restano 3. S3 non chiusa (C-4 rosso), sostituita da S3b.

# Dati (polmone 10x, CG000709, 1 blocco FFPE, 1 donatore)

- H1 (`..._Experiment1`, T) e H2 (`..._Experiment2`, V) in `/mnt/micron/geo_spatialtrans/data_lung/H1`, `/H2`; md5 in `results/sonda/S3/MD5_download_H.txt`. HLCA core in `/mnt/micron/geo_spatialtrans/ref_lung/`.
- MALAT1, NEAT1, KRT18 non sono nel pannello di sonde.
- Metà tenuta da parte (Luca, 2026-11-23): solo H2 (x < m − 50 µm, m = 3372 µm), X1, X2. H1 letta intera, mai in taratura (griglie di H1 e H2 ruotate di 180°).

# Metro (S3c, `docs/sonda/REGOLA_S3c.md`, commit a93641d prima del calcolo)

- Tabella e h di S3b (`results/sonda/S3b/area_values.csv`, `h.csv`; md5 invariati, `results/sonda/S3c/md5_*.txt`). Aree di T = immagini registrate in H1 delle aree di V (S3b).
- Stimatore appaiato (`tools/sonda/S3c_D.py`): per cella d_k = braccio − V sull'area k, q = (media(d)² − var(d)/5)/h², D = radice con segno della media di q sulle 12 celle; tutti i bracci (T, B, P, B0, B1).
- Incertezza: bootstrap per coppie stratificato per tessuto, B = 2000, seme 20261123, semi-ampiezze ed se × f = 1,584 (√(5/4) t(0,975; 4)/1,96). Scelta da COV: copertura minima (i) 0,807, (ii) 0,938.
- Cancelli passati: G-ID (permutazioni = stimatore indipendente, 1,7e-14), C-1p (0,068; 0,999; 1,996), C-2, C-4, C-5; CP-perm rango 1/120 in tumore e stroma.
- Diagnosi se (`tools/sonda/S3c_diag_se.py`, `results/sonda/S3c/diag_se.json`): sui sintetici con la struttura osservata la deviazione standard vera di D(B) è 0,11; il metro parametrico di S3b stimava 0,234 (circa il doppio), (ii) 0,164 (prudente). Il se di S3b era gonfiato.
- Revisione indipendente (sotto-agente, seme 777): 126 voci su 127 concordano; unica differenza f con 1,959964 invece di 1,96 (1,8e-5).

# Numeri (h, `results/sonda/S3c/D.json`)

Tumore + stroma: pavimento D(B) = 3,151 (2,908; 3,440), se 0,141; tetto D(T) = 0,336 (−0,263; 0,821); soglia 1,336 (δ = 1; 0,836 con δ = 0,5; 2,336 con δ = 2); D(P) = 4,019. B − T 2,81 (w 0,57) diverso; P − T 3,68 (w 1,29) diverso; B0 − B1 −0,11 (w 0,12) uguale. Pavimento − tetto 2,8 h > δ. In `progetto.json`: pavimento 3,151, tetto 0,336, soglia 1,336. A1 oggi 3,151 ± 0,141: non regge (è il banale).

Osservazioni (non decidono nulla): D(B1) > D(B0) in entrambi i tessuti (0,11 h, entro δ). La correlazione area per area fra B e V è bassa (mediana 0,11; T-V 0,94): il banale non riproduce le differenze fra aree.

# Check

`check smoke` = `tools/sonda/S3c_smoke.sh` (circa 6 min): S1_test_cc, catena ridotta, C-3, metro S3c con i suoi cancelli, identità con `results/sonda/S3c/D.json`, G-ID, CP-perm, metodo d'intervallo. `S3_smoke.sh` ha l'interruttore `S3_DSCRIPT` (predefinito S3_D.py). `check full` da scrivere.

# Bilancio dell'errore (`results/sonda/S3c/per_measure.csv`: quota di q_B − q_T, T appaiato)

Tumore M-4 (singoletti RCTD) 38 %; stroma M-4 22 %; tumore M-1 (UMI assegnati) 12 %; stroma M-6 9 %; tumore M-6 7 %; stroma M-1 7 %; stroma M-5 5 %; il resto sotto l'1 %. Singoletti RCTD: 60 %.

# Minacce alla validità

1. HLCA è tessuto sano: le cellule tumorali non hanno un tipo (pesa su M-4, la voce più grande).
2. Registrazione H&E -> quadratini di Space Ranger non misurata su H1 e H2.
3. Un donatore: tetto fra sezioni seriali; h è la precisione della media di H2, non un margine fra pazienti.
4. A5 (X1 contro X2) può avere lo stesso problema delle metà opposte: controllare l'orientamento delle griglie prima di usarle.
5. Il se scelto è prudente (circa 1,5 volte il vero per D(B)): rende più difficile «avanzare» in costruzione.

# Strade fallite

- Generatore dello strato cellulare (colon): 7,83 contro 7,77 del banale.
- Tetto di S3 su aree non omologhe e κ del colon: C-4 rosso.
- Stimatore indipendente su aree omologhe (S3b): tetto −0,77, distorto.

# Prossima mossa

Proroga della sonda firmata (Luca, 2026-11-23: «va bene la proroga e chiudi»): budget da 5 a 7, restano 3 sessioni. Ordine proposto: (1) componenti veri sulle voci più grandi del bilancio, singoletti RCTD (60 %) e UMI assegnati, per sostituzione con la versione perfetta o spenta; (2) A3 e A5 sullo Xenium, prima l'orientamento delle griglie di X1 e X2; (3) G1. Atteso della proroga: D ≤ pavimento − δ = 2,15 h; altrimenti G1 con i numeri raggiunti (chiudere o riposizionare).

# Aperti

- `results/F3.2/sheets/` non tracciato: lasciato com'è.
- `tools/okf.py check` rosso da SC1 (esenzioni non modificate).
- Revisione indipendente dei numeri di S1, S2, S3b non fatta (prima di G1); quelli di S3c sì.
- X2 senza registrazione misurata.
