---
type: State
title: "Stato di geo_spatialtrans"
description: "Memoria fra le sessioni: fase, dati, numeri, bilancio dell'errore, minacce, strade fallite, prossima mossa."
status: draft
generated: { by: geo-agent/progetto-master, at: 2026-11-23T16:00:00+01:00 }
---
# Fase

Sonda, budget 5, firmata a G0 il 2026-11-20. Nel registro S1 e S2 (`infra`). S3 aperta nella stessa chat di S2, in deroga (Luca, 2026-11-21: «ok a sessione 3 consigliata in deroga»). S3 non chiusa: `check smoke` rosso per il controllo positivo C-4 (sotto), quindi `chiudi` rifiuta anche `--infra` (2026-11-23: «l'ultimo check e' rosso; manca la misura di A1»). Allarme «più di una sessione su tre senza misura» aperto (con S1 e S2 `infra` resterebbe acceso anche dopo una sessione con misura).

# Dati (polmone 10x, CG000709, 1 blocco FFPE, 1 donatore)

- H1 (`..._HD_Only_Experiment1`, T) e H2 (`..._Experiment2`, V) in `/mnt/micron/geo_spatialtrans/data_lung/H1`, `/H2`: 8 file, dimensioni uguali al CDN, md5 in `results/sonda/S3/MD5_download_H.txt`, ricontrollo in `md5_recheck_H.txt`. HLCA core in `/mnt/micron/geo_spatialtrans/ref_lung/` (md5 `MD5_HLCA.txt`).
- MALAT1, NEAT1 e KRT18 non sono nel pannello di sonde.

# S3: metro D e catena banale (`docs/sonda/REGOLA_S3.md`, commit 31b773c, prima di ogni lettura di espressione)

- Codice in `tools/sonda/S3_*`; uscite in `/mnt/micron/geo_spatialtrans/S3_sonda/` (una cartella per unità); tabelle in `results/sonda/S3/`.
- Aree: rilassamento 2 (lato 750 µm, copertura ≥ 0,9) in H1 e H2; in H2 le aree di tumore sono poco epiteliali (r da 0,59 a −0,02, tumour 5 segnalata), in H1 nettamente (r da 1,09 a 2,54): le due metà di sviluppo coprono tessuto diverso.
- Riferimento: HLCA core `ann_level_2`, 11 tipi, 27 066 cellule, 17 015 geni (`S3_sonda/ref/ref_L.json`).
- B tarato su H2: densità e profondità realizzate sui bersagli; nucleo realizzato più piccolo (9,1 contro 9,73 quadratini nel tumore, 6,1 contro 6,45 nello stroma).
- Corsa: 50 unità, 20 processi, 1 h 55 min; Tctl massimo 83,6 °C.

# Numeri (h, `results/sonda/S3/D.json`), provvisori

Tumore + stroma: D(T) = 0,521 (0,052; 1,055); D(B) = 1,573 (1,212; 1,955), se 0,187; D(P) = 1,874. B − T = 1,052, w 0,526: indeterminato. Soglia che ne verrebbe: 1,521. Attese scritte: tetto circa 2, pavimento circa 7-8: sbagliate di un fattore 4-5.

Controlli: C-1 passa (−0,056; 0,998; 2,012); C-2 passa (scarto 0,0 sul colon); C-3 passa (catena ridotta rifatta da capo, scarto 0); C-5 passa (B0 − B1 = −0,069, w 0,116, uguale); C-4 NON passa (P − T = 1,353, w 1,004, indeterminato). Per la regola i numeri sono provvisori: pavimento, tetto e soglia non sono in `progetto.json`; A1 sta in `results/sonda/S3/misura_provvisoria.json`, non in `results/misura.json`.

Diagnosi (non decide nulla, `/tmp`): con h senza i κ del colon (κ = 1) C-4 passerebbe (P − T = 3,04, w 1,77) e B − T sarebbe diverso (2,29, w 0,87; D(T) 0,85, D(B) 3,14). Il fallimento viene dalla scelta scritta di usare i κ del colon (tumore 1,66, stroma 3,19), che allarga h; con quella scelta w(P, T) ≥ δ, quindi «diverso» è quasi impossibile per costruzione.

# Check

`check smoke` = `tools/sonda/S3_smoke.sh` (circa 5 minuti): S1_test_cc, catena ridotta da capo a fondo (250 µm di H2 tumour 1 e un B), C-3, metro con C-1, C-2, C-4, C-5. Rosso per C-4. `check full` da scrivere.

# Bilancio dell'errore (provvisorio, `results/sonda/S3/per_measure.csv`, termine di B meno termine di T in h²)

Tumore singoletti RCTD 62 %; tumore UMI assegnati 20 %; tumore vicini dello stesso tipo 12 %; stroma singoletti 10 %; il resto sotto il 5 %. Tumore marcatori estranei e concordanza col nucleo: T più lontano di B.

# Minacce alla validità

1. Metro senza controllo positivo superato: h dipende da κ trasferiti dal colon.
2. Aree di tumore non omologhe fra H1 e H2 (metà di sviluppo su tessuto diverso).
3. HLCA è tessuto sano: le cellule tumorali non hanno un tipo.
4. Registrazione H&E -> quadratini di Space Ranger non misurata su H1 e H2; il metro non vede uno spostamento dei nuclei (le etichette trasferite vengono dagli stessi quadratini).
5. Un donatore: tetto fra sezioni seriali.

# Strade fallite

- Generatore dello strato cellulare (colon): 7,83 contro 7,77 del banale.

# Prossima mossa

Decisione di Luca, 2026-11-23 («1»): riparare il metro con una regola nuova scritta prima (S3b, chat nuova). h del polmone senza κ del colon (κ = 1: sola precisione della media di H2), dichiarando la diagnosi κ = 1 già vista; aree omologhe fra H1 e H2 per posizione registrata (nuova corsa, circa 2 h); C-4 resta cancello. Poi check verde e chiusura come misura. Accettato: la soglia diventa un margine di precisione, non un margine fra pazienti.

# Aperti

- `results/F3.2/sheets/` non tracciato: lasciato com'è.
- `tools/okf.py check` rosso da SC1 (esenzioni non modificate).
- Revisione indipendente dei numeri di S1, S2, S3 non fatta (prima di G1).
- X2 senza registrazione misurata.
- Idea: aree omologhe fra H1 e H2 per posizione registrata invece che per punteggio.
