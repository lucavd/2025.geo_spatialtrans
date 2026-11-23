---
type: Rule
title: "S3b, metro riparato: h senza i κ del colon, tetto su aree omologhe di H1, metà tenuta da parte solo in H2"
description: "Regola scritta e committata prima di leggere l'espressione di H1 fuori dalle aree di S3: registrazione H1 -> H2, aree omologhe, h con κ = 1, controlli del metro, controprova, attese."
status: stable
generated: { by: geo-agent/progetto-master, at: 2026-11-23T20:30:00+01:00 }
sources:
  - id: s3
    path: docs/sonda/REGOLA_S3.md
  - id: paper
    path: docs/PAPER.md
  - id: expl
    path: results/sonda/S3b/overlap.json
---
# Perché

S3 (`REGOLA_S3.md`[^s3]) ha lasciato il controllo positivo C-4 rosso. Decisione di Luca, 2026-11-23 («1»): h senza i κ del colon e aree omologhe fra H1 e H2 per posizione registrata; C-4 resta cancello. Questa regola sostituisce in S3 solo le parti scritte qui; il resto di REGOLA_S3 vale invariato (catena, riferimento HLCA, braccio banale B, controllo positivo P, M-4, gruppi, δ = 1,0 h, bootstrap e semi).

# Già noto prima del commit (esplorazione senza espressione, commit 506f0df)

- Le griglie di array di H1 e H2 sono ruotate di 180° l'una rispetto all'altra nell'H&E; la regola della metà (x < m − 50 µm, per sezione) prende metà opposte del tessuto: in comune 0,02 mm² su 26,4 mm² (`results/sonda/S3b/overlap.json`[^expl]). Tutte le 20 aree di S3 cadono nella metà tenuta da parte dell'altra sezione.
- Registrazione globale H1 -> H2 (ECC affine di OpenCV sull'H&E hires, `tools/sonda/S3b_explore_overlap.py`): Dice delle maschere 0,98, rotazione circa 1,6° oltre ai 180°; residuo locale su 75 tessere da 1 mm (ECC di sola traslazione, `S3b_explore_local.py`): mediana 15 µm, 95° percentile 37 µm, massimo 42 µm; controllo positivo 40 µm ritrovato a 39,5 µm.
- Dalla diagnosi di S3 (non decide nulla, `/tmp`): con κ = 1 sulle aree di S3, D(T) 0,85, D(B) 3,14, P − T 3,04 (w 1,77). V, B e P non cambiano in S3b: per loro questi numeri sono già visti; il solo braccio nuovo è T.

# Metà tenuta da parte (decisione di Luca, 2026-11-23: «Sì, nascondi metà solo di H2»)

- H2: metà di REGOLA_S1, invariata (x < m − 50 µm, m = 3372 µm); il controllo nel codice resta.
- H1: letta per intero; solo riferimento, nessuna taratura usa H1. Il controllo della metà si spegne solo per H1 (variabile `S3_NOGUARD=H1`).
- Conferma: le aree di V' (metà tenuta da parte di H2) e le loro immagini in H1, che cadono nella metà di H1 già letta in S3 (10 aree, nessuna taratura): si dichiara.

# Aree di T (omologhe)

- Le 10 aree di V sono quelle di S3 (`results/sonda/S3/areas_H2.csv`), invariate; quindi V, B e P non si rifanno.
- Per ogni area di V: centro (x0 + lato/2, y0 + lato/2) in µm del riferimento di array di H2 -> pixel hires di H2 (affine array -> pixel per minimi quadrati, come nell'esplorazione) -> H1 con la registrazione globale di `overlap.json` (invariata) -> riferimento di array di H1. Area di T = quadrato di lato 750 µm allineato agli assi di H1, centrato nel punto riportato, angolo arrotondato a 2 µm. Nome dell'unità uguale a S3 (`H1_<tessuto>_<k>`, con k = numero dell'area di V), radice `/mnt/micron/geo_spatialtrans/S3b_sonda`: le uscite di S3 non si toccano.
- Uscita: `results/sonda/S3b/areas_H1_hom.csv`, con la frazione del quadrato vero riportato (ruotato) coperta dal quadrato allineato e la copertura in_tissue in H1.

# h (decisione di Luca: «κ = 1: sola precisione della media di H2»)

Costruzione di `tools/sonda/S3_h.py` con κ = 1 e fattore di asimmetria = 1 (anch'esso tarato sul colon, R7): h = (min(1, hi) − max(0, lo)) / 2 per le quote, (hi − max(0, lo)) / 2 per M-2, con lo e hi i percentili 2,5 e 97,5 del bootstrap della media sulle 5 aree di V (B = 2000, `default_rng([20261115, i, j])`, come S3). M-4 almeno 0,01. Uscita `results/sonda/S3b/h.csv`. Conseguenza accettata da Luca: la soglia è un margine di precisione della media di H2, non un margine fra pazienti.

# Cancelli (tecnici; se uno fallisce D non si legge)

| id | controllo | esito richiesto |
|---|---|---|
| G-R | registrazione H1 -> H2: 95° percentile del residuo locale; controllo positivo di 40 µm | ≤ 75 µm (un decimo del lato); 40 ± 5 µm |
| G-O | omologia geometrica: frazione del quadrato vero riportato coperta dal quadrato di T, per area | ≥ 0,95 per tutte le 10 |
| G-C | copertura in_tissue delle aree di T | ≥ 0,8 per tutte; sotto 0,9 si segnala |
| C-1, C-2, C-3, C-5 | come REGOLA_S3, sulla tabella nuova | come REGOLA_S3 |
| C-4 | controllo positivo P come REGOLA_S3, con h e T nuovi | D(P) − D(T) «diverso», P peggio |

# Controprova

- CP-1 (l'omologia conta): con lo stesso h nuovo, D(T omologo) < D(T di S3). Se no, la minaccia 2 di STATO non era la causa: si scrive.
- CP-2 (le aree sono davvero lo stesso tessuto): punteggio epitelio/stroma r di REGOLA_S3 (geni E e S, matrice filtrata da 8 µm) per le 10 aree di T omologhe contro le 10 di V: correlazione di Spearman ≥ 0,7 e segno di r concorde per tessuto in almeno 8 aree su 10. Le aree di T di S3 contro V come confronto (attesa: più bassa).

# Attese (scritte prima)

- D(T omologo) fra 0 e 0,85 h (sotto la diagnosi di S3, 0,85).
- D(B) circa 3,1 h (noto), quindi pavimento − tetto > 1 h: G1 non dice «chiudere» per la prima condizione.
- C-4 passa (P − T circa 3, w circa 1,8).

# Calcolo

8 processi a un thread (sul server gira un altro lavoro di circa 16 processi), file `WORKERS` della radice di S3b. Catena di REGOLA_S3 su T: ritaglio H&E e StarDist, esportazione, bin2cell, RCTD, misure. Tabella = righe V, B, P di `results/sonda/S3/area_values.csv` più T nuovo (`results/sonda/S3b/area_values.csv`); D con `tools/sonda/S3_D.py` (uscite in `results/sonda/S3b/`).

# Esito e conseguenze

- Tutti i cancelli passano: pavimento = D(B), tetto = D(T), soglia = tetto + 1 h in `progetto.json`; A1 = D(B) in `results/misura.json`; il check smoke usa il metro nuovo (identità con `results/sonda/S3b/D.json`); S3b si chiude come misura.
- Un cancello fallisce: numeri provvisori, fuori da `progetto.json`; la sessione lo dichiara.

[^s3]: docs/sonda/REGOLA_S3.md
[^expl]: results/sonda/S3b/overlap.json
