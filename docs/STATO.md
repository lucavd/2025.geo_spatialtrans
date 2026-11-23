---
type: State
title: "Stato di geo_spatialtrans"
description: "Memoria fra le sessioni: fase, dati, numeri, bilancio dell'errore, minacce, strade fallite, prossima mossa."
status: draft
generated: { by: geo-agent/progetto-master, at: 2026-11-23T22:00:00+01:00 }
---
# Fase

Sonda, budget 5, firmata a G0 il 2026-11-20. Registro: S1, S2 (`infra`), S3b (`misura`, A1 = 3,14). Allarme «più di una sessione su tre senza misura» ancora acceso (2 su 3). S3 non chiusa (C-4 rosso), sostituita da S3b.

# Dati (polmone 10x, CG000709, 1 blocco FFPE, 1 donatore)

- H1 (`..._Experiment1`, T) e H2 (`..._Experiment2`, V) in `/mnt/micron/geo_spatialtrans/data_lung/H1`, `/H2`; md5 in `results/sonda/S3/MD5_download_H.txt`. HLCA core in `/mnt/micron/geo_spatialtrans/ref_lung/`.
- MALAT1, NEAT1, KRT18 non sono nel pannello di sonde.
- Metà tenuta da parte (Luca, 2026-11-23, S3b): solo H2 (x < m − 50 µm, m = 3372 µm), X1, X2 come REGOLA_S1; H1 letta intera, mai in taratura. Motivo: le griglie di H1 e H2 sono ruotate di 180° nell'H&E, la regola per sezione prendeva metà opposte del tessuto (0,02 mm² in comune su 26,4).

# S3b: metro riparato (`docs/sonda/REGOLA_S3b.md`, commit 21bdd0a prima della corsa)

- h con κ = 1 e fattore di asimmetria 1 (`tools/sonda/S3b_h.py`): sola precisione della media delle 5 aree di V.
- T su aree omologhe: aree di V di S3 riportate in H1 (registrazione H&E hires, ECC affine; residuo locale p95 37 µm, controllo positivo 39,5 su 40). G-O ≥ 0,988, copertura ≥ 0,90. CP-2: punteggio epitelio/stroma omologo contro V, Spearman 0,90, segno concorde 8/10.
- V, B, P invariati da S3; T rifatto (10 aree, 8 processi, 48 min). Uscite `/mnt/micron/geo_spatialtrans/S3b_sonda/`, tabelle `results/sonda/S3b/`.
- Codice: `tools/sonda/S3b_*`; interruttori in S3_common/S3_export/S3_D/S3_smoke (S3_AREAS_H1, S3_NOGUARD, --provv, S3_TABLE, S3_H, S3_DREF), predefiniti invariati (C-3 identico).

# Numeri (h, `results/sonda/S3b/D.json`)

Tumore + stroma: D(B) = 3,139 (2,683; 3,636), se 0,241 = pavimento = A1 oggi; D(P) = 3,889; D(T) = −0,769 (−1,988; 0,014). Controlli: C-1 (0,05; 0,99; 2,04), C-2 (0), C-3 (0), C-4 passa (P − T 4,66, w 1,81), C-5 passa (−0,11, w 0,17). Check smoke verde.

Tetto e soglia NON in `progetto.json` (deviazione dichiarata dalla regola): lo stimatore (`tools/riesame_D.py`) sottrae il rumore come se T e V fossero aree indipendenti; le aree omologhe sono appaiate (correlazione mediana 0,94 sulle 12 voci, `results/sonda/S3b/diag_paired.csv`). Diagnosi: D(T) appaiato +0,336 (dentro l'attesa 0-0,85). B non appaiato (r mediana 0,11): pavimento 3,15 in entrambi i modi; P appaiato 4,02. C-1 simula insiemi indipendenti e non vede il difetto. Errore di pre-registrazione: il confronto di CP-2 con le aree di S3 è banale (aree numerate per punteggio in entrambe le sezioni, Spearman 1).

# Check

`check smoke` = `tools/sonda/S3b_smoke.sh` (circa 5 min): S1_test_cc, catena ridotta, C-3, metro S3b con C-1, C-2, C-4, C-5, identità con `results/sonda/S3b/D.json`. `check full` da scrivere.

# Bilancio dell'errore (provvisorio, `results/sonda/S3b/per_measure.csv`: quota dello scarto B − T in h²; T con la formula indipendente)

Tumore M-4 (singoletti RCTD) 36 %; stroma M-4 21 %; tumore M-1 (UMI assegnati) 12 %; stroma M-6 9 %; tumore M-6 7 %; stroma M-1 7 %; stroma M-5 5 %; il resto sotto il 2 %. Singoletti RCTD: 57 % in tutto.

# Minacce alla validità

1. Tetto non ancora stimato correttamente: serve lo stimatore appaiato.
2. HLCA è tessuto sano: le cellule tumorali non hanno un tipo.
3. Registrazione H&E -> quadratini di Space Ranger non misurata su H1 e H2.
4. Un donatore: tetto fra sezioni seriali.
5. A5 (X1 contro X2) può avere lo stesso problema delle metà opposte: controllare l'orientamento delle griglie di X1 e X2 prima di usarle.

# Strade fallite

- Generatore dello strato cellulare (colon): 7,83 contro 7,77 del banale.
- Tetto di S3 su aree non omologhe e κ del colon: C-4 rosso.

# Prossima mossa

S3c (chat nuova): stimatore appaiato per T e P (rumore = var(T_k − V_k)/n, bootstrap per coppie), regola scritta prima, con un controllo sintetico appaiato (C-1 con la correlazione osservata); poi tetto e soglia in `progetto.json`. Nessuna corsa nuova di catena: tabelle di S3b. Budget della sonda: 3 sessioni usate su 5 (S1, S2, S3b); dopo S3c ne resta 1 per componenti veri, A3 e A5: probabile richiesta di proroga, da portare a Luca con i numeri.

# Aperti

- `results/F3.2/sheets/` non tracciato: lasciato com'è.
- `tools/okf.py check` rosso da SC1 (esenzioni non modificate).
- Revisione indipendente dei numeri di S1, S2, S3b non fatta (prima di G1).
- X2 senza registrazione misurata.
