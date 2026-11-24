---
type: State
title: "Stato di geo_spatialtrans"
description: "Memoria fra le sessioni: fase, dati, numeri, bilancio dell'errore, minacce, strade fallite, prossima mossa."
status: draft
generated: { by: geo-agent/progetto-master, at: 2026-11-24T11:45:00+01:00 }
---
# Fase

Sonda, budget 7. Registro: S1, S2, S5 (`infra`), S3b, S3c, S4 (`misura`). Usate 6, resta 1. Luca 2026-11-24: «Oggi solo A3; proroga della sonda da 7 a 9» (A5, poi A4 con A2, poi G1). Seconda proroga: lo strumento chiede prima un riesame (`chiudi --riesame`).

# Dati (polmone 10x, CG000709, 1 blocco FFPE, 1 donatore)

- H1 (T), H2 (V) in `/mnt/micron/geo_spatialtrans/data_lung/`; md5 in `results/sonda/S3/MD5_download_H.txt`. HLCA core in `ref_lung/`.
- Metà tenuta da parte: solo H2 (x < m − 50 µm, m = 3372), X1 (m = 3300), X2. H1 letta intera.
- Orientamento (S5, `results/sonda/S5/orient.json`, `local_ecc.json`): X1 come H2 (rotazione −0,03°), H1 ruotata di 180°, X2 di circa 5,7°. Le 10 aree di V riportate in X1 stanno nella metà di sviluppo di X1. Residuo X1 -> H2 (ECC): mediana 8 µm, 95° percentile 19 µm.

# Metro (S3c, `docs/sonda/REGOLA_S3c.md`)

Stimatore appaiato `tools/sonda/S3c_D.py`, bootstrap per coppie stratificato, B = 2000, seme 20261123, intervallo (ii) con f = 1,584. Cancelli G-ID, C-1p, C-2, C-4, C-5, CP-perm passati.

# Numeri (h)

Pavimento D(B) = 3,151 (2,908; 3,440), tetto D(T) = 0,336, soglia 1,336. A1 oggi 3,151 ± 0,141: non regge.

Componenti veri (S4, `results/sonda/S4/D.json`): D(Y, mappa dei tipi) 1,771; D(G, nuclei) 3,014; D(U, UMI) 2,636; D(A) 0,961; D(AE) 0,509. Quote del divario chiuse: 0,69; 0,09; 0,30; 0,92; 0,99.

# S5, A3 (`docs/sonda/REGOLA_S5.md`, commit a9f782d e 8f8c0a4)

- Braccio XB costruito (Xenium X1 in quadratini, nuclei Xenium, correzione affine per area): G-XR 10/10 (scacchiera ≤ 1,81 µm, dopo la correzione ≤ 0,41 µm), G-SUB 50/50, G-HALF, G-MD5 passati. Anche XV (Visium HD di X1) esportato.
- A3 non misurabile come scritta: sui 284 geni del pannello H2 ha il 2 % degli UMI, 1-7 UMI per blocco da 8 µm (mediana); blocchi con ≥ 20 UMI 0-4 % in V, 0 % in B. RCTD (UMI_min 20) non parte in B 20/20, P 9/10, V 3/10 (`results/sonda/S5/diag_summary.csv`; ricalcolo indipendente concorde).
- P1 smentita di poco: f tumore 0,094, stroma 0,113 (Xenium circa 10 volte più profondo di H2 sui geni del pannello).
- Decisione aperta per Luca (al riesame): A3 su X2 (5K, serve registrarla), oppure A3 tolta e il fatto dichiarato. A5 non dipende da questo (Visium HD di X1 su tutti i geni, Xenium come verità).

# Check

`check smoke` = `tools/sonda/S3c_smoke.sh` (circa 5 min). Non contiene S4 né S5. `check full` da scrivere.

# Bilancio dell'errore (per componente, S4)

1. Disposizione dei tipi: 69 % del divario.
2. UMI concentrati sui nuclei (M-1): torna solo con nuclei e UMI veri insieme (A); limite superiore.
3. Espressione della singola cellula: 7 %. Nessuna sessione.
4. Geometria dei nuclei: 9 %. Nessuna sessione.

Con i componenti veri D = 0,509 (AE) e 0,961 (A), sotto la soglia.

# Minacce alla validità

1. HLCA è tessuto sano: le cellule tumorali non hanno un tipo (M-4).
2. Y porta anche la composizione locale dell'area: un generatore impara le statistiche, non la mappa.
3. Un donatore: tetto fra sezioni seriali.
4. A3 senza misura: il confronto con lo Xenium in quadratini, l'alternativa più semplice, manca.
5. U chiude il 30 % del divario senza tipi veri (intervallo 2-60 %): non spiegato.

# Strade fallite

- Generatore dello strato cellulare (colon): 7,83 contro 7,77 del banale.
- Tetto di S3 su aree non omologhe e κ del colon: C-4 rosso.
- Stimatore indipendente su aree omologhe (S3b): tetto −0,77.
- A3 sul pannello di X1 alla profondità di H2: troppo pochi UMI per RCTD.
- Correlazione di fase su tessere da 1 mm per X1: 240 µm, falsa; usare l'ECC.

# Prossima mossa

1. Riesame (nuova chat, senza codice): fatti di S1-S5, decisione su A3, proroga a 9.
2. A5 sullo Xenium di X1 (tetto da riscrivere senza X2). 3. A4 con A2. 4. G1.

# Aperti

- `results/F3.2/sheets/` non tracciato: lasciato com'è.
- `tools/okf.py check` rosso da SC1.
- Revisione indipendente dei numeri di S1, S2, S3b non fatta (prima di G1).
