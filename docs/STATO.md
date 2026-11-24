---
type: State
title: "Stato di geo_spatialtrans"
description: "Memoria fra le sessioni: fase, dati, numeri, bilancio dell'errore, minacce, strade fallite, prossima mossa."
status: draft
generated: { by: geo-agent/progetto-master, at: 2026-11-24T13:45:00+01:00 }
---
# Fase

Sonda, budget 7. Registro: S1, S2, S5 (`infra`), S3b, S3c, S4 (`misura`), RC (`riesame`, non conta). Usate 6, resta 1. Proroga 7 -> 9 chiesta da Luca il 2026-11-24; il riesame RC consiglia di legarla all'esito della sessione 7 (`docs/riesame/RIESAME_C.md` §6).

# Dati (polmone 10x, CG000709, 1 blocco FFPE, 1 donatore)

- H1 (T), H2 (V) in `/mnt/micron/geo_spatialtrans/data_lung/`; md5 in `results/sonda/S3/MD5_download_H.txt`. HLCA core in `ref_lung/`.
- Metà tenuta da parte: solo H2 (x < m − 50 µm, m = 3372), X1 (m = 3300), X2. H1 letta intera.
- Orientamento (S5, `results/sonda/S5/orient.json`): X1 come H2, H1 ruotata di 180°, X2 di circa 5,7°. Residuo X1 -> H2 (ECC): mediana 8 µm, p95 19 µm.

# Metro (S3c, `docs/sonda/REGOLA_S3c.md`)

Stimatore appaiato `tools/sonda/S3c_D.py`, bootstrap per coppie stratificato, B = 2000, seme 20261123, metodo (ii) f = 1,584. Ricalcolo indipendente in RC: 30 D e 15 quote con scarto ≤ 4e-15; estremi superiori di B e A bassi di 0,03 h rispetto al bootstrap esatto (seme, trascurabile).

# Numeri (h)

Pavimento D(B) 3,151 (2,908; 3,440), tetto D(T) 0,336 (−0,263; 0,821), soglia 1,336 (δ = 1). A1 = 3,151: non regge.

Componenti veri (S4, `results/sonda/S4/D.json`), tutti copiati dalle aree misurate di H2: Y (tipi) 1,771; G (nuclei) 3,014; U (UMI) 2,636; A 0,961 (estremo superiore esatto 1,346, a cavallo della soglia); AE 0,509. Quote chiuse 0,69; 0,09; 0,30; 0,92; 0,99. Voci: M-4 60% di D² in B; M-1 63% in Y.

# A3 (S5)

Non misurabile su X1 alla profondità di H2: pannello 2% degli UMI di H2, RCTD (UMI_min 20) non parte in B 20/20 (`results/sonda/S5/diag_summary.csv`). Braccio XB costruito (G-XR 10/10). RC consiglia di toglierla e dichiarare il fatto; decide Luca.

# Check

`check smoke` = `tools/sonda/S3c_smoke.sh` (circa 5 min). Non contiene S4 né S5. `check full` da scrivere.

# Bilancio dell'errore (per componente, S4; diagnosi, non fattibilità)

1. Disposizione dei tipi: 69% del divario (copiata dal bersaglio).
2. UMI concentrati sui nuclei (M-1): torna solo con nuclei e UMI veri insieme.
3. Espressione: 7%. Nessuna sessione.
4. Geometria dei nuclei: 9%. Nessuna sessione.

# Minacce alla validità

1. HLCA è tessuto sano: le cellule tumorali non hanno un tipo (M-4).
2. I bracci sotto la soglia copiano il bersaglio: G1 letta su di loro è circolare (RC §1.2).
3. Y copia le etichette trasferite dalla stessa catena (RCTD con HLCA) su V: M-4 e M-6 in parte costruiti; A4 con una tipizzazione diversa (RC §1.3).
4. δ senza fonte valida dopo S3b; con δ = 2 G1 si contraddice (RC §1.4-1.5).
5. Un donatore: tetto fra sezioni seriali.

# Strade fallite

- Generatore dello strato cellulare (colon): 7,83 contro 7,77 del banale.
- Tetto di S3 su aree non omologhe e κ del colon: C-4 rosso.
- Stimatore indipendente su aree omologhe (S3b): tetto −0,77.
- A3 sul pannello di X1 alla profondità di H2: troppo pochi UMI per RCTD.
- X2: registrazione non giudicabile (S1, S2): un colpo.
- Correlazione di fase su tessere da 1 mm per X1: falsa; usare l'ECC.

# Riesame RC (2026-11-24)

Fatti `docs/riesame/FATTI_C.md`, giudizio `docs/riesame/RIESAME_C.md`. Nessuna condizione di chiusura vale per D (pavimento − tetto 2,81 h > δ). Verdetto: continuare, messaggio invariato, condizionato alla prova dei tipi non copiati. Decisioni aperte per Luca: (1) G1 solo su bracci non copiati e proroga legata all'esito; (2) δ = 1 h con fonte nuova (precisione dell'area vera), sensibilità solo 0,5; (3) A3 tolta.

# Prossima mossa

1. Sessione 7 (sonda): bracci Y-alt (tipi dall'area successiva dello stesso tessuto, permutazione ciclica fra le 5 aree di V) e Y-H1 (tipi dall'area omologa di H1). Regola di esito già scritta in RIESAME_C §6, soglia 2,151 h.
2. Se passa: proroga 7 -> 9, poi A5 su X1 (tetto da riscrivere senza X2), A4 con A2 e tipizzazione diversa, G1.
3. Se non passa: G1 senza proroga.

# Aperti

- `results/F3.2/sheets/` non tracciato: lasciato com'è.
- `tools/okf.py check` rosso da SC1.
- Revisione indipendente dei numeri di S1, S2, S3b non fatta (prima di G1).
- PAPER «Dati tenuti da parte» dice «H1 intera (S3b)», STATO «H1 letta intera»: chiarire a G1.
