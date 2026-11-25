---
type: State
title: "Stato di geo_spatialtrans"
description: "Memoria fra le sessioni: fase, dati, numeri, bilancio dell'errore, minacce, strade fallite, prossima mossa."
status: draft
generated: { by: geo-agent/progetto-master, at: 2026-11-25T20:30:00+01:00 }
---
# Fase

Sonda, budget 9 (proroga 7 -> 9, Luca 2026-11-25), usate 8: S1, S2, S5 (`infra`), S3b, S3c, S4, S6, S8 (`misura`); RC (`riesame`) non conta. Resta S9: A4 con A2. Poi G1; nessuna terza proroga.

# Dati (polmone 10x, CG000709, 1 blocco FFPE, 1 donatore)

- H1 (T), H2 (V) in `/mnt/micron/geo_spatialtrans/data_lung/`; md5 in `results/sonda/S3/MD5_download_H.txt`. HLCA core in `ref_lung/`. Aree omologhe di H1 in `S3b_sonda/runs/T`.
- Metà tenuta da parte: solo H2 (x < m − 50 µm, m = 3372), X1 (m = 3300), X2. H1 letta intera.
- Orientamento (S5): X1 come H2, H1 ruotata di 180°, X2 di circa 5,7°.

# Metro (S3c, `docs/sonda/REGOLA_S3c.md`)

`tools/sonda/S3c_D.py`, bootstrap per coppie stratificato, B = 2000, seme 20261123, metodo (ii) f = 1,584.

# Numeri (h, tumore + stroma)

Pavimento D(B) 3,151 (2,908; 3,440), tetto D(T) 0,336 (−0,263; 0,821), soglia 1,336 (δ = 1). A1 = 3,151: non regge.

Componenti veri (S4): Y 1,771; G 3,014; U 2,636; A 0,961; AE 0,509. Tipi non copiati (S6): Yalt 1,735 (1,190; 2,489), YH1 1,786, Yshuf 2,637; regola RIESAME_C §6: caso 4 con δ = 1.

# A3 tolta (Luca, 2026-11-25, S8)

Non misurabile su X1 alla profondità di H2 (S5). Fra i limiti del contratto.

# A5 (S8, `docs/sonda/REGOLA_S8.md`, `results/sonda/S8/A5.json`)

Metodi (Luca): ENACT 0.2.5 naive, area, gene, cluster; bin2cell. Coseno per cellula sui 284 geni di G_p contro la verità (B: territori; X1: cellula Xenium abbinata); tau di Kendall su metà disgiunte delle aree; margine 0,2 (Luca).

- T5 (X1 con sé) 0,996 ± 0,079; S5c (B contro X1) 0,398 ± 0,038; A5 −0,598 ± 0,067 (regola: non regge).
- CP-shuf fallita: T5 con verità permutata 0,872 ± 0,176. Per la regola A5 non si legge come fedeltà: non scritta in `results/misura.json`.
- Esplorativo, non pre-registrato: guadagno = coseno − coseno permutato: T5 1,000, S5c 0,400, A5 −0,600. Disaccordo: b2c ultimo in B, terzo in X1; E-area ed E-cluster scambiati (3 coppie su 10).
- Cancelli: G-ENACT 30/30 (E-cluster: peso 1/n di ENACT alle cellule senza quadratini propri, errore di pre-registrazione corretto prima della corsa), G-TRUTH 20/20, G-MATCH 10/10, G-ID 2/2, G-REV (max 7,8e-16), md5 invariati.
- X1: abbinati 176-984 nuclei StarDist su 1288-2781 per area; tenute 151-837 cellule.

# Check

`check smoke` = `tools/sonda/S3c_smoke.sh` (circa 5 min). Non contiene S4, S5, S6, S8. `check full` da scrivere.

# Bilancio dell'errore (S4 + S6)

1. Disposizione dei tipi: 69-71 % del divario, non viene dal copiare (Yalt = Y). Agisce su M-4.
2. UMI concentrati sui nuclei (M-1): con M-5 dello stroma quasi tutto il resto. Coerente con A5: b2c raccoglie meno UMI del pannello per cellula in B (mediana 7 contro 9-11 di ENACT) che in X1 (15 contro 14-19,5).
3. Espressione 7 %, geometria dei nuclei 9 %: nessuna sessione.

# Minacce alla validità

1. HLCA è tessuto sano: le cellule tumorali non hanno un tipo (M-4).
2. Le tipizzazioni di V e di H1 vengono dalla stessa catena che misura; A4 con una tipizzazione diversa.
3. δ = 1 h; con 5 aree per tessuto e una sola permutazione (S6) molte prove restano indeterminate.
4. Un donatore: tetto fra sezioni seriali.
5. A5: su X1 si abbina circa un nucleo su cinque; il coseno sul pannello non è specifico (CP-shuf).

# Strade fallite

- Generatore dello strato cellulare (colon): 7,83 contro 7,77 del banale.
- Tetto di S3 su aree non omologhe; stimatore indipendente su aree omologhe (S3b).
- A3 sul pannello di X1 alla profondità di H2.
- X2: registrazione non giudicabile (S1, S2).
- Correlazione di fase su tessere da 1 mm per X1: usare l'ECC.
- A5 con il coseno grezzo sul pannello: controprova fallita (S8).

# Decisioni aperte per Luca

1. Punteggio di A5 = guadagno sulla verità permutata, pre-registrato e calcolato sui file di S8 (minuti).

# Prossima mossa

Sessione 9: A4 con A2 (catena con tipizzazione diversa da RCTD con HLCA); in apertura, se Luca approva, regola di A5 con il guadagno. Poi G1: D(Yalt) 1,735 nella fascia «riposizionare» (intervallo a cavallo di 2,151); A5 sul banale oltre il margine.

# Aperti

- `results/F3.2/sheets/` non tracciato: lasciato com'è.
- `tools/okf.py check` rosso da SC1.
- Revisione indipendente dei numeri di S1, S2, S3b non fatta (prima di G1).
- REGOLA_S3c: 1,96 contro z = 1,959964 in f; nessun esito cambia.
- PAPER «H1 intera (S3b)» contro STATO «H1 letta intera»: chiarire a G1.
- `.venv-enact` (ENACT 0.2.5, setuptools 79 per StarDist, pyarrow) non ricostruibile da script: da aggiungere a `tools/`.
