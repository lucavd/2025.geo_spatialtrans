---
type: State
title: "Stato di geo_spatialtrans"
description: "Memoria fra le sessioni: fase, dati, numeri, bilancio dell'errore, minacce, strade fallite, prossima mossa."
status: draft
generated: { by: geo-agent/progetto-master, at: 2026-11-27T19:30:00+01:00 }
---
# Fase

Sonda, budget 9 (proroga 7 -> 9, Luca 2026-11-25), usate 9: S1, S2, S5 (`infra`), S3b, S3c, S4, S6, S8, S9 (`misura`); RC (`riesame`) non conta. Prossima: G1. Nessuna terza proroga.

# Dati (polmone 10x, CG000709, 1 blocco FFPE, 1 donatore)

- H1 (T), H2 (V) in `/mnt/micron/geo_spatialtrans/data_lung/`; md5 in `results/sonda/S3/MD5_download_H.txt`. HLCA core e LuCA core (md5 058d7069…) in `ref_lung/`.
- Metà tenuta da parte: solo H2 (x < m − 50 µm, m = 3372), X1 (m = 3300), X2. H1 letta intera.
- Orientamento (S5): X1 come H2, H1 ruotata di 180°, X2 di circa 5,7°.

# Metro (S3c, `docs/sonda/REGOLA_S3c.md`)

`tools/sonda/S3c_D.py`, bootstrap per coppie stratificato, B = 2000, seme 20261123, metodo (ii) f = 1,584.

# Numeri (h, tumore + stroma)

Catena di oggi (bin2cell, RCTD e classificatore con HLCA): pavimento D(B) 3,151, tetto 0,336, soglia 1,336 (δ = 1). A1 = 3,151: non regge.

S4: Y 1,771; G 3,014; U 2,636; A 0,961; AE 0,509. S6: Yalt 1,735 (1,190; 2,489), YH1 1,786, Yshuf 2,637.

# S9 (`docs/sonda/REGOLA_S9.md`, `results/sonda/S9/`)

- A4, catena U (RCTD e classificatore con LuCA, 13 tipi con «Tumor cells»): tetto 0,328, D(B) 4,436, e = +3,11 ± 0,40: non regge.
- A4, catena C (CellTypist, modello HLCA): tetto 0,186, D(B) 2,293, e = +1,11 ± 0,49: non regge (e − 2 se = +0,14). A4 non regge.
- CP-chain, C-4, C-5 passati con U e C.
- A2 (catena di oggi): T_max 0,806 (stroma M-1); B oltre in 9 voci su 12 (6 nette, 6 indeterminate, 0 dentro), se 1,95: non regge.
- A5 col guadagno (Luca 2026-11-26): T5 1,000 ± 0,060, S5c 0,400 ± 0,031, A5 −0,600 ± 0,066: non regge. CP-gain: tau con gli UMI assegnati 0,6 (< 0,8). Guadagno in X1 0,19-0,21, in B 0,58-0,70.
- Cancelli passati; G-REF e G-B2C solo per errore di pre-registrazione (REGOLA_S9).

# Check

`check smoke` = `S3c_smoke.sh` + `S8_to_misura.py` + `S9_to_misura.py` (circa 5 min). Non ricalcola S4-S9. `check full` da scrivere.

# Bilancio dell'errore (S4 + S6 + S9)

1. Disposizione dei tipi (M-4, M-6): 69-71 % del divario con la catena di oggi. Dipende dalla catena: con LuCA i singoletti fanno l'85 % (tumore 11,2 h), con CellTypist quasi spariscono (tumore 0,47 h) e domina M-6 dello stroma (55 %).
2. UMI concentrati sui nuclei (M-1): uguale in tutte le catene; con C è il 35 %.
3. Espressione 7 %, geometria dei nuclei 9 %: nessuna sessione.

# Minacce alla validità

1. Tipi tumorali: LuCA non li risolve. RCTD sceglie «Tumor cells» come primo tipo nell'1,5 % dei blocchi delle aree tumorali di V (Epithelial cell 56 %).
2. Il valore di D e la sua voce più grande dipendono dalla catena di misura (4,44; 3,15; 2,29).
3. δ = 1 h; con 5 aree per tessuto molte prove restano indeterminate.
4. Un donatore: tetto fra sezioni seriali.
5. A5: su X1 si abbina circa un nucleo su cinque; il guadagno in X1 è piccolo.

# Strade fallite

- Generatore dello strato cellulare (colon): 7,83 contro 7,77 del banale.
- Tetto di S3 su aree non omologhe; stimatore indipendente su aree omologhe (S3b).
- A3 sul pannello di X1 alla profondità di H2.
- X2: registrazione non giudicabile (S1, S2).
- Correlazione di fase su tessere da 1 mm per X1: usare l'ECC.
- A5 con il coseno grezzo (S8); LuCA come rimedio ai tipi tumorali (S9).

# G1 (2026-11-27, `docs/sonda/G1.md`, `results/G1/`)

- Revisione indipendente: D di 13 bracci e M-1…M-6 dal basso su 10/10 aree, scarto 0; registrazione di X1 rifatta dal grezzo, p95 0,97/0,87/0,80 µm, passa. Testo: 284 geni, non 286; quota UMI di X1 42-90 %.
- Yalt tocca tutte e tre le fasce di G1 (1,190; 2,489); caso 4 con δ = 1.
- Proposta a Luca: costruzione, messaggio ristretto (disposizione dei tipi da una sezione vera non misurata, nessun generatore della disposizione), budget 5. Prima prova: Yalt + UMI addensati sul nucleo stimati dai dati, D e A5. Uscita: estremo inferiore di D > 2,151 → proposta di chiudere.

# Prossima mossa

Risposta di Luca su G1. Poi `chiudi --sessione G1 --riesame`, `cancello` con le sue parole, PAPER (messaggio, budget, decisione), commit.

# Aperti

- `results/F3.2/sheets/` non tracciato: lasciato com'è.
- `tools/okf.py check` rosso da SC1.
- REGOLA_S3c: 1,96 contro z = 1,959964 in f; nessun esito cambia.
- `.venv-enact` (ENACT, CellTypist) non ricostruibile da script; modelli CellTypist in `~/.celltypist` (61, scaricati in S9).
- `results/misure.tsv`: la riga S9 registra solo A4; A5 col guadagno e A2 solo in `misura.json`.
