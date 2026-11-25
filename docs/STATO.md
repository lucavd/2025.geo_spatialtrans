---
type: State
title: "Stato di geo_spatialtrans"
description: "Memoria fra le sessioni: fase, dati, numeri, bilancio dell'errore, minacce, strade fallite, prossima mossa."
status: draft
generated: { by: geo-agent/progetto-master, at: 2026-11-24T17:30:00+01:00 }
---
# Fase

Sonda, budget 7, usate 7: S1, S2, S5 (`infra`), S3b, S3c, S4, S6 (`misura`); RC (`riesame`) non conta. Proroga 7 -> 9 approvata da Luca il 2026-11-25 dopo l'esito di S6 (caso 4), registrata in `progetto.json`. Scopo: prima misura di A5 e di A4 con A2. Altrimenti G1 con i numeri di oggi; nessuna terza proroga.

# Dati (polmone 10x, CG000709, 1 blocco FFPE, 1 donatore)

- H1 (T), H2 (V) in `/mnt/micron/geo_spatialtrans/data_lung/`; md5 in `results/sonda/S3/MD5_download_H.txt`. HLCA core in `ref_lung/`. Aree omologhe di H1 in `S3b_sonda/runs/T` (non `S3_sonda/runs/T`, aree vecchie).
- Metà tenuta da parte: solo H2 (x < m − 50 µm, m = 3372), X1 (m = 3300), X2. H1 letta intera.
- Orientamento (S5): X1 come H2, H1 ruotata di 180°, X2 di circa 5,7°.

# Metro (S3c, `docs/sonda/REGOLA_S3c.md`)

`tools/sonda/S3c_D.py`, bootstrap per coppie stratificato, B = 2000, seme 20261123, metodo (ii) f = 1,584. Ricalcolato in RC e in S6 da sotto-agenti.

# Numeri (h, tumore + stroma)

Pavimento D(B) 3,151 (2,908; 3,440), tetto D(T) 0,336 (−0,263; 0,821), soglia 1,336 (δ = 1). A1 = 3,151: non regge.

Componenti veri copiati da H2 (S4, `results/sonda/S4/D.json`): Y (tipi) 1,771 (1,514; 2,200); G 3,014; U 2,636; A 0,961; AE 0,509.

Tipi non copiati (S6, `results/sonda/S6/D.json`, `docs/sonda/REGOLA_S6.md`): Yalt (tipi dall'area successiva di V) 1,735 (1,190; 2,489), quota chiusa 0,71; YH1 (tipi dall'area omologa di H1) 1,786 (1,504; 2,259), 0,69; Yshuf (tipi dell'area permutati) 2,637 (2,324; 3,076), 0,30. Yalt e YH1 contro Y: «uguali» (Δ −0,04 e +0,02). Regola RIESAME_C §6: caso 4 con δ = 1 (estremi superiori 2,489 e 2,259 > 2,151), caso 1 con δ = 0,5.

# A3 (S5)

Non misurabile su X1 alla profondità di H2: pannello 2 % degli UMI, RCTD non parte in B 20/20. RC consiglia di toglierla; decide Luca.

# Check

`check smoke` = `tools/sonda/S3c_smoke.sh` (circa 5 min). Non contiene S4, S5, S6 (`S4_test.py` 16/16 e `S6_gates_sim.py` a parte). `check full` da scrivere.

# Bilancio dell'errore (S4 + S6)

1. Disposizione dei tipi: 69-71 % del divario, e non viene dal copiare (Yalt = Y). Agisce su M-4: tumore q 45,6 in B, 4,2 in Yalt, 31,1 in Yshuf. Si impara da un'altra area dello stesso tessuto.
2. UMI concentrati sui nuclei (M-1): in Yalt 76 % di D² nel tumore, 54 % nello stroma; con M-5 dello stroma (48 %) è quasi tutto il resto. Torna solo con nuclei e UMI veri insieme (S4).
3. Espressione 7 %, geometria dei nuclei 9 %: nessuna sessione.

# Minacce alla validità

1. HLCA è tessuto sano: le cellule tumorali non hanno un tipo (M-4).
2. Le tipizzazioni di V e di H1 vengono dalla stessa catena (RCTD con HLCA) che misura: M-4 e M-6 in parte costruiti; A4 con una tipizzazione diversa.
3. δ senza fonte valida dopo S3b; la precisione non basta a decidere con δ = 1 (anche Y ha estremo superiore 2,200 > 2,151).
4. Un donatore: tetto fra sezioni seriali.
5. Una sola permutazione e un solo ciclo di aree in S6 (5 aree per tessuto).

# Strade fallite

- Generatore dello strato cellulare (colon): 7,83 contro 7,77 del banale.
- Tetto di S3 su aree non omologhe e κ del colon; stimatore indipendente su aree omologhe (S3b): tetto −0,77.
- A3 sul pannello di X1 alla profondità di H2.
- X2: registrazione non giudicabile (S1, S2).
- Correlazione di fase su tessere da 1 mm per X1: usare l'ECC.

# Decisioni aperte per Luca (G1)

1. δ = 1 h con la fonte nuova (precisione dell'area vera), sensibilità 0,5 (RIESAME_C §7).
2. A3 tolta con il fatto dichiarato.

# Prossima mossa

Sessione 8: A5 su X1 (tetto riscritto senza X2). Sessione 9: A4 con A2. Poi G1. Numero del simulatore senza copia: D(Yalt) 1,735, nella fascia «riposizionare» del contratto (fra soglia 1,336 e pavimento − δ 2,151) con intervallo a cavallo di 2,151. In costruzione la prima voce sarebbe M-1 (UMI sui nuclei), poi un generatore della mappa dei tipi che impari la disposizione da un'altra area. A3, A4, A5 non misurate: lo strumento non lascia entrare in costruzione senza una prima misura.

# Aperti

- `results/F3.2/sheets/` non tracciato: lasciato com'è.
- `tools/okf.py check` rosso da SC1.
- Revisione indipendente dei numeri di S1, S2, S3b non fatta (prima di G1).
- REGOLA_S3c scrive 1,96 in f, il codice z = 1,959964 (f 1,5837842 contro 1,5837551): estremi spostati di 1e-5 h, nessun esito cambia (revisione S6).
- PAPER «Dati tenuti da parte» dice «H1 intera (S3b)», STATO «H1 letta intera»: chiarire a G1.
