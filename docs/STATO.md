---
type: State
title: "Stato di geo_spatialtrans"
description: "Memoria fra le sessioni: fase, numeri, bilancio dell'errore, minacce, strade fallite, prossima mossa."
status: draft
generated: { by: geo-agent/skill-v2, at: 2026-11-20T19:00:00+01:00 }
---
# Fase

Scommessa. Sessione SC1 (2026-11-20): adozione di progetto-master, contratto in bozza in `docs/PAPER.md`, revisore simulato fatto. Attesa: scelta del tessuto primario (SC2), poi firma (G0).

# Metro

D = radice con segno della media su 12 voci (tumore, stroma x M-1...M-6) di (Delta^2 - Var Delta)/h^2; Delta = braccio meno P2; h = precisione della media delle 5 aree di P2. `tools/riesame_D.py`; controlli in `docs/riesame/FATTI_B.md`.

# Numeri (h)

Pavimento 7,77 (banale); tetto 2,23 (P1, P5 contro P2); simulatore ereditato 7,83; soglia proposta 3,23 (tetto + 1 h).

# Bilancio dell'errore (per misura, non ancora per causa)

Quota di D^2(simulatore) - D^2(tetto), `results/riesame/B_per_measure_common.csv`, senza incertezza.

| voce | tetto | simulatore | banale | quota |
|---|---|---|---|---|
| tumore, concordanza con il nucleo | 1,9 | 17,1 | 17,7 | 43 % |
| tumore, vicini dello stesso tipo | 0,2 | 10,9 | 12,3 | 18 % |
| stroma, singoletti RCTD | 1,9 | 10,4 | 10,6 | 15 % |
| tumore, singoletti RCTD | 1,3 | 10,2 | 10,2 | 15 % |
| tumore, marcatori estranei | 2,3 | 8,1 | 3,3 | 9 % |
| stroma, vicini dello stesso tipo | 3,6 | 5,6 | 5,4 | 3 % |

Raggiungibilita' (ricalcolata dallo stesso file): le 4 voci maggiori portate a 2,0 h danno D = 3,28 (sopra 3,23); le 5 maggiori a 2,0 h danno 2,37; a 3,6 h danno 3,71 e 3,05.

# Minacce alla validita'

1. Circolarita': il metro usa bin2cell + RCTD, che il simulatore deve aiutare a provare (A4: catena mai usata in taratura).
2. Una sola sezione di taratura (P2): A1 si conferma con P1 e P5 come riferimento.
3. h dipende dal disegno a 5 aree: tutti i numeri valgono con quel disegno.
4. Realismo non vuol dire classifica fedele (Crowell 2023); nessuna verita' vera nel colon.
5. Sul pannello Xenium i singoletti RCTD non hanno spazio (99,9 %): A3 puo' restare indecidibile.

# Strade fallite

- Generatore dello strato cellulare: 7,83 contro 7,77 del banale, uguale (`docs/riesame/CHIUSURA_GENERATORE.md`).
- Taratura attuale dello spostamento: marcatori estranei del tumore 8,1 con, 3,3 senza.
- Cattura contro Xenium in quadratini: indecidibile, i pacchetti Xenium registrati (circa 1,1 milioni di cellule) non sono su GEO.

# Prossima mossa

Sessione SC2 (scommessa, ultima del budget): scelta del tessuto primario, polmone o colon. Luca (2026-11-20): una sola decisione presa soppesando pro e contro, non due contratti. Nessun codice; si scaricano solo metadati e file piccoli.

1. Inventario del polmone dopo Xenium (10x): registrazione Xenium -> Visium HD fornita?, geni (Xenium v1, Prime 5K), dimensioni dei file, regioni disgiunte utilizzabili (sviluppo, tetto, dati tenuti da parte), qualita' del Visium HD.
2. Colon di Long et al. 2025 (Zenodo 15722536): molecole o solo file Loupe?
3. Tabella unica pro e contro, sui criteri scritti qui prima dell'inventario: verita' per cellula sulla stessa sezione; tetto (fra pazienti o solo fra aree della stessa sezione); geni; qualita' del dato; costo di rifare pavimento, tetto e taratura; dati tenuti da parte per la conferma; misura di utilita' possibile (classifica dei metodi contro la verita').
Regola SC2 (scritta 2026-11-20 prima dell'inventario). Eliminatori: (E1) il tessuto deve permettere un tetto misurato con lo stesso metro, fra pazienti o, in mancanza, fra regioni disgiunte della stessa sezione dichiarate tali; (E2) servono dati tenuti da parte per A1 che non siano la sezione di taratura. Un tessuto che fallisce E1 o E2 si scarta. Fra i tessuti rimasti decide la misura di utilita' (verita' per cellula sulla stessa sezione, registrazione con residuo al 95 percentile <= 2 um come XD-4); a parita', il costo di rifare pavimento, tetto e taratura.

4. Una raccomandazione, decisa con Luca; poi il contratto si aggiorna e si firma (G0).

Dopo la firma, sonda: check che scrive `results/misura.json`; sostituzioni con i componenti veri (regola di G1 in `docs/PAPER.md`); braccio ENACT e seconda catena. Pulizie: `CLAUDE.md` sotto 100 righe, storici in `docs/archivio/`, controllo OKF su PAPER e STATO.

# Aperti

- SpatialSimBench letto per intero (riga in `docs/ESISTENTE.md`).
- `results/F3.2/sheets/` non tracciato (PNG del 18-19 novembre, di F3.2): lasciato com'e'.
