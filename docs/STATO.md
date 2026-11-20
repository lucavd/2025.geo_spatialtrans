---
type: State
title: "Stato per la ripartenza di geo_spatialtrans"
description: "Metro, pavimento, tetto, primo bilancio dell'errore, che cosa manca per un contratto."
status: draft
generated: { by: geo-agent/skill-v2, at: 2026-11-20T17:11:38+01:00 }
---
# Fase

Scommessa, in una chat nuova. Chiusi generatore dello strato cellulare e taratura attuale dello spostamento (`docs/riesame/CHIUSURA_GENERATORE.md`).

# Metro

D = radice con segno della media su 12 voci (tumore, stroma x M-1...M-6) di (Delta^2 - Var Delta)/h^2; Delta = braccio meno P2; h = precisione della media delle 5 aree di P2. M-1 UMI assegnati, M-2 UMI per cellula, M-3 marcatori estranei, M-4 singoletti RCTD, M-5 concordanza con il nucleo, M-6 vicini dello stesso tipo. `tools/riesame_D.py`; controlli in `docs/riesame/FATTI_B.md`.

# Numeri (h)

Pavimento 7,77 (braccio banale); tetto 2,23 (P1, P5 contro P2); simulatore 7,83.

# Primo bilancio dell'errore

Per misura, non per causa: quota di D^2(simulatore) - D^2(tetto), senza incertezza (`results/riesame/B_per_measure_common.csv`).

| voce | tetto | simulatore | banale | quota |
|---|---|---|---|---|
| tumore, concordanza con il nucleo | 1,9 | 17,1 | 17,7 | 43 % |
| tumore, vicini dello stesso tipo | 0,2 | 10,9 | 12,3 | 18 % |
| stroma, singoletti RCTD | 1,9 | 10,4 | 10,6 | 15 % |
| tumore, singoletti RCTD | 1,3 | 10,2 | 10,2 | 15 % |
| tumore, marcatori estranei | 2,3 | 8,1 | 3,3 | 9 % |
| stroma, vicini dello stesso tipo | 3,6 | 5,6 | 5,4 | 3 % |

Altre 6 voci: da -2 % a 0 %.

# Manca per un contratto

- Messaggio, lettore, rivista; soglia con fonte fra 2,23 e 7,77.
- 3-5 affermazioni con figura e numero; esiti A, B, C; budget.
- Esistente: nessun simulatore esistente passato nel metro.
- Tenuti da parte: rene e fegato sigillati; P1, P5, 6,5 mm letti 2 volte.
- Cattura: pacchetti Xenium da circa 1,1 milioni di cellule.
