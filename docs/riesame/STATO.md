---
type: State
title: "Riesame: stato prima della Parte B"
description: "Dove stanno script e dati, i passi della Parte B in ordine, le trappole viste."
status: draft
generated: { by: geo-agent/skill-v2, at: 2026-11-20T13:00:00+01:00 }
---
# Dove

- Branch `riesame`. Regola: `docs/riesame/REGOLA_B.md` (approvata `45fdbbe`; G-3 corretto `cf4d6b5`). Lettura dei dati di verifica annotata: `38f640c`.
- Analisi: `tools/riesame_D.py`, provato sui dati di F3.3. Tetto: `tools/riesame_tetto.py`.
- Catena del pannello e braccio banale: `tools/riesame_panel_*`, `tools/riesame_B_*`, scelte in `docs/riesame/DECISIONI_B.md` (`9bc87a7`).
- Bracci Xenium: `tools/riesame_xen_*`, `docs/riesame/DECISIONI_B_XENIUM.md`: bozza di un sotto-agente interrotto, non rivista.
- Dati: `/mnt/micron/geo_spatialtrans/riesame/{xenium,B}`; zip Xenium di P2 parziale (5,5 su 15,9 GB) in `data_xenium/P2/`.

# Passi della Parte B

1. Rileggere la regola e le due tabelle di scelte; rivedere quella Xenium prima di usarla.
2. Finire lo scaricamento (`curl -C -`), controllare l'integrità, estrarre.
3. G-4: cancellare `B/G4/runs/A3_r1_common` (interrotto), rilanciare `riesame_panel_run.sh B/jobs/g4.tsv`.
4. V, T, S, B sul pannello; B anche sui geni comuni.
5. Registrazione (regola XD-4), esportazioni X, X1, XC, XC2, catena.
6. `riesame_D.py`: G-1, G-2, G-3, G-6; revisione G-5; due difensori; verdetto; `PAPER.md` o `CHIUSURA.md`; esito in `HOLDOUT_LOG.tsv`.

Nessun braccio calcolato e nessun dato di verifica letto finora.

# Trappole

- I processi lanciati sul server da un sotto-agente sopravvivono alla chat: fermarli per gruppo di processi.
- Il falso confronto 3 + 2 di G-3 era vuoto per costruzione: vedi la correzione nella regola.
- Sul pannello i singoletti di V superano 0,99: M-4 esce da D.
- `call_command`: 60 s al massimo, uscita 64 KB; invii a pezzi di 5 000 caratteri.
- Server condiviso: al più 4 processi.
- `results/F3.2/sheets/` non è nostro: non aggiungerlo.
