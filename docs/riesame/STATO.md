---
type: State
title: "Riesame: stato dopo la Parte B"
description: "Parte B chiusa: verdetto chiudere in attesa dell'ok di Luca; dove stanno numeri, script e dati."
status: draft
generated: { by: geo-agent/skill-v2, at: 2026-11-20T17:30:00+01:00 }
---
# Dove

- Branch `riesame`. Regola `docs/riesame/REGOLA_B.md`; scelte `DECISIONI_B.md` (fino a B-13, C-10), `DECISIONI_B_XENIUM.md` (fino a XD-26).
- Fatti: `docs/riesame/FATTI_B.md`. Verdetto: `docs/CHIUSURA.md` (pi_approval pending).
- Numeri: `results/riesame/B_D_common.json`, `B_area_values.csv`, `B_per_measure_common.csv`, `B_G4.csv`, `B_h_panel.csv`, `xen_match.csv`, `xen_same_acquisition.json`, `G5/`.
- Dati: `/mnt/micron/geo_spatialtrans/riesame/B` (bracci V, T, S, B sui geni comuni; simulazioni banali); Xenium P2 completo e `cells.parquet` di P1, P5 in `/mnt/micron/geo_spatialtrans/data_xenium/`.

# Esito

- Generatore: S uguale a B (+0,06, w 0,29), S peggio del tetto (+5,60, w 0,50): si toglie.
- Cattura: indecidibile con i dati disponibili (XD-26). Bracci Xenium e bracci del pannello non calcolati.
- Verdetto: chiudere. Candidati: A (cattura) sospeso sul vincolo dei pacchetti Xenium; B cade; C vale.

# Prossima mossa

Ok di Luca su `CHIUSURA.md`. Nessun altro calcolo.

# Trappole

- La pagina dei dataset 10x risponde 429 / pagina anti-bot ai client automatici.
- Simulazione banale A2_r2 seme 2: triangolazione del pacchetto fallita (B-13).
- `results/F3.2/sheets/` non e' nostro: non aggiungerlo.
