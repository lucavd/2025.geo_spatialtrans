---
type: Closure
title: "Chiusura di geo_spatialtrans dopo il riesame"
description: "Che cosa abbiamo misurato, perche' non basta, che cosa si riusa, quali vincoli riaprono il progetto."
status: draft
generated: { by: geo-agent/skill-v2, at: 2026-11-20T17:30:00+01:00 }
pi_approval: pending
sources:
  - id: regola
    resource: "/riesame/REGOLA_B.md"
    title: "Riesame, Parte B: regola di esito scritta prima del calcolo"
  - id: fatti
    resource: "/riesame/FATTI_B.md"
    title: "Riesame, Parte B: i fatti sul generatore"
---
# Verdetto del riesame: chiudere

Regola di esito scritta prima (`45fdbbe`, approvata da Luca il 2026-11-20)[^regola]; fatti riletti dai file[^fatti].

# Che cosa abbiamo misurato

Sul tumore del colon e sul suo stroma, 18 024 geni, stessa catena di analisi (bin2cell, RCTD, etichette trasferite), distanza D dalla sezione di taratura in unita' di h:

| | tumore + stroma | tumore | stroma |
|---|---|---|---|
| tetto: altri due pazienti veri | 2,23 | 2,10 | 2,35 |
| simulatore di oggi | 7,83 | 9,89 | 4,99 |
| braccio banale (Poisson, Voronoi, tipi a caso, stessa cattura, spostamento spento) | 7,77 | 9,83 | 4,91 |
| simulatore meno banale (w) | +0,06 (0,29), uguale | +0,05 (0,44), uguale | +0,09 (0,22), uguale |
| simulatore meno tetto (w) | +5,60 (0,50), diverso | +7,78 (0,83), diverso | +2,65 (0,38), diverso |

Cancelli superati: ricalcolo della catena (7 su 7 entro 1e-9), due pazienti veri non diversi, stimatore senza distorsione, revisione indipendente (G-5), md5 dei 550 file letti invariati.

# Perche' non basta

1. Il generatore dello strato cellulare non batte la soluzione banale, in nessuno dei due tessuti. La regola dice: si toglie.
2. Simulatore e banale restano a circa 3,5 volte il tetto, con gli stessi scarti sulle misure piu' lontane (singoletti circa 10 h in tutti e due i tessuti, concordanza con il nucleo nel tumore circa 17 h): l'errore sta in cio' che hanno in comune (profili di espressione, cattura, profondita'), che questo confronto non separa.
3. La domanda sulla cattura e' indecidibile con i dati disponibili (decisione di Luca, XD-26): le registrazioni Xenium pubblicate da 10x appartengono a pacchetti che non sono su GEO.

# Che cosa si riusa

- Il quadro di validazione: tetto fra pazienti veri, braccio banale come pavimento, stimatore quadratico corretto con bootstrap per aree, margine di equivalenza e controlli scritti prima (`tools/riesame_D.py`, `tools/riesame_tetto.py`).
- Il braccio banale (`tools/riesame_B_sim.R`), piu' semplice del simulatore e non peggiore.
- Le misure di riferimento sui Visium HD veri dei sei tessuti e i dati di verifica: rene e fegato restano sigillati.
- Il codice della cattura (strisce, fondo, profondita'), da considerare non validato.

# Che cosa riapre il progetto

- Per la cattura: i pacchetti Xenium da circa 1,1 milioni di cellule (P1 1 123 632, P2 1 105 615, P5 1 029 090) a cui appartengono le registrazioni di 10x; con quelli la Parte B si fa come scritta (bracci X1, XC, XC2, cancelli G-1 e G-6).
- Per un simulatore: un nuovo contratto con una misura principale, il banale come pavimento e il tetto fra pazienti, che attacchi la voce comune a simulatore e banale (profili e cattura). Non e' una terza prova di questo generatore.

[^regola]: [Riesame, Parte B: regola di esito scritta prima del calcolo](/riesame/REGOLA_B.md).
[^fatti]: [Riesame, Parte B: i fatti sul generatore](/riesame/FATTI_B.md).
