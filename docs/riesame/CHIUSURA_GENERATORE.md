---
type: Closure
title: "Chiusura del generatore dello strato cellulare e della taratura attuale dello spostamento"
description: "Il riesame chiude due componenti del simulatore di oggi, non l'obiettivo: che cosa abbiamo misurato, perche' quei componenti non bastano, che cosa si riusa, che cosa resta aperto."
status: stable
generated: { by: geo-agent/skill-v2, at: 2026-11-20T17:11:38+01:00 }
pi_approval: { by: human:luca-vedovelli, at: 2026-11-20, evidence: "in chat: «Non approvo la chiusura del progetto: il riesame ha giudicato il simulatore di oggi, non l'obiettivo. [...] si chiudono il generatore dello strato cellulare e la taratura attuale dello spostamento. Nessuna terza prova su quel codice. L'obiettivo resta aperto.»" }
sources:
  - id: regola
    resource: "/riesame/REGOLA_B.md"
    title: "Riesame, Parte B: regola di esito scritta prima del calcolo"
  - id: fatti
    resource: "/riesame/FATTI_B.md"
    title: "Riesame, Parte B: i fatti sul generatore"
---
# Decisione

Si chiudono il generatore dello strato cellulare (regioni, semina, tessellazione, nuclei, disposizione dei tipi) e la taratura attuale dello spostamento dei trascritti. Nessuna terza prova su quel codice. L'obiettivo, un simulatore di dati Visium HD per provare metodi di analisi, resta aperto: la ripartenza comincia dalla fase di scommessa (`docs/STATO.md`). Decisione di Luca del 2026-11-20, dopo la regola di esito scritta prima[^regola] e i fatti riletti dai file[^fatti].

# Che cosa abbiamo misurato

Tumore del colon e suo stroma, 18 024 geni, stessa catena di analisi (bin2cell, RCTD, etichette trasferite); distanza D dalla sezione di taratura in unita' di h:

| | tumore + stroma | tumore | stroma |
|---|---|---|---|
| tetto: altri due pazienti veri | 2,23 | 2,10 | 2,35 |
| simulatore di oggi | 7,83 | 9,89 | 4,99 |
| braccio banale (Poisson, Voronoi, tipi a caso, stessa cattura, spostamento spento) | 7,77 | 9,83 | 4,91 |
| simulatore meno banale (w) | +0,06 (0,29), uguale | +0,05 (0,44), uguale | +0,09 (0,22), uguale |
| simulatore meno tetto (w) | +5,60 (0,50), diverso | +7,78 (0,83), diverso | +2,65 (0,38), diverso |

Cancelli superati: ricalcolo della catena (7 su 7 entro 1e-9); due pazienti veri non distinguibili come richiede il cancello (P1 meno P5 = -0,42, w 0,68: esito indeterminato, non «diverso»; non vuol dire «uguali»); stimatore senza distorsione; revisione indipendente (G-5); md5 dei 550 file letti invariati.

# Perche' quei componenti non bastano

1. Il generatore non batte la soluzione banale, in nessuno dei due tessuti. La regola dice: si toglie.
2. Sostituire geometria e disposizione dei tipi con la versione banale non cambia la distanza dal vero (7,83 contro 7,77). I vicini dello stesso tipo restano lontani dal vero in tutti e due i bracci: nel tumore 10,9 con il simulatore e 12,3 con il banale, contro 0,2 fra due pazienti veri.
3. Nel tumore i marcatori estranei sono a 8,1 con lo spostamento (simulatore) e a 3,3 senza (banale), con il tetto a 2,3. Valori per misura, senza incertezza (`results/riesame/B_per_measure_common.csv`).
4. La domanda sulla cattura e' indecidibile con i dati disponibili (decisione di Luca, XD-26): le registrazioni Xenium pubblicate da 10x appartengono a pacchetti che non sono su GEO.

# Che cosa si riusa

- Il metro: D con lo stimatore quadratico corretto, bootstrap per aree, margine di equivalenza, controlli positivo e negativo (`tools/riesame_D.py`, `tools/riesame_tetto.py`, ricalcolo indipendente in `results/riesame/G5/`).
- Il braccio banale (`tools/riesame_B_sim.R`) come pavimento.
- Le misure di riferimento sui Visium HD veri dei sei tessuti; rene e fegato restano sigillati.
- Il codice della cattura (strisce, fondo, profondita'), non validato.

# Che cosa resta aperto

- L'obiettivo: ripartenza dalla scommessa, con un contratto nuovo.
- La cattura: si riapre con i pacchetti Xenium da circa 1,1 milioni di cellule (P1 1 123 632, P2 1 105 615, P5 1 029 090) a cui appartengono le registrazioni di 10x.

[^regola]: [Riesame, Parte B: regola di esito scritta prima del calcolo](/riesame/REGOLA_B.md).
[^fatti]: [Riesame, Parte B: i fatti sul generatore](/riesame/FATTI_B.md).
