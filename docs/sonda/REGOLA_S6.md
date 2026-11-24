---
type: Rule
title: "S6, tipi non copiati: il braccio Y con la mappa dei tipi presa da un'altra area o dall'altra sezione"
description: "Regola scritta e committata prima di simulare: bracci, codice, metro invariato, cancelli, controprova, attese, esito (regola di RIESAME_C par. 6)."
status: stable
generated: { by: geo-agent/progetto-master, at: 2026-11-24T13:40:00+01:00 }
sources:
  - id: rc
    path: docs/riesame/RIESAME_C.md
  - id: r4
    path: docs/sonda/REGOLA_S4.md
  - id: d4
    path: results/sonda/S4/D.json
  - id: stato
    path: docs/STATO.md
---
# Perché

In S4 il braccio Y (il banale B con il tipo di ogni cellula preso dalla cellula vera tipizzata più vicina) ha chiuso il 69 % del divario fra pavimento e tetto: D(Y) = 1,771 h contro D(B) = 3,151[^d4]. Ma Y prende i tipi dall'area che poi si misura: il guadagno può venire dal copiare il bersaglio (minaccia 2 di `STATO.md`[^stato]). Il riesame RC ha scritto la regola di esito di questa prova (`RIESAME_C.md` par. 6[^rc]); Luca ha dato l'ok alla prova il 2026-11-24. Affermazione servita: A1, attraverso G1. Sessione 7 della sonda, l'ultima del budget firmato.

# Bracci

Tutti sono Y di S4 (`REGOLA_S4.md`[^r4]): geometria, densità di UMI ed espressione del banale; tipo = tipo della cellula sorgente tipizzata più vicina. Cambia solo la sorgente dei tipi.

| braccio | sorgente dei tipi | che cosa prova |
|---|---|---|
| Yalt | area successiva dello stesso tessuto in H2 (k → k mod 5 + 1), stesso riferimento locale | se bastano le statistiche del tessuto, cioè ciò che un generatore può imparare |
| YH1 | area omologa di H1 (`S3b_sonda/runs/T`, tipizzata dalla sua catena), baricentri portati in H2 con la catena di coordinate di S3b | un simulatore che prende in ingresso una sezione vera |
| Yshuf | cellule tipizzate dell'area misurata, tipi permutati a caso | controprova: composizione dell'area giusta, disposizione distrutta |
| Yself | l'area misurata (come Y) | solo cancello: deve rigenerare Y bit per bit; non passa dalla catena |

Seme della simulazione = seme di Y in S4, `default_rng([20261129, 0, i, 0])`: nuclei e UMI per cellula sono identici a Y, cambia solo il tipo. Permutazione di Yshuf `default_rng([20261124, i])`. Un seme, 10 aree per braccio, 30 unità nella catena. Codice: `tools/sonda/S6_sim_types.py` (con l'argomento `type_source` aggiunto a `S4_sim_arms.simulate`, valido solo con geometria B), `S6_gates_sim.py`, `S6_D.py`, `S6_run.sh`, `S6_md5.sh`.

# Catena e metro (invariati)

Catena di S3 (`tools/sonda/S3_chain.sh`, bin2cell + RCTD con HLCA, seme RCTD invariato) con `S3_ROOT=/mnt/micron/geo_spatialtrans/S6_sonda`. Tabella con `S3_table.py` in `results/sonda/S6/area_values.csv`. D con `S6_D.py`: funzioni di `S3c_D.py` senza modifiche, stessi indici bootstrap (seme 20261123, B = 2000), metodo (ii) di S3c (f = 1,584). Confronti con la regola di S3c (uguale se |Δ| + w ≤ 1 h, diverso se |Δ| − w > 1 h).

# Cancelli (se uno fallisce, i numeri non si leggono)

| id | controllo | esito richiesto |
|---|---|---|
| G-CODE | `S4_test.py` dopo la modifica; Yself rigenera Y di S4 (conteggi per quadratino e gene, nuclei, tipi) | 16/16 e 10/10 |
| G-MAP | la mappa H2 → H1 ritrova i centri riportati di `areas_H1_hom.csv`; andata e ritorno; trasformazione affine | scarti ≤ 1e-6 µm |
| G-ID | Yalt, YH1, Yshuf: nuclei, UMI totali e numero di cellule identici a Y; tipo = tipo della sorgente | 30/30 |
| G-INT | `S6_D.py` ritrova D e intervalli di S3c e S4 (`results/sonda/S4/D.json`); md5 dei file letti prima e dopo | ≤ 1e-12; uguali |
| G-REV | un sotto-agente ricalcola dalla tabella D dei tre bracci nuovi e il caso della regola con codice proprio | entro 1e-9, stesso caso |

# Regola di esito (copiata da RIESAME_C par. 6, livello: obiettivo, attraverso G1)

Tumore + stroma, metodo (ii). Con δ = 1 h, pavimento − δ = 2,151:

1. Estremo superiore di D(Yalt) ≤ 2,151: la disposizione si impara dalle statistiche. Proroga da 7 a 9 (A5 con il tetto riscritto su X1, poi A4 con A2), poi G1.
2. Altrimenti, se l'estremo superiore di D(YH1) ≤ 2,151: serve una sezione vera in ingresso. A G1 si propone di restringere il messaggio a un simulatore che parte da una sezione reale; proroga come in 1.
3. Estremi inferiori di D(Yalt) e di D(YH1) entrambi > 2,151: il guadagno dei tipi veniva dal copiare il bersaglio. G1 letta sui bracci non copiati dice «chiudere»; si va a G1 senza proroga.
4. Ogni altro caso: indeterminato; a G1 con i numeri, senza proroga automatica.

Sensibilità con δ = 0,5 (soglia 2,651), riportata accanto, non decide. La proroga resta una decisione di Luca: la regola dice che cosa si propone.

# Nota scritta prima del calcolo: la precisione e la soglia

A tabella vuota `S6_D.py` ritrova S3c e S4 con scarto 0 e dà l'intervallo di Y: 1,514; 2,200. L'estremo superiore di Y supera già 2,151. Quindi un braccio buono quanto Y cade nel caso 2 o nel caso 4, non nel caso 1: il caso 1 chiede che Yalt faccia meglio di Y. Un esito indeterminato è probabile per costruzione. L'intervallo viene dalla variabilità fra le 5 aree per tessuto, non dal seme: più semi non lo stringono. La regola non si cambia (l'ha scritta il riesame, Luca l'ha approvata). Per leggerla si aggiunge un confronto appaiato, che non decide l'esito:

- P4: Yalt contro Y e YH1 contro Y con la regola di S3c (uguale se |Δ| + w ≤ 1 h, diverso se |Δ| − w > 1 h). «Uguale» vuol dire che, alla precisione delle 5 aree, togliere la copia non costa più di δ; «diverso» con il braccio peggiore vuol dire che la copia contava.

# Controprova

CP-SHUF: Yshuf ha la composizione esatta dell'area misurata ma i tipi messi a caso. Se la disposizione conta, Yshuf deve stare lontano da Y: Yshuf contro Y «diverso», con Yshuf peggiore. Se fallisce (Yshuf uguale o vicino a Y), il guadagno di Y viene dalla composizione dell'area e non dalla disposizione: la voce del bilancio «disposizione dei tipi» si rinomina «composizione locale», e un esito 1 o 2 si legge come composizione, non come disposizione.

# Attese (scritte prima)

| id | attesa | giudizio |
|---|---|---|
| P1 | D(Yalt) circa 2,3 h, fra D(Y) e D(B) | numero per il registro (atteso della sessione); l'esito è quello della regola |
| P2 | D(YH1) ≤ D(Yalt): la sezione seriale ha composizione locale più vicina all'area misurata | confermata se YH1 contro Yalt non è «diverso» con YH1 peggiore e D(YH1) < D(Yalt) |
| P3 | CP-SHUF passa: D(Yshuf) ≥ 2,5 h | confermata se l'estremo inferiore ≥ 2,5; smentita se l'estremo superiore < 2,5 |

A1 non cambia (il simulatore resta B, A1 = 3,151). Esito della sessione nel registro: «misura».

# Calcolo

Sul server gira un'altra corsa (circa 6 processi R da 3 core). Al più 10 processi a un thread per questa catena (file `WORKERS`), sotto il tetto di 50 core e il limite termico (Tctl controllata durante la corsa). 30 unità da circa 37 minuti: circa 2 ore. Prima Yalt e YH1, poi Yshuf. Uscite in `/mnt/micron/geo_spatialtrans/S6_sonda/`, un file per unità, rilanciabile.

[^rc]: docs/riesame/RIESAME_C.md
[^r4]: docs/sonda/REGOLA_S4.md
[^d4]: results/sonda/S4/D.json
[^stato]: docs/STATO.md
