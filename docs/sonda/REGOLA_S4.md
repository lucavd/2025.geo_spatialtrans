---
type: Rule
title: "S4, componenti veri: il banale con un componente alla volta preso da H2, per dividere il divario dal tetto per causa"
description: "Regola scritta e committata prima di simulare i bracci a componenti veri: bracci, codice, metro invariato, cancelli, controprova, attese, esito."
status: stable
generated: { by: geo-agent/progetto-master, at: 2026-11-23T22:40:00+01:00 }
sources:
  - id: s3
    path: docs/sonda/REGOLA_S3.md
  - id: s3c
    path: docs/sonda/REGOLA_S3c.md
  - id: per
    path: results/sonda/S3c/per_measure.csv
  - id: stato
    path: docs/STATO.md
---
# Perché

Oggi il simulatore è il banale B di `REGOLA_S3.md`[^s3]: D(B) = 3,151 h, tetto 0,336, soglia 1,336 (`REGOLA_S3c.md`[^s3c]). Il bilancio per voce (`per_measure.csv`[^per]) dice dove sta il divario (singoletti RCTD 60 %, UMI assegnati alle cellule 19 %, vicini dello stesso tipo 16 %) ma non quale pezzo del simulatore lo produce. Questa sessione fa il bilancio per componente con il primo metodo di progetto-master: si sostituisce un componente di B con la sua versione perfetta, presa dall'area omologa di H2, e si misura quanto scende D. Affermazione servita: A1. La prossima mossa approvata in `STATO.md`[^stato] è questa.

Stato al commit: nessun braccio nuovo simulato. Codice scritto e provato: `tools/sonda/S4_sim_arms.py`, `S4_test.py` (16/16), `S4_D.py` (a tabella vuota ritrova S3c con scarto 0).

# Bracci

Ogni braccio è B (`tools/sonda/S3_sim_B.py`) con alcuni interruttori su «vero». Componenti:

- geometria: B = centri di Poisson, Voronoi, nucleo a disco; vero = nuclei StarDist di V, territori di Voronoi dei loro baricentri, nucleo vero dentro la sua cellula;
- tipi: B = estratti dalla composizione, indipendenti dalla posizione; veri = tipo trasferito della cellula vera tipizzata più vicina (con la geometria vera è la cellula stessa);
- densità di UMI: B = Poisson(ρ × territorio), trascritti a caso nel territorio; vera = ogni quadratino riceve esattamente gli UMI di V;
- espressione: B = profilo medio del tipo; vera = profilo della cellula vera sorgente (conteggi del suo nucleo allargato in V), profilo del tipo se ha meno di 20 UMI.

| braccio | geometria | tipi | densità | espressione |
|---|---|---|---|---|
| Y | B | vera | B | B |
| G | vera | B | B | B |
| U | B | B | vera | B |
| A | vera | vera | vera | B |
| AE | vera | vera | vera | vera |

Y, G, U sono le sostituzioni singole. A e AE sono la catena con i componenti veri migliori, quella che il contratto chiede di leggere a G1. In Y nessun altro pezzo è allineato a V, quindi il guadagno di Y viene dalla struttura della mappa dei tipi (raggruppamento e composizione locale), che un generatore può imparare. In A e AE M-1 e M-2 coincidono con V per costruzione (stessi nuclei, stessi UMI per quadratino): lo si scrive, non è un risultato.

Aree: le 10 di V (5 tumore, 5 stroma); un seme (k = 0), `default_rng([20261129, braccio, i, 0])`; 50 unità. Un seme invece di due: in S3c la differenza fra i due semi di B è 0,11 h. Confronti contro B (media dei 2 semi) e contro B0 (un seme, come i bracci nuovi).

# Catena e metro (invariati)

Catena di `REGOLA_S3.md`: `tools/sonda/S3_chain.sh` con `S3_ROOT=/mnt/micron/geo_spatialtrans/S4_sonda`, stesso riferimento HLCA, stesso seme RCTD. Tabella con `tools/sonda/S3_table.py` in `results/sonda/S4/area_values.csv`. D con `tools/sonda/S4_D.py`: funzioni di `S3c_D.py` senza modifiche, stessi indici bootstrap (seme 20261123, B = 2000), metodo (ii) di S3c (f = 1,584). Confronti con la regola di S3c: uguale se |Δ| + w ≤ 1 h, diverso se |Δ| − w > 1 h. Quota del divario chiusa dal braccio a: (q_B − q_a) / (q_B − q_T), q = D con segno al quadrato; per voce, la stessa quota cella per cella.

# Cancelli (se uno fallisce, i numeri non si leggono)

| id | controllo | esito richiesto |
|---|---|---|
| G-T | `tools/sonda/S4_test.py`: identità dei componenti veri su una griglia piccola, determinismo, e tutti gli interruttori su B con il seme di B rigenerano bit per bit la simulazione B esistente | 16/16 PASS |
| G-ID | in ogni unità (`S4_check.json`): UMI per quadratino = V dove la densità è vera; nuclei = nuclei veri dove la geometria è vera; tipo = tipo della sorgente dove i tipi sono veri | 50/50 |
| G-S3c | `S4_D.py` ritrova D e intervalli di T, B, P, B0, B1 di `results/sonda/S3c/D.json` | scarto ≤ 1e-12 |
| G-MD5 | md5 dei file letti (`tools/sonda/S4_md5.sh`) prima e dopo | uguali |
| G-REV | un sotto-agente ricalcola dalla tabella D dei bracci nuovi e le quote chiuse con codice proprio | entro 1e-9 |

# Controprova

CP-AE: AE ha geometria, tipi, densità ed espressione veri; resta solo come sono disposti i geni dentro la cellula e fra cellule vicine. Deve arrivare al tetto: D(AE) ≤ soglia e AE contro T non «diverso». Se fallisce, una parte del divario non sta in nessuno dei componenti sostituiti (disposizione dei geni, spillover fra cellule, o la catena stessa), e il bilancio per componente è incompleto: la si dichiara come voce residua, pari a D(AE).

# Attese (scritte prima; tumore + stroma, metodo (ii))

| id | attesa | giudizio |
|---|---|---|
| P1 | Y: D(Y) fra 1,5 e 2,6 h; la quota q di tumore M-4 e di M-6 (entrambi i tessuti) scende almeno della metà rispetto a B | confermata se l'intervallo di D(Y) sta dentro [1,5; 2,6] e le tre voci scendono della metà; smentita se l'intervallo sta tutto fuori |
| P2 | G: la geometria vera da sola non sposta D: \|D(G) − D(B)\| ≤ 0,5 h, G contro B non «diverso» | come sopra, su \|Δ\| |
| P3 | U: la densità vera da sola non aiuta (i suoi picchi non cadono sui nuclei simulati): D(U) ≥ D(B) − 0,5 h | come sopra |
| P4 | A: D(A) ≤ pavimento − δ = 2,15 h | confermata se l'estremo superiore ≤ 2,15 |
| P5 | AE (controprova): D(AE) ≤ 1,336 e AE contro T non «diverso» | confermata se l'estremo superiore ≤ 1,336 |

Precisione attesa: semi-ampiezza dell'intervallo di D circa 0,3 h (D(B): 2,908; 3,440), più larga per i bracci con un seme.

# Esito e conseguenze

- Livello bilancio: il bilancio dell'errore in `STATO.md` si rifà per componente con le quote chiuse da Y, G, U (sostituzioni singole) e con la voce residua D(AE). Le sessioni di costruzione andranno al componente con la quota singola più grande; un componente sotto il 10 % non riceve sessioni (progetto-master).
- Livello obiettivo, solo come numero per G1 (non si decide qui): D con i componenti veri migliori = min(D(A), D(AE)). Il contratto lo legge a G1 insieme ad A3 e A5. Se anche AE resta sopra pavimento − δ (2,15 h), lo si scrive come condizione di chiusura possibile.
- A1 non cambia: il simulatore resta B, A1 = 3,151. Esito della sessione nel registro: «misura».

# Calcolo

Sul server convivono altri utenti: al più 12 processi a un thread in tutto per questa corsa (file `WORKERS`), sotto il tetto di 50 core e il limite termico (circa 28). Pilota: un'unità per braccio (tumore 1), con i tempi; se la corsa completa supera 4 ore con 12 processi ci si ferma e si chiede a Luca. Uscite in `/mnt/micron/geo_spatialtrans/S4_sonda/`, un file per unità, rilanciabile.

[^s3]: docs/sonda/REGOLA_S3.md
[^s3c]: docs/sonda/REGOLA_S3c.md
[^per]: results/sonda/S3c/per_measure.csv
[^stato]: docs/STATO.md
