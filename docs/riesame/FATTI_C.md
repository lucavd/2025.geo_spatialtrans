---
type: Facts
title: "Riesame C: i fatti della sonda sul polmone (S1-S5), riletti dai file"
description: "Solo numeri riletti dai file del repo al 2026-11-24 (commit 7ab28d0): misura, pavimento, tetto, soglia, componenti veri, A3, budget, esistente."
status: draft
generated: { by: geo-agent/progetto-master, at: 2026-11-24T13:10:00+01:00 }
sources:
  - { id: paper, path: docs/PAPER.md }
  - { id: pj, path: progetto.json }
  - { id: d3c, path: results/sonda/S3c/D.json }
  - { id: d4, path: results/sonda/S4/D.json }
  - { id: pm4, path: results/sonda/S4/per_measure.csv }
  - { id: r4, path: docs/sonda/REGOLA_S4.md }
  - { id: s5, path: results/sonda/S5/diag_summary.csv }
  - { id: mis, path: results/misure.tsv }
  - { id: es, path: docs/ESISTENTE.md }
  - { id: rb, path: docs/riesame/REGOLA_B.md }
---
Non verificato da Luca.

# 1. Obiettivo e misura

- Messaggio del contratto: «Un simulatore di Visium HD con verità per cellula riproduce un tumore polmonare quanto una seconda sezione e ordina i metodi come la verità»[^paper]. Lettore: chi sviluppa o confronta metodi Visium HD.
- Dati: un blocco FFPE di polmone 10x, un donatore. H1 e H2 = due sezioni seriali solo Visium HD. X1 = Xenium v1 (289 geni) con Visium HD sulla stessa sezione; X2 = Xenium 5K con Visium HD sulla stessa sezione.
- Misura D (unità h, meglio verso il basso): si fa girare la stessa catena (bin2cell, poi RCTD con il riferimento HLCA del polmone sano) sul simulato e su H2 (V). Per 2 tessuti (tumore, stroma) × 6 misure (M-1 quota di UMI assegnati alle cellule, M-2 UMI per cellula, M-3 marcatori estranei, M-4 quota di singoletti RCTD, M-5 concordanza del tipo RCTD, M-6 vicini dello stesso tipo) si calcola, su 5 aree appaiate per tessuto, q = (media(d)² − var(d)/n)/h², d = braccio − V per area. D = radice con segno della media delle 12 q. h = semiampiezza dell'intervallo bootstrap al 95% della media di V sulle 5 aree (quindi 1 h è circa il rumore di campionamento del bersaglio stesso). Unità indipendente: l'area; bootstrap per coppie, B = 2000[^d3c].
- δ = 1 h. Fonte scritta: «circa metà del tetto di tumore e stroma (2,12 h) misurato in Parte A»[^rb], cioè con una versione precedente di h (con i κ tarati sul colon) e sui geni comuni. Da S3b h non ha più i κ e il tetto attuale è 0,336 h. Sensibilità prevista dal contratto: δ = 0,5 e 2.
- Soglia = tetto + δ = 1,336 h[^pj].

# 2. Pavimento, tetto, controllo positivo (S3c)[^d3c]

| braccio | che cos'è | D (h) | intervallo 95% |
|---|---|---|---|
| T (tetto) | H1, aree omologhe a quelle di V | 0,336 | −0,263; 0,821 |
| B (pavimento) | simulatore banale: centri di Poisson, Voronoi, nucleo a disco, tipi estratti dalla composizione senza posizione, UMI Poisson uniformi nel territorio, profilo medio del tipo | 3,151 | 2,908; 3,440 |
| B0, B1 | B con i due semi separati | 3,096; 3,205 | B0 contro B1: «uguale» (Δ −0,11, w 0,12) |
| P (controllo positivo) | banale peggiorato apposta | 4,019 | 2,983; 5,472 |

- Pavimento − tetto = 2,81 h, cioè circa 2,8 δ.
- Voci di B (q, h²): M-4 tumore 45,6; M-4 stroma 25,7; M-1 tumore 14,4; M-6 stroma 10,3; M-6 tumore 8,27; M-1 stroma 8,26; M-5 stroma 5,51; le altre 5 sotto 1. Somma 119,1. La quota di singoletti RCTD (M-4) pesa il 60% di D²; M-1 il 19%[^pm4].
- Quota di singoletti in V: tumore 0,812, stroma 0,775; in B 0,555 e 0,485[^pm4].

# 3. Componenti veri (S4)[^r4][^d4]

Ogni braccio è B con uno o più componenti presi da H2, cioè dallo stesso bersaglio contro cui si misura D. Nel braccio Y ogni cellula simulata riceve il tipo della cellula vera tipizzata più vicina della stessa area di V: copia la disposizione dei tipi del bersaglio, non solo le sue statistiche.

| braccio | componenti veri | D (h) | intervallo 95% | quota del divario chiusa |
|---|---|---|---|---|
| Y | tipi | 1,771 | 1,514; 2,200 | 0,69 (0,56; 0,74) |
| G | geometria dei nuclei | 3,014 | 2,766; 3,351 | 0,09 (0,00; 0,17) |
| U | UMI per quadratino | 2,636 | 1,565; 3,582 | 0,30 (0,02; 0,60) |
| A | geometria + tipi + UMI | 0,961 | 0,691; 1,315 | 0,92 (0,88; 0,96) |
| AE | A + espressione della cellula sorgente | 0,509 | 0,359; 0,652 | 0,99 (0,97; 1,03) |

- Per tessuto: Y tumore 1,853, stroma 1,686; A tumore 0,548, stroma 1,244; AE 0,528 e 0,488.
- Y azzera M-4 e M-6 (q da 45,6/25,7/10,3/8,27 a 3,4/−0,07/0,19/−0,06) e lascia M-1 com'era (14,9 e 8,78): in Y M-1 pesa il 63% di D². Somma delle q di Y 37,65; senza le due M-1 resterebbe 14,0, cioè D ≈ 1,08 h (calcolo aritmetico, non un braccio simulato)[^pm4].
- In A e AE M-1 e M-2 coincidono con V per costruzione (stessi nuclei, stessi UMI per quadratino)[^r4].
- Revisione: un sotto-agente ha ricalcolato D dei bracci in S4 (`results/sonda/S4/review_S4.csv`, 235 righe, scarto 0). In questo riesame D di ogni braccio è stato riottenuto dalle 12 q di `per_measure.csv` (stesso valore alla quarta cifra).

# 4. Regola di G1 scritta nel contratto[^paper]

«G1, con i componenti veri o spenti migliori: pavimento e tetto distanti meno di δ: chiudere; D ≤ soglia: costruzione; soglia < D ≤ pavimento − δ: riposizionare sull'esito B; D > pavimento − δ: chiudere.» Letta oggi con A (0,961) o AE (0,509), dice «costruzione».

# 5. Le altre affermazioni

- A2 (nessuna voce oltre il tetto più largo per voce): non misurata.
- A3 (il simulatore batte lo Xenium diviso in quadratini, stessa sezione, geni del pannello): non misurabile come scritta su X1. Sui 284 geni del pannello H2 ha il 2% degli UMI, 1-7 UMI per blocco da 8 µm (mediana). RCTD (UMI minimo 20) non parte in B in 20 unità su 20, nel controllo positivo in 9 su 10, in V in 3 su 10 (tutte stroma)[^s5]. Il braccio Xenium (XB) è costruito e registrato (G-XR 10/10, residuo dopo correzione ≤ 0,41 µm). Lo Xenium è circa 10 volte più profondo di H2 sui geni del pannello.
- A3 su X2: la registrazione di X2 non è giudicabile (S1: errore standard 1,1-3,4 µm; S2: nessun metro con errore standard ≤ 0,5 µm; il filtro per gene tiene 0-34% delle UMI di X2): primo colpo su quell'ipotesi[^mis].
- A4 (catena mai usata in taratura): tetto non misurato.
- A5 (classifica dei metodi): tetto definito come concordanza X1 contro X2; con X2 non registrata il tetto non è definito.

# 6. Storia e budget

- Prima del polmone (colon): il generatore dello strato cellulare costruito in 43 sessioni dava 7,83 contro 7,77 del banale, tetto 2,23; chiuso con decisione di Luca (`docs/riesame/CHIUSURA_GENERATORE.md`). L'obiettivo è rimasto aperto; ripartenza dal pavimento nel polmone.
- Registro dall'adozione: SC1, SC2 (scommessa), S1, S2, S3b, S3c, S4, S5 (sonda), 2026-11-20 → 2026-11-24[^mis].
- Sonda: budget 5, prorogato a 7 il 2026-11-23 (atteso nella proroga: D ≤ 2,15 h; altrimenti G1 con i numeri raggiunti); usate 6. Luca il 2026-11-24 ha chiesto una seconda proroga, da 7 a 9; lo strumento chiede prima questo riesame[^pj].
- A1 (il simulatore oggi è B): 3,151. Nessun pezzo costruito nel polmone.
- Costruzione proposta nel contratto: 10 sessioni.

# 7. Esistente[^es]

- Ricerca del 2026-11-20 (Europe PMC, arXiv, GitHub, PyPI, Bioconductor, web, Zotero): non trovato un simulatore con griglia 2 µm, forme cellulari reali, spostamento, cattura e verità quadratino → cellula. Parti separate: ENACT (Xenium in quadratini da 2 µm), ALBIS (molecole con cellula d'origine, cellule sferiche, nessuna cattura), STHD, SPARKLE, modello delle strisce di Malsot. SpatialSimBench: solo spot, scala relativa senza tetto vero contro vero.
- Controllo leggero del 2026-11-24 (2 ricerche web): nessun lavoro nuovo con lo stesso messaggio. Un preprint del 2026-10-27 su Visium HD nel Crohn osserva che il supporto di RCTD non valida confini né identità delle cellule.

# 8. Minacce scritte in STATO

1. HLCA è tessuto sano: le cellule tumorali non hanno un tipo (pesa su M-4).
2. Y porta anche la composizione locale dell'area: un generatore impara le statistiche, non la mappa.
3. Un donatore: il tetto è fra sezioni seriali.
4. A3 senza misura.
5. U chiude il 30% del divario senza tipi veri (intervallo 2-60%): non spiegato.

[^paper]: docs/PAPER.md
[^pj]: progetto.json
[^d3c]: results/sonda/S3c/D.json
[^d4]: results/sonda/S4/D.json
[^pm4]: results/sonda/S4/per_measure.csv
[^r4]: docs/sonda/REGOLA_S4.md
[^s5]: results/sonda/S5/diag_summary.csv
[^mis]: results/misure.tsv
[^es]: docs/ESISTENTE.md
[^rb]: docs/riesame/REGOLA_B.md
