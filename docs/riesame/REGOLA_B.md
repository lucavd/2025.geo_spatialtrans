---
type: Preregistration
title: "Riesame, Parte B: regola di esito scritta prima del calcolo"
description: "Bracci, misura, soglie di «uguale» e «diverso», cancello sul pannello Xenium e corrispondenza fra esiti e verdetto, scritti prima di ogni calcolo della Parte B."
status: draft
generated: { by: geo-agent/skill-v2, at: 2026-11-20T13:30:00+01:00 }
pi_approval: { by: human:luca-vedovelli, at: 2026-11-20, evidence: "in chat: «1. Sì: Parte B come scritta, con il braccio X1. 2. Sì: margine di equivalenza 1,0 h e onere della prova sui nostri componenti. Aggiungi un controllo positivo per la cattura (G-6)»" }
sources:
  - id: dec098
    resource: "/decisions/DEC-098.md"
    title: "F3.3: spostamento e degrado con lo stimatore quadratico corretto"
  - id: dec093
    resource: "/decisions/DEC-093.md"
    title: "Regole contro il sovradattamento"
  - id: fatti
    resource: "/riesame/FATTI.md"
    title: "Riesame, Parte A: i fatti riletti dai file"
---
# Regola di esito della Parte B

Scritta e committata prima di ogni calcolo della Parte B (prima versione `ca73e42`). Approvata da Luca il 2026-11-20 in chat: «1. Sì: Parte B come scritta, con il braccio X1. 2. Sì: margine di equivalenza 1,0 h e onere della prova sui nostri componenti. Aggiungi un controllo positivo per la cattura (G-6): XC con lo spostamento raddoppiato deve risultare «diverso» da X1. Se non lo è, sul pannello la domanda sulla cattura è «indecidibile», non «non aggiunge».» Questa versione aggiunge G-6 e precisa i bracci B, X, XC prima di ogni calcolo.

## Domanda

Che cosa dà questo simulatore che non si ottiene dividendo in quadratini da 2 µm un dato Xenium vero?

## Tessuti, unità, geni

- Tumore del colon e suo stroma; il riferimento V è la sezione di taratura P2.
- Unità indipendente: l'area. Per i bracci simulati, prima la media dei 2 semi per area (come F3.3).
- Geni principali: i 408 geni del pannello Xenium dentro i 18 024 comuni, uguali per tutti i bracci. Secondari: i 18 024 comuni, per i bracci senza Xenium (V, T, S, B).
- Le misure dipendono dai geni, quindi su ogni insieme di geni si ricalcolano V, T e il margine h (costruzione di F3.1 sulle 5 aree di V).

## Bracci

| braccio | che cos'è | aree |
|---|---|---|
| V | Visium HD P2 | le 5 aree di taratura per tessuto |
| T | Visium HD P1 e P5 (tumore da 6,5 mm come sensibilità, paziente incerto, DEC-093 R-8) | le aree di verifica di F3.3 |
| S | il simulatore di oggi (simulazioni di F3.1, nessuna nuova) | 5 aree × 2 semi |
| B | banale: centri da un processo di Poisson omogeneo con la densità del preset; territori di Voronoi; nucleo = cerchio con l'area nucleare media del preset, ritagliato nel territorio; tipi estratti a caso dalla composizione del preset; trascritti per cellula dal profilo medio del tipo con la profondità del preset, messi a caso dentro la cellula; cattura dello Step 2 (profondità, strisce, fondo) con lo spostamento spento | stesse 5 aree (geometria) × 2 semi |
| X | Xenium P2: trascritti con qualità almeno 20 sommati in quadratini da 2 µm sulla griglia Visium HD (registrazione pubblicata da 10x); nuclei Xenium dati a bin2cell | 5 aree nelle posizioni delle aree di V |
| X1 | come X, con un diradamento binomiale unico che porta gli UMI per quadratino in tessuto alla media di V | come X |
| XC | trascritti Xenium P2 con qualità almeno 20, poi lo spostamento dello Step 2 con i parametri del preset del tessuto, poi la cattura dello Step 2 (strisce, fondo) e un fattore di profondità unico che porta gli UMI per quadratino in tessuto alla media di V; nuclei Xenium dati a bin2cell | come X, × 2 semi; solo se costa meno di una sessione |
| XC2 | come XC, con la scala dello spostamento raddoppiata (controllo positivo G-6) | come XC |

X1 è stato aggiunto dal riesame: lo Xenium grezzo perderebbe contro V per la sola profondità, che chiunque corregge con un numero. Il modello della cattura deve battere X1.

Scelte di attuazione che il codice impone (per esempio come spegnere lo spostamento, come leggere i parametri del preset) si registrano come decisioni tecniche, con il commit, prima del calcolo del braccio interessato; nessuna si decide guardando D.

## Misura

D(braccio) = radice con segno della media, sulle celle tessuto × misura M-1...M-6, di (Δ² − Var Δ) / h², con Δ = media del braccio meno media di V (stimatore quadratico corretto di DEC-098). Più basso è meglio. D(T) è il tetto, media su P1 e P5.

- Incertezza: bootstrap parametrico per aree (modello t con la deviazione standard combinata; i bracci simulati ridisegnati con il loro errore standard fra aree), B = 2 000, seme 20261123.
- Confronto fra due bracci a e b: Δab = D(a) − D(b), con w = semi-ampiezza al 95 % della sua distribuzione bootstrap.
- La misura dei singoletti (M-4) si toglie da D prima di vedere gli altri bracci se V ha almeno 0,99 di singoletti sui geni del pannello: lì non ha spazio (in Parte A: 0,9994 su un'area di tumore).

## Soglie

- Margine di equivalenza δ = 1,0 h. Fonte: circa metà del tetto di tumore e stroma (2,12 h, intervallo 1,78-2,55), misurato in Parte A. Rumore atteso di una differenza fra due bracci con 5 aree: w da 0,42 a 0,67 h (Parte A, dati di F3.3). Con δ = 1 «uguale» richiede |Δab| sotto circa 0,4-0,6 e «diverso» sopra circa 1,4-1,7: tutti e due possibili.
- Uguale: |Δab| + w ≤ δ. Diverso: |Δab| − w > δ (meglio o peggio secondo il segno). Altrimenti indeterminato.
- Se per un confronto che decide w supera δ, l'esito è indeterminato per costruzione e lo si scrive.
- Onere della prova: un nostro componente vale solo con un esito «diverso» a suo favore. Un esito indeterminato conta come «non aggiunge» per quel componente.

## Cancelli e controlli (prima di confrontare i bracci)

| id | controllo | esito richiesto | se fallisce |
|---|---|---|---|
| G-1 | controllo positivo sul pannello: D(S) oltre D(T) di più di δ, come sui geni completi (8,0 contro 2,1) | diverso, S peggio | il pannello ha perso i difetti noti: la domanda sul generatore si decide sui geni completi con V, T, S, B, e la domanda sulla cattura è indecidibile sul pannello |
| G-2 | controllo negativo: P1 contro P5, D(P1) − D(P5) | non diverso | il rumore è sottostimato: nessun esito «uguale» vale |
| G-3 | falso confronto: le 5 aree di V divise 3 + 2 in tutti i 10 modi; D della parte da 2 contro quella da 3 | media entro ±0,5 h | lo stimatore è distorto: la Parte B non si legge |
| G-4 | ricalcolo: i bracci V, T, S sui geni completi ritrovano i valori di F3.3 entro 1e-9 relativo | identità | errore di catena: si corregge prima di tutto |
| G-5 | revisione avversariale: un sotto-agente ricalcola D e Δ di ogni braccio dalle tabelle per area | entro l'errore Monte Carlo | ogni discrepanza è un dubbio da chiudere prima del verdetto |
| G-6 | controllo positivo della cattura (Luca, 2026-11-20): D(XC2) − D(X1) | diverso | sul pannello la domanda sulla cattura è «indecidibile», non «non aggiunge» |

## Esiti e conseguenze

| domanda | risultato | conseguenza |
|---|---|---|
| Il modello della cattura serve? (solo se G-1 e G-6 passano) | D(X1) non diverso da D(T), cioè D(X1) − D(T) non oltre δ con il suo w | la strada corta basta: il simulatore non ha un lettore come sostituto del dato vero |
| | D(X1) oltre D(T), e D(XC) meglio di D(X1) (diverso) | il modello della cattura vale: è il candidato articolo |
| | D(X1) oltre D(T), e D(XC) uguale o indeterminato rispetto a D(X1) | il modello della cattura non aggiunge |
| Il generatore serve? | D(S) uguale a D(XC) | il generatore regge il confronto con la disposizione vera: si tiene |
| | D(S) peggio di D(XC) e meglio di D(B) (tutti e due diversi) | si tiene solo per prove controllate, con lo scarto dichiarato |
| | D(S) uguale, indeterminato o peggio rispetto a D(B) | il generatore non batte il banale: si toglie |

Se G-1 o G-6 non passano, la domanda sulla cattura è «indecidibile» sul pannello. Se XC non si fa, al suo posto si usa X1, e il verdetto lo dice. Se lo Xenium non è usabile, la Parte B si fa con V, T, S, B sui geni completi e decide solo sul generatore (D(T) al posto di D(XC)).

Verdetto, uno solo: chiudere, se la strada corta basta oppure se il modello della cattura non aggiunge; cambiare strada, se la cattura vale e il generatore no; continuare, se valgono tutti e due (contratto in `docs/PAPER.md`); restringere, se gli esiti del tumore e dello stroma, calcolati anche separatamente, portano a verdetti diversi: vale quello del tessuto dove la cattura vale.

Prima del verdetto, due sotto-agenti che non conoscono la storia del progetto difendono «continuare» e «chiudere» sui fatti completi; un terzo ricalcola pavimento, tetto e misura. Al più due sessioni: senza risposta dopo la seconda, il verdetto è «chiudere», con il motivo. Una domanda sulla cattura «indecidibile» alla fine della seconda sessione è una domanda senza risposta: il motivo e il vincolo che la riaprirebbe (per esempio uno Xenium con molti più geni) vanno in `docs/CHIUSURA.md`.

Calcolo: al più 4 processi in totale sul server condiviso (carico circa 23 e 79 °C il 2026-11-20 alle 12:18).

## Attese (scritte prima)

- D(T) sul pannello: circa 2 h, come sui geni completi (2,12).
- D(S): oltre D(T) di più di δ (G-1), come sui geni completi (8,0).
- D(X): oltre D(T), per la profondità. D(X1) e D(XC): nessuna attesa numerica; sono la domanda.
- D(B): fra D(T) e D(S) sulle misure di mescolanza, perché il banale non sposta i trascritti; nessuna attesa numerica.

## Dati di verifica

La Parte B legge per la seconda volta gli insiemi di P1 e P5 del tumore e dello stroma e del tumore da 6,5 mm (DEC-093, R-6: ogni insieme giudica al più 3 cambiamenti). Si annota in `data/holdout/HOLDOUT_LOG.tsv` prima della lettura, con il commit di questa regola. Rene e fegato restano sigillati. Lo Xenium non è un dato di verifica: è un dato nuovo, usato solo dentro i bracci X, X1, XC.
