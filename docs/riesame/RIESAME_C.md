---
type: Review
title: "Riesame C: la sonda del polmone regge? Criterio, pezzi, obiettivo e verdetto"
description: "Riesame prima della seconda proroga della sonda: critica del criterio, giudizio dei pezzi, quattro condizioni sull'obiettivo, avvocati e ricalcolo indipendente, verdetto e regola di esito della prossima prova."
status: draft
generated: { by: geo-agent/progetto-master, at: 2026-11-24T13:40:00+01:00 }
sources:
  - { id: fatti, path: docs/riesame/FATTI_C.md }
  - { id: r4, path: docs/sonda/REGOLA_S4.md }
  - { id: sim, path: tools/sonda/S4_sim_arms.py }
---
Non verificato da Luca. Fatti e numeri: `FATTI_C.md`[^fatti]. Nessun codice toccato.

# 1. Il criterio

1. Lo spazio c'è: pavimento − tetto = 2,81 h, circa 2,8 δ. Il pavimento è riproducibile (due semi «uguali») e il controllo positivo peggiora (4,02). Il metro funziona.
2. La lettura di G1 con i componenti veri è circolare. A e AE prendono nuclei, tipi e UMI dalla stessa area di H2 contro cui si misura D; AE è H2 risimulata da H2. Che stia al tetto è un controllo di coerenza della catena, non una prova che un generatore ci arrivi. G1 «costruzione» oggi non dice nulla sulla fattibilità. Correzione proposta: G1 si legge solo su bracci i cui componenti non vengono dalle aree misurate.
3. Circolarità anche nei tipi. Y copia l'etichetta trasferita prodotta dalla catena stessa (RCTD con HLCA) su V (`labels_transferred.parquet`[^sim]), cioè l'uscita da cui si calcolano M-4 (singoletti) e M-6 (vicini dello stesso tipo). Parte della caduta di M-4 e M-6 in Y è costruita. Il controllo giusto è A4 con una tipizzazione diversa da RCTD con HLCA.
4. δ = 1 h non ha più una fonte valida. La fonte scritta (metà del tetto di 2,12 h) è di una versione precedente di h; con la stessa logica oggi δ sarebbe 0,17 h, soglia 0,50 h, dentro l'incertezza del tetto stesso (estremo superiore 0,82): criterio irraggiungibile per costruzione, quindi quella logica non si può tenere. Una fonte che regge: 1 h è la semiampiezza dell'intervallo della media del bersaglio sulle 5 aree, quindi δ = 1 h vuol dire «un errore grande quanto la precisione con cui conosciamo l'area vera». Decide Luca.
5. Sensibilità del contratto: con δ = 2 la regola di G1 si contraddice (soglia 2,336 sopra pavimento − δ = 1,151: «costruzione» e «chiudere» si sovrappongono per ogni δ > 1,41 h). Con δ = 0,5 la soglia è 0,836 e A (0,961) non passa. La sensibilità ammissibile è solo verso il basso.
6. La misura è dominata da una voce: nel pavimento M-4 pesa il 60% di D², in Y M-1 pesa il 63%. Il quadro deve sempre mostrare le voci.

# 2. I pezzi

Nel polmone non c'è un pezzo costruito: il simulatore è il banale B (A1 = 3,151). Niente da togliere. Il generatore del colon è già stato tolto. Il bilancio per componente resta valido come diagnosi (dove sta il divario: disposizione dei tipi, UMI concentrati sui nuclei), non come prova che lo si possa chiudere.

# 3. L'obiettivo, sulle quattro condizioni

| condizione | vale? | fatto |
|---|---|---|
| 1. pavimento e tetto non distinguibili | no | 2,81 h contro δ = 1 h |
| 2. l'articolo esiste o non ha lettore | no | nessun simulatore con griglia 2 µm, forme reali e verità per cellula (ricerca del 2026-11-20; controllo leggero del 2026-11-24) |
| 3. i dati per misurarlo non esistono | in parte | per D sì (H1, H2); per A3 no (X1 sul pannello: RCTD non parte); per A5 il tetto scritto (X1 contro X2) non è definito, ma X1 da solo ha Visium HD e Xenium sulla stessa sezione: il tetto si può riscrivere |
| 4. budget finito senza avvicinarsi alla soglia | no | 6 sessioni su 7; nessun braccio non copiato misurato |

Nessuna condizione vale per la misura principale. Il costo già speso (43 sessioni nel colon) non entra nel giudizio.

# 4. Avvocati e ricalcolo

- «Continuare» (sotto-agente senza storia): continuare, condizionato a un braccio con tipi non copiati da H2 nella sessione 7 già nel budget; la proroga si decide dopo quel numero.
- «Chiudere» (sotto-agente senza storia): nessuna condizione di chiusura per D; propone di restringere il messaggio togliendo A3 e A5 e di non concedere proroga né costruzione prima di un braccio non copiato.
- Ricalcolo indipendente (sotto-agente, codice proprio, sola lettura): 30 valori di D (10 bracci × 3 ambiti) e 15 quote chiuse ritrovati con scarto ≤ 4·10⁻¹⁵. Intervalli riprodotti esattamente dallo schema di semi dichiarato; rispetto al bootstrap esatto (5⁵ ricampionamenti per tessuto) gli estremi superiori di B e A sono bassi di 0,025 e 0,03 h (circa 2 errori Monte Carlo): trascurabile rispetto a δ. Con il bootstrap esatto l'estremo superiore di A (1,346) supera la soglia (1,336): A è al limite, non sotto.

I due avvocati concordano su tre punti: G1 letta su bracci copiati è circolare; δ va riancorato; A3 come scritta va tolta.

# 5. Verdetto

**Continuare**, con il messaggio invariato e una condizione: la sessione 7 (l'ultima nel budget della sonda) misura la disposizione dei tipi non copiata dal bersaglio. La seconda proroga (da 7 a 9) si concede solo se questa prova passa.

# 6. Regola di esito della prossima prova (scritta ora, livello: obiettivo, attraverso G1)

Bracci: B con i tipi presi, per vicino più prossimo, da cellule tipizzate che non sono quelle dell'area misurata.

- Y-alt: tipi dall'area successiva dello stesso tessuto fra le 5 di V (permutazione ciclica, nessun dato nuovo). Prova se bastano le statistiche del tessuto, cioè ciò che un generatore può imparare.
- Y-H1: tipi dall'area omologa di H1 (l'altra sezione, tipizzata dalla sua catena). Prova un simulatore che prende in ingresso una sezione vera.

Stesso metro, stessi indici bootstrap, seme di B come in S4, intervallo con il metodo di S3c. Con δ = 1 h (pavimento − δ = 2,151):

1. Estremo superiore di D(Y-alt) ≤ 2,151: la disposizione si impara dalle statistiche. Proroga da 7 a 9 (A5 con il tetto riscritto su X1, poi A4 con A2), poi G1.
2. Altrimenti, se l'estremo superiore di D(Y-H1) ≤ 2,151: serve una sezione vera in ingresso. A G1 si propone di restringere il messaggio a un simulatore che parte da una sezione reale; proroga come in 1.
3. Estremo inferiore di D(Y-alt) e di D(Y-H1) entrambi > 2,151: il guadagno dei tipi veniva dal copiare il bersaglio. G1 letta sui bracci non copiati dice «chiudere» (D > pavimento − δ); si va a G1 senza proroga.
4. Ogni altro caso (intervalli a cavallo di 2,151): indeterminato; a G1 con i numeri, senza proroga automatica.

Le soglie cambiano solo se Luca cambia δ; con δ = 0,5 pavimento − δ = 2,651.

# 7. Decisioni per Luca

1. G1 letta solo su bracci non copiati dal bersaglio, prossima prova con Y-alt e Y-H1, proroga legata all'esito (consiglio: sì).
2. δ = 1 h con la fonte nuova (la precisione dell'area vera), sensibilità solo con δ = 0,5 (consiglio: sì).
3. A3 tolta, con il fatto dichiarato nell'articolo: lo Xenium in quadratini copre i soli geni del pannello, il 2% degli UMI di Visium HD, e la catena standard non ci gira (consiglio: sì; alternativa: secondo tentativo di registrare X2, almeno 2 sessioni, esito incerto).

[^fatti]: docs/riesame/FATTI_C.md
[^r4]: docs/sonda/REGOLA_S4.md
[^sim]: tools/sonda/S4_sim_arms.py
