
Commit dei cancelli: `fc79131` («F3.1: perimetro, esplorazioni, misure, cancelli»). Ok di Luca sulle previsioni, sulle soglie (intervalli corretti del vero), sulla regola R-1 e sulle controprove, e sul calcolo (circa 5 ore con 6 processi): 2026-11-12, «Ok, come proposte (Consigliato)». Margine minimo per M-4: 2026-11-12, «Margine minimo di un punto percentuale (Consigliato)» ([F3.1: margine minimo di un punto percentuale per la quota di singoletti di RCTD (DEC-086)](../docs/decisions/DEC-086.md)).

### Cancelli (al massimo 5)

| id | controllo | tolleranza | comando |
|---|---|---|---|
| G-1 | **batteria di prove** `tools/F3.1_test.py`: (a) esportazione simulata: per ogni quadratino UMI = origine + fondo; esportazione vera: conteggi uguali alla matrice grezza su 1 000 quadratini estratti; (b) quadratini in tessuto e UMI medio per quadratino uguali a `S2.5_real_areas.csv` sulle 30 aree vere; (c) somma per cellula della catena uguale a `bin_to_cell()` di bin2cell su un'area per tessuto; (d) misure M-1...M-6 e verità su una griglia costruita con risposte note; (e) mutanti: nuclei spostati (purezza vera < 0.1), tipi rimescolati (concordanza al livello del caso) | identità; (b) entro 1e-9 | `.venv-b2c/bin/python tools/F3.1_test.py` |
| G-2 | **invarianti per area**: UMI assegnati + non assegnati = totale; blocchi = somma dei 16 figli; blocchi passati a RCTD = blocchi con almeno 8 figli in tessuto; riferimento con al più 3 000 e almeno 25 cellule per tipo | identità | `tools/F3.1_check.py` |
| G-3 | **riproducibilità**: un'area per tessuto, simulata e vera, rifatta in un'altra cartella con un altro numero di processi | JSON identici senza i tempi | `tools/F3.1_repro.sh` |
| G-4 | **integrità**: md5 dei file letti prima e dopo (matrici e posizioni delle 5 sezioni, maschere e nuclei di R2 e R3b, conteggi ed etichette di R6, livelli di S1.6, profili di S2.2, tabelle di S2.5 e R7) | identità | `tools/F3.1_md5.sh before / after` |
| G-5 | **revisione avversariale**: sotto-agenti con codice proprio ricalcolano M-1...M-6 su un'area per tessuto per lato dai file esportati e dalle uscite dei metodi, gli intervalli e gli esiti dalle tabelle per area, la regola R-1 | quantità deterministiche entro 1e-9 relativo; ogni discrepanza è un dubbio | sotto-agenti |

### Regola di decisione dei dati (scritta prima del calcolo)

**R-1, che cosa sistemare per primo** (versione del 2026-11-12, dopo la fase A e prima del simulato; sostituisce la tabella della bozza `fc79131`, che legava ogni misura a un solo limite). Ogni giudizio smentito (misura × tessuto, cuore escluso; solo le misure che passano CP-1 in quel tessuto, per M-3 e M-5) si attribuisce così, nell'ordine:

1. **profili presi dai nuclei con la corona** (contaminazione contata due volte, E-3): se il braccio di CP-3 sposta la misura verso il vero di almeno metà dello scarto (valore del tessuto, simulato principale contro vero). Solo intestino, tumore, linfonodo, cervello (dove i profili vengono dai nuclei).
2. **densità dei preset** ([DEC-069](../docs/decisions/DEC-069.md), S2.2 punto 4): M-1 e M-2, se lo scarto ha lo stesso segno dello scarto fra cellule simulate e nuclei veri per area, e quel rapporto è fuori da 0.9-1.1.
3. **tipi disposti senza guardare la posizione e semina uniforme** (S2.2 punto 2; limiti della guida d'uso): M-4 e M-6, se lo scarto è per difetto.
4. altrimenti **non attribuito**: diventa un dubbio.

Il limite da sistemare per primo è quello con più giudizi attribuiti; a pari numero, quello con lo scarto medio più grande in unità della semi-ampiezza dell'intervallo del vero. Il risultato si registra come decisione tecnica.

### Previsioni

Unità indipendente: l'area (5 per tessuto e per lato). Intervallo di rilevanza = intervallo corretto del valore vero del tessuto (fonte (2) del protocollo: variabilità fra aree vere; costruzione di [DEC-060](../docs/decisions/DEC-060.md) e [DEC-061](../docs/decisions/DEC-061.md), B = 2 000, seme 20261115; tabella `results/F3.1/F3.1_intervals.csv`). Intervallo del simulato: stessa costruzione con il fattore per poche unità del lato simulato (1.632). Esito per tessuto: confermata, smentita, indeterminata.

**Valore atteso**: valore vero del tessuto più lo scarto simulato meno vero visto sull'area dell'esplorazione; **2 errori standard tipici del simulato**: semi-ampiezza del bootstrap del vero per 1.632, cioè supponendo che il simulato vari fra aree quanto il vero. La semi-ampiezza dell'intervallo del vero è 1.0-2.4 volte questi 2 errori standard: una conferma è possibile solo se il simulato varia fra aree meno del vero (in S2.5 accadeva); le smentite con scarti grandi restano nette.

**Margine minimo per M-4** (deciso da Luca il 2026-11-12, [DEC-086](../docs/decisions/DEC-086.md)): dove nel vero RCTD chiama singoletto quasi ogni blocco (intestino, linfonodo, cervello, cuore) l'intervallo del vero è largo meno di un millesimo; proposta: intervallo di rilevanza almeno [vero − 0.01, vero + 0.01] (un punto percentuale di blocchi chiamati doppietto). Con il margine, l'esito atteso di M-4 diventa: smentita nell'intestino, nel tumore e nello stroma; confermata o indeterminata nel linfonodo, nel cervello e nel cuore.

| misura | tessuto | vero | intervallo di rilevanza | atteso simulato ± 2 ES | esito atteso (senza margine minimo) |
|---|---|---|---|---|---|
| M-1 | A1 | 0.676 | 0.667-0.684 | 0.753 ± 0.007 | smentita |
| M-1 | A2 | 0.744 | 0.717-0.776 | 0.784 ± 0.029 | indeterminata |
| M-1 | A3 | 0.586 | 0.389-0.713 | 0.400 ± 0.083 | indeterminata |
| M-1 | A4 | 0.984 | 0.969-0.999 | 0.990 ± 0.010 | indeterminata |
| M-1 | A5 | 0.387 | 0.326-0.453 | 0.281 ± 0.057 | indeterminata |
| M-1 | A6 | 0.177 | 0.139-0.211 | 0.196 ± 0.036 | indeterminata |
| M-2 | A1 | 470 | 414-526 | 504 ± 49 | indeterminata |
| M-2 | A2 | 1136 | 923-1352 | 1084 ± 210 | indeterminata |
| M-2 | A3 | 669 | 403-837 | 631 ± 126 | confermata |
| M-2 | A4 | 102 | 56-149 | 130 ± 25 | indeterminata |
| M-2 | A5 | 4514 | 1534-6564 | 4782 ± 1062 | confermata |
| M-2 | A6 | 150 | 80-226 | 64 ± 73 | indeterminata |
| M-3 | A1 | 0.241 | 0.226-0.255 | 0.374 ± 0.013 | smentita |
| M-3 | A2 | 0.045 | 0.035-0.057 | 0.144 ± 0.011 | smentita |
| M-3 | A3 | 0.222 | 0.118-0.334 | 0.346 ± 0.063 | indeterminata |
| M-3 | A4 | 0.672 | 0.630-0.708 | 0.790 ± 0.021 | smentita |
| M-3 | A5 | 0.232 | 0.179-0.285 | 0.294 ± 0.022 | indeterminata |
| M-3 | A6 | 0.164 | 0.136-0.192 | 0.199 ± 0.028 | indeterminata |
| M-4 | A1 | 1.000 | 1.000-1.000 | 0.973 ± 0.000 | smentita |
| M-4 | A2 | 0.954 | 0.940-0.968 | 0.816 ± 0.014 | smentita |
| M-4 | A3 | 0.985 | 0.971-0.999 | 0.840 ± 0.007 | smentita |
| M-4 | A4 | 0.999 | 0.998-1.000 | 0.992 ± 0.001 | smentita |
| M-4 | A5 | 1.000 | 0.999-1.000 | 1.000 ± 0.000 | smentita |
| M-4 | A6 | 1.000 | 1.000-1.000 | 1.000 ± 0.000 | confermata |
| M-5 | A1 | 0.864 | 0.843-0.882 | 0.825 ± 0.017 | smentita |
| M-5 | A2 | 0.990 | 0.987-0.993 | 0.922 ± 0.003 | smentita |
| M-5 | A3 | 0.890 | 0.835-0.941 | 0.788 ± 0.027 | smentita |
| M-5 | A4 | 0.602 | 0.459-0.699 | 0.374 ± 0.076 | smentita |
| M-5 | A5 | 0.904 | 0.874-0.941 | 0.950 ± 0.030 | indeterminata |
| M-5 | A6 | 0.809 | 0.783-0.832 | 0.808 ± 0.025 | indeterminata |
| M-6 | A1 | 0.909 | 0.899-0.919 | 0.825 ± 0.009 | smentita |
| M-6 | A2 | 0.936 | 0.920-0.951 | 0.790 ± 0.015 | smentita |
| M-6 | A3 | 0.886 | 0.846-0.929 | 0.621 ± 0.021 | smentita |
| M-6 | A4 | 0.472 | 0.416-0.520 | 0.354 ± 0.033 | smentita |
| M-6 | A5 | 0.775 | 0.755-0.796 | 0.811 ± 0.019 | indeterminata |
| M-6 | A6 | 0.750 | 0.685-0.816 | 0.649 ± 0.066 | indeterminata |

In sintesi (P-1...P-6 = M-1...M-6, giudicate per tessuto):

- **P-1, UMI assegnati**: smentita per eccesso nell'intestino (più cellule nel preset); indeterminata negli altri.
- **P-2, UMI per cellula**: confermata o indeterminata ovunque (intervalli larghi).
- **P-3, marcatori estranei**: smentita per eccesso nell'intestino, nel tumore e nel linfonodo; indeterminata nello stroma, nel cervello e nel cuore.
- **P-4, singoletti**: smentita per difetto nell'intestino, nel tumore e nello stroma (con il margine minimo).
- **P-5, concordanza con il nucleo**: smentita per difetto nell'intestino, nel tumore, nello stroma e nel linfonodo; indeterminata nel cervello e nel cuore.
- **P-6, coerenza spaziale**: smentita per difetto nell'intestino, nel tumore, nello stroma e nel linfonodo; indeterminata nel cervello; nel cuore indeterminata (vedi CP-4).

### Controprove

- **CP-1, le misure vedono gli errori di assegnazione?** Nuclei spostati di 20 µm in diagonale (bin2cell) e tipo del nucleo preso 24 µm più in là (RCTD), sui due lati: M-3 deve salire e M-5 scendere di più di 2 errori standard in ogni tessuto. Dove non succede, quella misura in quel tessuto non dice nulla sull'assegnazione e il suo esito non entra nella regola R-1. **Lato vero, fase A: passa in tutti e sei i tessuti per M-3 e M-5** (`results/F3.1/F3.1_CP1_real.csv`; la più debole: M-5 nel cuore, 0.093 ± 0.037).
- **CP-2, le misure senza verità seguono la verità?** Sulle 60 esecuzioni simulate: correlazione fra aree della quota di marcatori estranei (M-3) con la quota vera di UMI da cellule di altro tipo, e della concordanza con il nucleo (M-5) con l'accuratezza vera di RCTD; attesa positiva in almeno 4 tessuti su 6.
- **CP-3, la contaminazione contata due volte.** Braccio con i profili presi dal solo nucleo (stessi sottotipi, stessi trascritti per cellula, stesso fattore di profondità; cambia solo la composizione in geni), 30 aree × 1 seme: se l'ipotesi di E-3 è giusta, M-3 scende verso il vero di almeno metà dello scarto in almeno 3 dei 4 tessuti con profili dai nuclei (intestino, tumore, linfonodo, cervello; nello stroma i profili vengono dall'atlante, il cuore è il controllo negativo). Profili: stessi file di S2.2 con le somme per sottotipo ricalcolate sui conteggi del solo nucleo degli stessi nuclei (`tools/F3.1_profiles_nuc.py`); trascritti per cellula, dispersioni e fattore di profondità invariati. L'UMI medio per quadratino si riporta come controllo.
- **CP-4, il cuore come controllo negativo.** Cellule isotrope nel simulato, fibre allungate nel vero: M-6 smentita per difetto nel cuore. Se non lo è, la misura non vede un limite noto. Nota dopo la fase A: l'esplorazione fa attendere un esito indeterminato (0.65 ± 0.07 contro 0.69-0.82), quindi questa controprova può fallire per mancanza di precisione, non solo di sensibilità; lo si dirà.

