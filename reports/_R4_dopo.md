# Dopo il risultato

## Esiti come scritti

- G-1 PASS: `tools/R4_test.py` 45/45 (32 test e 13 mutanti rilevati; i 3 mutanti pre-registrati più 10 aggiunti dopo la revisione). G-2 PASS: 238 file con md5 identico prima e dopo (i 5 file dei nomi dei geni aggiunti dopo la revisione, F14). G-3 PASS 45/45 (`tools/R4_check.py`: medie per bin a 2 e 8 µm uguali ai riassunti 10x a meno di 2e-13 nelle 5 sezioni; 0 bin a 8 µm diversi dalla somma dei loro 16 bin a 2 µm; impronte identiche a R6 in 30/30 ROI). G-4 PASS: 14 tabelle identiche fra due esecuzioni complete (stadi, analisi e post hoc). G-5: due revisori, sotto.
- P-1: confermata in intestino, stroma, cervello, cuore; smentita in tumore del colon (0.872) e linfonodo (0.801), nella direzione opposta all'attesa (nucleo meno denso). Regola R-1: uniforme in intestino e stroma, due densità negli altri quattro.
- CP-1 non superata in intestino (1.166), tumore (1.091), stroma (1.259), cervello (1.052): il rapporto sui bin raggruppati contiene un effetto di compartimento. Spiegata da PH-1 (sotto).
- P-2: confermata in cervello e cuore; indeterminata in intestino e stroma; smentita in tumore del colon (-0.026, zeri in difetto); non valutabile nel linfonodo (3 geni). CP-2 superata nei 5 tessuti valutabili.
- P-3: confermata in tutti (0 cellule sopra 50 000; massimo 28 126 nel cervello). Insieme delle cellule: tutte le cellule tenute, anche al bordo (scelta conservativa: le impronte al bordo sono troncate; revisione F05).
- P-4: confermata in cervello (3.76) e cuore (1.67); sopra 1.25 in intestino (1.83) e stroma (2.47); a cavallo in tumore del colon (1.449 [1.219, 1.676]) e linfonodo (1.344 [1.214, 1.455]). Regola R-3: un valore per tipo in quattro tessuti, uno solo in tumore e linfonodo (decisione 1 a Luca).

## Scelte decise dai dati e deviazioni dichiarate

- Riferimento del rapporto nucleo/citoplasma per lo Step 2: il rapporto pesato per cellula di PH-1 (nullo 0.98-1.00 in tutti i tessuti), con il rapporto pre-registrato come sensibilità: [Rapporto nucleo/citoplasma pesato per cellula come riferimento (DEC-034)](../docs/decisions/DEC-034.md). Le classi della regola R-1 non cambiano.
- P-4: la prima esecuzione calcolava la quota del 5 % sulle sole cellule tipizzate; la pre-registrazione dice «dei nuclei». Corretto dopo la revisione (RV-1, F03); la versione precedente è in `results/R4/R4_sensitivity.csv` (tumore del colon 3.48, sopra 1.25).
- Il comando di G-3 (`tools/R4_check.py`) mancava alla prima esecuzione: il controllo era calcolato dagli stadi; aggiunto (F04).

## Analisi post hoc (dichiarate dopo i risultati)

- PH-1 (CP-1): rapporto nucleo/corona pesato per cellula (media per cellula della densità nel nucleo diviso media per cellula della densità nella corona) e mediana dei rapporti per cellula, ciascuno diviso per lo stesso con la maschera spostata. Il nullo pesato per cellula vale 0.98-1.00; le direzioni sono quelle di P-1 in tutti i tessuti (intestino 0.986, tumore 0.847, stroma 0.892, linfonodo 0.788, cervello 1.383, cuore 1.784). Il nome «appaiato» della prima versione era sbagliato (F08).
- PH-2 (P-2): binomiale negativa con fattore di taglia per cellula, sugli stessi geni e sullo stesso dropout di P-2 e CP-2. Zeri in eccesso ≈ 0 in tumore, stroma, cervello, cuore; +0.04 nell'intestino, dai marcatori di sottotipi (Anpep, Apoa1, Fabp6, Defa21, Lyz1, Reg3b; Igkc, Igha, Jchain fra le immunitarie; stromali 0.0001; `results/R4/R4_PH2_A1_genes.csv`). Il revisore del codice ha mostrato per simulazione che la statistica pre-registrata è spostata verso il negativo quando le cellule hanno quantità diverse di trascritti (-0.034...-0.005 senza alcun dropout) e che lo stimatore con fattore di taglia non lo è (F06). Bootstrap sulle sole ROI: le classi non si confrontano con i verdetti di P-2.
- PH-3 (linfonodo): soglia 0.5 UMI per cellula. 13 geni, ma il dropout di controllo non ha più geni sufficienti (7): non valutabile.
- PH-4 (D-5): alternanza della cattura per parità di riga e colonna della griglia a 2 µm, sezione per sezione e ROI per ROI (media raggruppata).
- Sensibilità aggiunte dopo la revisione: S-4 (bin a ≥ 8 µm dal bordo della ROI, F15) e intervallo ± 2 e ± t(4) errori standard per P-1 (F09). Cambiano classe solo il tumore del colon con ± t(4) (indeterminata) e lo stroma con ± 2 e ± t(4) (indeterminata, uniforme per difetto): nessuna decisione cambia.

## Revisione avversariale

Due revisori indipendenti, in sola lettura sul repo. Il revisore numerico ha ricalcolato da zero (codice proprio, senza le funzioni del progetto) 935 numeri: identici tutti salvo P-4 nel tumore del colon (RV-1, deviazione, corretta), le quote «altro nucleo» e «acellulare» (definizione diversa, somma identica, RV-3) e i rapporti di parità (peso dei quadranti, scarti ≤ 2.2e-5, ora media raggruppata, RV-2). Il revisore del codice ha trovato 2 bug (lo spettro non guardava la frequenza di Nyquist, F01; lo script di riproduzione non includeva il post hoc, F02), 3 deviazioni (F03-F05), 1 discrepanza (F07) e 11 note; i suoi limiti dichiarati: strati ricalcolati da zero su 5 ROI su 30, simulazione degli zeri provvisoria (300 geni, 3 repliche), spettri delle ROI non ricalcolati, G-4 non osservato. Correzioni: F01 (Nyquist incluso, ricerca fuori dagli assi per il reticolo diagonale, F16), F02, F03, F04, F07, F08, F10 (10 mutanti nuovi, tutti rilevati), F14, F15; dopo le correzioni l'intero calcolo è stato rifatto due volte (G-4). I numeri corretti coincidono con quelli dei revisori: P-4 tumore 1.449 [1.219, 1.676]; sfondo dell'intestino a Nyquist 24.1.

## Errori di pre-registrazione

- La statistica di P-2 (binomiale negativa con una sola media per gene) non vale 0 senza dropout quando le cellule hanno quantità diverse di trascritti: la soglia ±0.02 confonde il difetto di zeri dovuto alle taglie con l'assenza di dropout. P-2 resta giudicata come scritta; la lettura pertinente è PH-2.
- La regola di D-5 (un solo picco, il più alto) non poteva separare un reticolo diagonale dall'alternanza di righe e colonne; dopo la revisione si riporta anche il picco fuori dagli assi.
