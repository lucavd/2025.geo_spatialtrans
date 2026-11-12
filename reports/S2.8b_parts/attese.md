Commit dei cancelli: `4d2c637` (perimetro, scelte di Luca all'apertura, cancelli e controprove). Ok di Luca sull'assenza di previsioni, sui cancelli, sulle controprove e sul calcolo di circa 2 ore-processo: 2026-11-12, «Approvo (Recommended)».

### Cancelli (al massimo 5)

T = commit provato (documenti, guida resa e script dei cancelli; dopo T cambiano solo `docs/`, `reports/`, `CHANGELOG.md`).

| id | controllo | tolleranza | comando |
|---|---|---|---|
| G-1 | **I documenti nominano solo cose che esistono, e la guida copre l'ingresso.** (a) Ogni percorso del repository fra apici inversi e ogni funzione R scritta come `nome()` o chiamata nei blocchi di codice R di `README.md`, `CLAUDE.md`, della guida d'uso e della guida per Daniele (con i file che include) esiste nell'albero di T (come in S1.6b). (b) Ogni argomento di `simulate_visium_hd_pieces()` (letto dal sorgente di T) compare fra apici inversi nella guida d'uso; ogni funzione dei file `R/05b*.R` marcata `@export` compare nella guida d'uso o nel README. **Mutanti** (3), su copie dei documenti: un percorso inventato, una funzione inventata, la guida senza la riga di un argomento: rilevati tutti. | 0 mancanti; mutanti 3/3 | `python3 -B tools/S2.8b_doc_check.py --tree T`; `--mutants` |
| G-2 | **Gli esempi girano come sono scritti.** In un clone nuovo di T da GitHub, libreria dal lockfile, librerie dell'utente escluse: (a) l'avvio rapido dello Step 1 del README, due esecuzioni: 0 ERROR nei controlli e md5 di `cell_df` = `9f361ea9e1a56d8f46eef8a0a7fba842` (come in S1.6b) in entrambe; (b) l'avvio rapido dello Step 2 del README, estratto senza modifiche, due esecuzioni: md5 di `counts_2um` e `counts_8um` uguali fra le due; (c) `quarto render` della guida d'uso due volte: nessun errore, le impronte stampate dalla guida uguali fra le due rese, e i controlli che la guida stessa esegue (stesso risultato con pezzi diversi e con più processi; ripresa che riusa i file; ripresa rifiutata con argomenti diversi) tutti `TRUE`; nessun file versionato modificato | esatta | `bash R/testing/run_S2.8b.sh clone T`, `quickstart`, `guide` |
| G-3 | **Le batterie passano nel clone nuovo.** Nel clone di T: `bash R/testing/run_S1.6.sh test` stampa `GATE G-1 PASS` con i soli quattro FAIL noti (come S1.6b); `R/testing/test_S2.1.R` ... `test_S2.8.R` con le stesse asserzioni di S2.8 (21, 36, 39, 37, 40, 33, 78, 46), 0 FAIL | esatta | `bash R/testing/run_S2.8b.sh battery` |
| G-4 | **L'unione non perde e non cambia nulla del simulatore.** A fine step: `git ls-remote origin dev` = `git ls-remote origin step2-visium-grid`; il vecchio `dev` (`9250749`) è antenato (avanzamento); `main` su GitHub = `ec29fd6` come all'apertura; `git diff 48a5f18 dev -- R/ inst/ renv.lock` vuoto salvo `R/testing/run_S2.8b.sh`; `DESCRIPTION` diverso da `48a5f18` solo nei campi `Version`, `Date`, `Description`; dopo T solo `docs/`, `reports/`, `CHANGELOG.md` | esatta | `bash R/testing/run_S2.8b.sh merge` |
| G-5 | **I numeri dei documenti tornano e il bundle è in regola.** Ogni numero della guida d'uso, del README e della guida per Daniele che non è prodotto dall'esecuzione della guida stessa ha una riga nella tabella «Da dove vengono i numeri» (file, riga o colonna, valore); due revisori indipendenti rileggono i numeri dai file e i documenti contro il codice: le discrepanze sono dubbi, chiusi prima del verdetto. Controllo OKF: 0 documenti del bundle in errore | 0 discrepanze aperte | revisione avversariale; `python3 -B tools/okf.py check --bundle docs` |

### Previsioni

Nessuna previsione con intervallo di rilevanza: lo step non cambia il modello generativo né un valore biologico (G-4: `R/`, `inst/`, `renv.lock` identici a S2.8), quindi non c'è una quantità simulata da giudicare. La correttezza dei documenti e dell'unione è binaria e sta nei cancelli. Regola applicata in forma leggera, come in S1.6b; unità indipendente, B e seme non si applicano. I tempi degli esempi della guida sono descrittivi.

### Controprove

- **CP-1, il controllo di copertura deve bocciare i documenti di oggi.** `tools/S2.8b_doc_check.py --coverage` sul `README.md` e sul `CLAUDE.md` di `48a5f18` (prima di questo step) deve trovare scoperti argomenti di `simulate_visium_hd_pieces()` e funzioni esportate dello Step 2. Se li promuovesse, G-1 (b) non distinguerebbe un documento vecchio da uno giusto.
- **CP-2, la guida si esegue davvero.** Una copia della guida d'uso in cui un esempio chiama `simulate_visium_hd_pieces()` con `piece_um = 500` (non multiplo di 56) deve far fallire `quarto render`. Se la resa passasse, G-2 (c) non proverebbe che gli esempi girano.

### Stima dei tempi (dell'agente)

Dai tempi di S1.6b e S2.8: clone e libreria circa 1 minuto; avvii rapidi circa 3 minuti; due rese della guida circa 10 minuti (esempi su aree da 0.16 a 0.25 mm², da 30 a 60 s l'uno); batteria dello Step 1 circa 36 minuti; batterie dello Step 2 circa 25 minuti con 8 processi. In tutto circa 2 ore-processo, circa un'ora di orologio, al più 9 processi insieme.
