---
name: progetto-master
description: Protocollo master dei progetti di ricerca destinati a un articolo. Caricare all'inizio di ogni sessione di progetto e per avviare, riprendere, riesaminare o chiudere un progetto.
---

# Scopo

Portare un progetto di ricerca a un articolo rigoroso con un messaggio semplice, spendendo pochi token e poche ore del PI. Se l'obiettivo non si può raggiungere, scoprirlo presto e dirlo. Il PI è Luca Vedovelli. Le regole proprie di un progetto stanno nel suo `CLAUDE.md`. Se una skill di progetto contraddice questa sul flusso di lavoro, vale questa.

# Le parole

Queste nove parole hanno un solo significato in tutti i progetti. Con il PI si usano senza rispiegarle. L'esempio è un programma che assegna a ogni cellula il suo tipo cellulare.

| parola | significato | nell'esempio |
|---|---|---|
| misura principale | Il numero che dice se l'articolo regge. Si calcola sul risultato finale, non sui singoli pezzi. | La quota di cellule a cui il programma dà lo stesso tipo dato dal patologo. |
| pavimento | Il valore della misura senza il nostro lavoro: con la soluzione banale o con il miglior strumento esistente. | Lo strumento più usato arriva al 71%. |
| tetto | Il massimo sensato della misura. Oltre non si va, perché nemmeno il riferimento concorda con se stesso. | Due patologi fra loro concordano sull'88% delle cellule. |
| soglia | Il valore, fra pavimento e tetto, che l'articolo deve raggiungere per sostenere il messaggio. Ha una fonte. | 80%: sotto questo valore il lettore non cambierebbe strumento. |
| metro | Il codice che calcola la misura principale. | Lo script che confronta i tipi assegnati con quelli del patologo. |
| catena | Il percorso completo dal dato grezzo alle figure dell'articolo. | Dai conteggi grezzi alle tre figure del contratto. |
| bilancio dell'errore | La tabella che divide per causa la distanza fra la misura di oggi e il tetto. Dice dove lavorare. | Siamo al 74%. Dei 14 punti che mancano, 8 vengono dalla normalizzazione, 4 dal riferimento, 2 dal resto. |
| colpo | Una sessione che non migliora la misura oltre la sua incertezza. | Una modifica attesa a +4 punti che ne dà +0,5, con un'incertezza di 1 punto. |
| quadro dell'articolo | La tabella, calcolata dallo strumento, con il messaggio e lo stato di ogni affermazione: regge, non regge, incerta, non misurata. | La prima affermazione regge (74% contro 70%), la seconda non regge, la terza non è misurata. |

# Il metodo in sette righe

1. Prima del lavoro si firma un contratto: un messaggio, le sue affermazioni, una misura principale con pavimento, tetto e soglia, un budget.
2. Il primo codice è il metro. Il secondo è la catena più banale che gira da capo a fondo: dà il pavimento e tutte le figure in forma grezza.
3. Un pezzo entra nella catena solo se la migliora oltre l'incertezza. Vale anche per il codice che esiste già.
4. Ogni sessione parte dal quadro dell'articolo, serve una sua affermazione e lavora sulla voce più grande del bilancio dell'errore.
5. Il rigore si spende sul metro e sulla corsa di conferma. Il resto è usa e getta finché non entra nell'articolo.
6. Due colpi chiudono un'ipotesi, non il progetto. L'obiettivo si chiude solo per una delle condizioni scritte in «Quando fermarsi».
7. A ogni cancello «chiudere» è fra le opzioni. Raccomandala quando quelle condizioni ci sono.

# Perché funziona così

In tre progetti del PI l'esito ha seguito il criterio di successo e il momento in cui lo si è misurato.

- MMM aveva un numero con una fonte, l'errore di I tipo dentro l'intervallo di Bradley. È arrivato alla rivista in cinque settimane.
- simulomicsr aveva «nessun errore si accetta». È durato sei mesi e mezzo. Il tetto, misurato solo dopo, era 15 studi su 21 fra due gruppi di curatori umani.
- geo_spatialtrans aveva «biologicamente solido», senza una misura, per 43 sessioni. Poi due sessioni con metro, tetto e versione banale hanno mostrato che il generatore costruito non batteva il banale: 7,83 contro 7,77, con il tetto a 2,23.

Quindi criterio, tetto e pavimento si misurano nelle prime sessioni, non alla fine.

# Il contratto

Il contratto dice che cosa vuol dire «riuscito» prima di cominciare. Sta in una pagina, in `docs/PAPER.md`. I suoi numeri stanno anche in `progetto.json`, dove li legge lo strumento. Contiene:

- il messaggio in una frase di al massimo 25 parole, il lettore, la rivista;
- la misura principale, con l'unità indipendente su cui si calcola l'incertezza (paziente, studio, replica);
- pavimento e tetto, misurati da noi e non supposti;
- la soglia, con la sua fonte;
- da 3 a 5 affermazioni, ognuna con la sua figura e con il numero che le serve per reggere;
- tre esiti scritti prima: A, B, e C, che è «non c'è articolo»;
- i limiti che accettiamo già adesso;
- il budget di ogni fase (sessioni, calcolo, ore del PI), proposto con il suo motivo e firmato dal PI;
- i dati tenuti da parte, che nessuno guarda fino alla conferma.

Il tetto si misura in uno di questi modi: l'accordo fra due esperti, un dato vero contro un altro dato vero, l'errore Monte Carlo, la versione perfetta di un componente.

Regole sul criterio:

- Un criterio senza misura o senza fonte si contesta prima di lavorare, anche se viene dal PI («nessun errore», «solido», «il migliore»). Proponi tu la versione misurabile.
- Se la soglia sta oltre il tetto, il criterio non si può raggiungere: riposiziona. Se pavimento e tetto non si distinguono, non c'è spazio per un articolo: non partire.
- In uno studio di simulazione il contratto è il protocollo ADEMP (Morris et al., 2019). La misura è una misura di prestazione, per esempio l'errore di I tipo. Il riferimento è il valore nominale, la soglia è il margine di rilevanza, e l'errore Monte Carlo dice quali differenze si riescono a vedere.
- In uno studio applicato, cioè su dati reali con una domanda clinica o biologica, la misura è la precisione della stima principale e la risposta dei controlli. Non è mai la direzione dell'effetto, altrimenti si lavora finché esce il risultato sperato. Il tetto è la precisione che i dati permettono, e si calcola prima. Fino alla conferma si lavora in cieco sull'esito.

# Fasi e cancelli

| fase | durata | che cosa produce | cancello del PI |
|---|---|---|---|
| Scommessa | 1 o 2 sessioni | Contratto, ricerca dell'esistente, obiezioni di un revisore simulato. Nessun codice. | G0: firma del contratto |
| Sonda | 3 sessioni, fino a 5 con dati pesanti | Metro con i suoi controlli, tetto, catena banale con tutte le figure in forma grezza (è il pavimento), primo bilancio dell'errore. | G1: partire, riposizionare o chiudere |
| Costruzione | firmata a G1 | La misura sale verso la soglia, una voce del bilancio per volta. | G2: soglia raggiunta sui dati di sviluppo, oppure budget finito |
| Conferma | 1 o 2 sessioni | Catena congelata, una sola pre-registrazione, corsa completa sui dati tenuti da parte. | G3: i numeri dell'articolo |
| Scrittura | firmata a G3 | Articolo scritto a partire dal contratto, revisione finale. | invio |

Il budget non è uguale per tutti i progetti. Scommessa e sonda restano corte comunque, perché rispondono a una domanda sola: c'è spazio per un articolo? Costruzione e scrittura crescono con il progetto. Per la costruzione proponi due sessioni per ogni voce grande del bilancio, perché ogni ipotesi ha due colpi. Il PI firma il numero al cancello, e lo strumento lo registra con le sue parole (`cancello --fase FASE --budget N --consenso TESTO`). Lo strumento non lascia entrare in costruzione senza pavimento, tetto, soglia e una prima misura di ogni affermazione.

A budget finito il PI può prorogare. La richiesta porta quattro cose: dove è arrivata la misura, quante sessioni in più, il valore atteso, che cosa si fa se non lo si raggiunge (`proroga`). La seconda proroga della stessa fase chiede prima un riesame: decide sempre il PI, ma con i fatti sul tavolo.

Il rigore pieno (test completi, pre-registrazione, verifica indipendente) si spende in due punti: sul metro e sulla corsa di conferma. Altrove il codice è usa e getta, e si testa solo ciò che resta nella catena. Nessun risultato si chiama «finale» prima della corsa di conferma.

A ogni cancello porta al massimo 3 opzioni, ognuna con costo e guadagno atteso. «Chiudere» è sempre una delle tre. Prima di G1 e di G3 un sotto-agente a contesto pulito ricalcola i numeri dai dati, perché chi ha fatto il lavoro non lo giudica da solo.

Il PI decide ai cancelli e su ciò che spetta solo a lui: messaggio, perimetro, rivista, autori, budget, accessi, azioni irreversibili. Fra un cancello e l'altro le sessioni non aspettano la sua approvazione: una scelta tecnica si fissa con la sua fonte, si scrive nel messaggio e vale finché il PI non la cambia.

Ogni aggiunta al perimetro ha un prezzo, anche se la chiede il PI: quale affermazione serve, quante sessioni e quanto calcolo costa, che cosa cambia nell'articolo nel caso migliore. Se cambia una frase, dillo e sconsigliala. Una corsa completa si fa solo per ciò che entra nell'articolo: prima risponde una corsa ridotta.

# La sessione

Una sessione è una chat, serve una sola affermazione del contratto e prova una sola ipotesi.

Apertura:

1. Esegui `python3 tools/progetto.py apri`. Stampa una riga di stato, il quadro dell'articolo, le ultime misure e gli allarmi. Copia la riga di stato in testa alla tua prima risposta: `progetto-master | fase | sessione n di budget`. Per il PI è la prova che la skill è caricata. Se c'è un allarme, la sessione serve a risolverlo: non si lavora sopra un allarme.
2. Leggi `docs/PAPER.md` e `docs/STATO.md`. Non leggere altro, salvo ciò che `STATO.md` indica.
3. Dalla sonda in poi esegui `python3 tools/progetto.py check smoke`. Se fallisce, la sessione serve a ripararlo.
4. Scrivi in tre righe l'affermazione che la sessione serve, l'ipotesi (quale voce del bilancio attacchi e come), l'effetto atteso come numero, il costo. Se la sessione non serve nessuna affermazione, o l'effetto atteso è più piccolo dell'incertezza della misura, non si fa.

Chiusura:

1. Esegui `check smoke` (non nella scommessa, dove non c'è codice), poi `python3 tools/progetto.py chiudi --sessione ID --affermazione ID --ipotesi TESTO --atteso NUMERO`. Lo strumento legge dal file scritto dal check il valore di ogni affermazione. Stabilisce se la sessione avanza o è un colpo sull'affermazione servita, e avvisa se un'altra affermazione che reggeva non regge più. Rifiuta la chiusura se il check è fallito, se il codice è cambiato dopo il check, se manca l'affermazione, o se un documento supera il suo limite di lunghezza.
2. Se è un colpo, togli la modifica dalla catena.
3. Riscrivi `docs/STATO.md`. Poi commit e push.
4. Scrivi al PI il messaggio di chiusura, secondo la sezione «Scrivere per il PI». Non produrre un report.

`docs/STATO.md` è l'unica memoria fra una sessione e l'altra. Contiene: che cosa esiste, la tabella del bilancio dell'errore, le minacce alla validità (5 al massimo), le strade fallite con il motivo, la prossima mossa.

# Il metro e il controllo dell'insieme

Il difetto da evitare è questo: ogni pezzo passa i suoi test, ma l'insieme non serve allo scopo, e lo si scopre alla fine. Per questo a ogni sessione si ricalcola su tutta la catena il quadro dell'articolo: la misura principale e ogni altra affermazione. Un'affermazione che reggeva e non regge più è una regressione, e si ripara prima di ogni altro lavoro.

Il metro è il primo codice che si scrive, e l'unico con test completi da subito. Ha un controllo positivo, cioè un caso in cui la risposta è nota e deve uscire, e un controllo negativo, in cui non deve uscire nulla. Un numero prodotto da un metro senza controlli non si usa per decidere.

Se la misura principale è una media di più voci, leggi prima le voci: una media dominata dalle voci più grandi nasconde le differenze sulle altre. Il check scrive anche le voci, e il quadro le mostra.

`check smoke` gira a ogni apertura e a ogni chiusura, in circa dieci minuti, su una configurazione ridotta e sui dati di sviluppo. Esegue i test rimasti, la catena da capo a fondo e il metro, e scrive in `results/misura.json` la misura principale, le sue voci e il valore di ogni affermazione. `check full` è la versione completa e gira ai cancelli.

Non si modifica un test, una soglia o un elenco di esenzioni per far passare un controllo.

# Il bilancio dell'errore

Il bilancio decide la prossima sessione. Si stima in uno di questi modi, dal più solido:

1. Sostituzione con la versione perfetta: al posto di un componente metti la verità nota o il dato reale, e misuri quanto sale la misura.
2. Ablazione: al posto di un componente metti la versione banale, e misuri quanto scende.
3. Lettura degli errori: classifichi per causa un campione di almeno 50 errori.
4. Dove non c'è una misura: le obiezioni di un revisore, in ordine di gravità.

Regole:

- Lavora sulla voce più grande. Una voce sotto il 10% non riceve sessioni.
- Rifai il bilancio dopo ogni sessione che avanza e dopo due colpi.
- Non esiste un backlog. Ciò che non è nel bilancio e non minaccia la validità si cancella.
- Un componente che il bilancio non nomina resta com'è, anche se è grezzo.

# Cercare prima di costruire

Si cerca nella fase di scommessa, prima di ogni componente che costa più di un'ora, e a ogni blocco. Si cercano quattro cose:

- l'articolo già scritto;
- lo strumento che fa già il lavoro;
- il metro e la verità già pronti: benchmark, banche dati curate, dataset annotati, linee guida;
- il modo in cui altri hanno superato lo stesso ostacolo.

L'ordine delle fonti è: repository e progetti passati del PI; Zotero (`zotero-library`); CRAN, Bioconductor, PyPI, GitHub; PubMed, bioRxiv, arXiv; web. Servono almeno tre ricerche con parole diverse per tipo di fonte. Le ricerche larghe le fa un sotto-agente, che restituisce dieci righe con i link.

Ogni ricerca lascia una riga in `docs/ESISTENTE.md`: bisogno, candidati con link, decisione, motivo. «Costruire» è ammesso in due casi: non c'è un candidato, oppure il migliore è stato provato sul nostro caso e ha fallito. Il candidato migliore diventa il pavimento. Si cita solo ciò che è stato aperto, e prima di un contratto i lavori più vicini si leggono per intero.

# Quando ti blocchi

Il PI non è un oracolo. Prima di scrivergli sali sei gradini e annota l'esito di ognuno:

1. Scrivi che cosa ti aspettavi, che cosa hai osservato e tre ipotesi sul perché.
2. Distingui le ipotesi con un calcolo sui dati che hai.
3. Cerca in `STATO.md` se la risposta c'è già.
4. Cerca nei dati e nei lavori di altri.
5. Fai l'esperimento più piccolo in cui la verità è nota.
6. Cambia strada: uno strumento esistente, un'affermazione più stretta, un limite dichiarato.

Solo dopo scrivi al PI una domanda sola, in forma di decisione, con la tua raccomandazione e con la scelta che applichi se non risponde.

Non chiedere al PI: lavoro manuale oltre 15 minuti, un valore che si può stimare o citare, una scelta che i dati possono decidere, la conferma di una cosa già decisa. Se una scelta si decide con meno di un'ora di calcolo, calcola tutte le opzioni invece di chiedere.

# Quando fermarsi

Fermarsi ha tre livelli. Chiudere un'ipotesi o togliere un pezzo non chiude il progetto.

| che cosa si ferma | quando | che cosa succede |
|---|---|---|
| Un'ipotesi | Due colpi di fila sulla stessa ipotesi. | La strada è chiusa. Non si propone una terza prova, neanche più piccola. Si passa a un'altra voce del bilancio. |
| Un pezzo | Non batte la sua versione banale, oltre l'incertezza. | Si toglie, e resta la versione banale. Vale anche per il codice ereditato, qualunque lavoro sia costato. |
| L'obiettivo | Vale una di queste quattro condizioni. Fra pavimento e tetto non c'è spazio. L'articolo esiste già, o non ha un lettore. I dati per misurarlo non esistono. Il budget firmato è finito, proroga compresa, e la misura non si è avvicinata alla soglia. | Riesame. Poi decide il PI. |

Lo stato dei pezzi non decide dell'obiettivo. Un progetto condotto finora senza metro non ha ancora avuto la sua prova: riparte dal pavimento, non si chiude. Ogni regola di esito scritta prima di una prova dice a quale dei tre livelli si applica.

Gli allarmi (due colpi, tre colpi, budget finito, criterio irraggiungibile, regressione) li calcola `tools/progetto.py` dal registro, non il giudizio di chi lavora. Un allarme dice quale livello guardare: non è un verdetto.

Il riesame si apre con tre colpi di fila, con un allarme sul criterio, alla seconda proroga, o quando esce un lavoro con lo stesso messaggio. Occupa una sessione, non tocca il codice e si registra con `chiudi --riesame`:

1. Attacca prima il criterio, poi il metodo. La soglia è sotto il tetto? Il tetto è misurato? La letteratura usa un altro riferimento?
2. Scrivi i fatti: solo numeri riletti dai file.
3. Giudica i pezzi uno per uno contro la loro versione banale: tenere, correggere, sostituire, togliere. Per ognuno scrivi l'osservazione che ribalterebbe il giudizio.
4. Giudica l'obiettivo a parte, sulle quattro condizioni. I costi già sostenuti non contano, in nessuno dei due sensi: il lavoro fatto non è un motivo per continuare, e il lavoro fatto male non è un motivo per chiudere.
5. Due sotto-agenti che non conoscono la storia del progetto ricevono i soli fatti. Uno difende «continuare», l'altro «chiudere». Un terzo ricalcola pavimento, tetto e misura.
6. Dai un verdetto sull'obiettivo, uno fra quattro: continuare; restringere il messaggio; cambiare strada, cioè stesso obiettivo e catena nuova che riparte dal pavimento; chiudere. «Chiudere» cita la condizione che vale, con il suo numero. Una strada nuova comincia dalla prova più economica della sua ipotesi chiave, con la regola di esito scritta prima, e ha due colpi come ogni altra.

Se le condizioni ci sono, scrivilo: «Ci abbiamo provato con il metodo giusto. Con i nostri mezzi non arriviamo a un risultato decente. Consiglio di chiudere.» La chiusura lascia `docs/CHIUSURA.md`, in una pagina: che cosa abbiamo misurato, quale condizione vale, che cosa si riusa, quale vincolo riapre il progetto se cade.

# Scrivere per il PI

Il PI decide su ciò che gli scrivi, quindi deve capire l'essenziale senza aprire un file. Conosce la domanda scientifica e il contesto. Non è uno specialista dei metodi: lavora con dati omici ma non è un bioinformatico, conosce la statistica medica ma non è uno statistico puro. Non spiegargli che cos'è un intervallo di confidenza o l'RNA-seq. Spiegagli che cosa fa uno strumento, un algoritmo o uno stimatore poco comune, e perché qui conta.

Ogni messaggio al PI comincia dal quadro dell'articolo, in due righe: il messaggio, quante affermazioni reggono e che cosa è cambiato oggi. Poi risponde a cinque domande, in quest'ordine:

1. Che cosa volevamo sapere, e perché serve all'articolo.
2. Che cosa abbiamo fatto, in una o due frasi che capirebbe un collega di un altro laboratorio.
3. Che cosa è uscito: il numero con unità e incertezza, accanto ai suoi termini di confronto (pavimento, tetto, soglia, valore di prima).
4. Che cosa significa: la frase con «quindi» che porta dal numero alla conseguenza per il progetto.
5. Che cosa proponi, quanto costa, che cosa si perde scegliendo altro.

Regole:

- Ogni termine tecnico, sigla o nome di strumento si spiega in una riga la prima volta che compare in un messaggio: che cosa fa e perché conta qui. Le nove parole della tabella non si rispiegano.
- Un codice interno (di sessione, decisione, misura, campione) non sta mai da solo: prima il nome per esteso, poi il codice fra parentesi.
- I legami logici restano nel testo: perché, quindi, però, invece. Le frasi sono brevi ma non telegrafiche.
- Un concetto nuovo si aggancia a uno che il PI conosce: «il tetto è come l'accordo fra due patologi».
- Si usa un termine per ogni concetto, sempre lo stesso. I numeri si rileggono dai file.
- La lunghezza è quella che serve per decidere: da 8 a 15 righe a fine sessione, una pagina a un cancello. Le note di servizio (commit, controlli superati, file toccati) stanno nei file, non nel messaggio.
- Prima di inviare fai questa prova: dopo una sola lettura, il PI saprebbe spiegare la decisione a un collega?

Un esempio con numeri inventati. Questo è troppo asciutto:

> S4 chiusa. ARI 0,61, atteso 0,70: colpo. Voce maggiore del bilancio: normalizzazione, 38%. Prossima: SCTransform.

Questo va bene:

> Delle tre affermazioni dell'articolo ne reggono due, come ieri. Non regge ancora la principale: i gruppi di cellule trovati dal nostro metodo devono coincidere con i tipi annotati dai patologi. Oggi abbiamo lavorato su quella. L'accordo si misura con l'indice ARI (0 è un accordo casuale, 1 un accordo perfetto). Lo strumento più usato arriva a 0,55, due patologi fra loro a 0,82, e all'articolo serve 0,70. La modifica di oggi doveva portarci a 0,70, ma siamo rimasti a 0,61: è un colpo. Ho poi sostituito un pezzo per volta con la sua versione perfetta, per vedere dove si perde di più. Più di un terzo della distanza dal tetto viene dalla normalizzazione, cioè da come si rendono confrontabili cellule lette con profondità diverse. Quindi la prossima sessione lavora lì, con un metodo già pubblicato (SCTransform). Costa una sessione.

Regole di forma:

- Ai cancelli il riassunto è un «In breve» completo, nella forma approvata dal PI (`mmm-session`, file `esempi/in_breve_S5.md`): a cosa serviva, cosa abbiamo trovato, cosa devi decidere. Il messaggio di fine sessione è la stessa cosa in piccolo.
- I file per l'agente (`STATO.md`, `ESISTENTE.md`, commit) restano telegrafici. I documenti nuovi hanno il frontmatter OKF.
- Niente lineette come punteggiatura. Nei documenti niente etichette in grassetto in testa ai punti di un elenco: sono ammesse solo nei messaggi in chat. Codice, tabelle e figure sono in inglese.
- Per l'articolo si usano `manuscript-style-review`, `revisione-e-piano-implementazione` e `zotero-library`.

# Economia dei token

Il costo di una sessione cresce con il numero di chiamate moltiplicato per la lunghezza del contesto: ciò che entra nel contesto si ripaga a ogni chiamata successiva.

- Tieni le sessioni corte: due sessioni da 40 chiamate costano molto meno di una da 80. Oltre le 60 chiamate chiudi alla prima occasione: finisci il passo in corso, riscrivi `STATO.md`, fai commit e fermati. Delegare a sotto-agenti non sostituisce la chiusura.
- In apertura leggi meno di tre pagine. Ogni comando stampa al massimo 30 righe: il resto va su file e si interroga con grep.
- La lettura pesante (letteratura, log, codice altrui) la fa un sotto-agente: il suo contesto si butta e torna una sintesi di dieci righe.
- Non scrivere un report per sessione né un documento per decisione: bastano una riga in `STATO.md` e il messaggio al PI.
- I modelli locali fanno il lavoro meccanico su molti casi, con l'errore misurato su un campione etichettato. Non giudicano e non scrivono per il PI.
- Scommessa, cancelli e riesame usano il modello più capace. Le sessioni di costruzione hanno un'ipotesi chiusa e un check, quindi possono usare un modello più economico, se il PI lo decide.

# I file

| file | contenuto | limite di lunghezza |
|---|---|---|
| `CLAUDE.md` | Scheda del progetto: dove sta tutto, regole scientifiche, limiti di calcolo. La prima riga dice di caricare `progetto-master` e di fermarsi se manca. | 100 righe |
| `.claude/skills/progetto-master/SKILL.md` | Copia di questa skill, per le sessioni che non vedono quelle dell'account. Si sostituisce a mano quando la skill cambia. | nessuno |
| `docs/PAPER.md` | Il contratto e le decisioni del PI, una riga ciascuna. | 500 parole |
| `progetto.json` | I numeri del contratto, le affermazioni, i comandi del check, i cancelli e le proroghe. | nessuno |
| `docs/STATO.md` | La memoria fra le sessioni, con il bilancio dell'errore. | 800 parole |
| `results/misure.tsv`, `results/quadro.tsv` | Una riga per sessione chiusa, e il valore di ogni affermazione a ogni chiusura. Li scrive lo strumento. | nessuno |
| `docs/ESISTENTE.md` | Registro di ciò che esiste già. | 1 riga per voce |
| `tools/progetto.py` | Comandi `init`, `apri`, `quadro`, `check`, `chiudi`, `cancello`, `proroga`. Si copia da un progetto che lo ha. | nessuno |

# Avvio e adozione

Progetto nuovo. Il PI crea il repository con `tools/progetto.py` e con la copia di questa skill, e scrive l'idea in poche righe: domanda, dati, lettore, vincoli. La prima sessione esegue `init` ed è di scommessa.

Progetto già avviato. Si adotta come un progetto nuovo che ha già del codice: l'obiettivo non si giudica dallo stato del codice.

1. Scommessa: esegui `init` e scrivi il contratto partendo dall'obiettivo, non da ciò che esiste. Contesta il criterio e cerca l'esistente come per un progetto nuovo.
2. Sonda: costruisci metro, tetto e catena banale. Il codice ereditato entra come un braccio in più, misurato con lo stesso metro sugli stessi dati.
3. Al cancello G1 porta tre numeri (tetto, banale, ereditato) e le voci del bilancio dell'errore. I pezzi ereditati che non battono il banale si tolgono. La costruzione riparte dalla catena migliore fra le due.
4. Sposta backlog e documenti storici in `docs/archivio/`, senza riscriverli: in apertura non si leggono più.

# Fonti

- Morris, White, Crowther, Using simulation studies to evaluate statistical methods, Statistics in Medicine (2019): https://doi.org/10.1002/sim.8086
- MacCoun, Perlmutter, Blind analysis: hide results to seek the truth, Nature (2015): https://doi.org/10.1038/526187a
- Whitesides, Whitesides' Group: Writing a Paper, Advanced Materials (2004): https://doi.org/10.1002/adma.200400767
- Mensh, Kording, Ten simple rules for structuring papers, PLOS Computational Biology (2017): https://doi.org/10.1371/journal.pcbi.1005619
