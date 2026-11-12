## In breve: la guida d'uso e i documenti descrivono la griglia di Visium HD com'è, e lo Step 2 è entrato nel ramo di sviluppo senza perdere nulla?

**A cosa serviva.** Nello step precedente avevamo rivisto e fissato il codice che porta le cellule simulate sulla griglia di Visium HD. Mancavano i documenti. Serviva una guida d'uso per chi userà il simulatore per mettere alla prova i metodi. Il README e le istruzioni del progetto andavano aggiornati, e a Daniele serviva una guida per rivedere il codice. Poi il lavoro andava unito al ramo di sviluppo. Il codice del simulatore non doveva cambiare, e non è cambiato.

**Cosa abbiamo trovato.**

1. **La guida d'uso, in inglese, si esegue davvero: ogni volta che la si produce, i suoi esempi girano sul server, e due esecuzioni danno conteggi identici.**
   - Un disco di intestino di topo simulato porta il lettore attraverso tutto il simulatore: la verità con cui si valutano i metodi, la simulazione a pezzi, gli strati tecnici da accendere e spegnere, le aree ritagliate da un tessuto più grande, la ripresa di un calcolo interrotto, i tempi, la memoria e i limiti.
   - Una copia della guida con un esempio sbagliato apposta (pezzi di un lato non consentito) non si produce. Quindi gli esempi della guida girano davvero, non sono solo testo.
2. **README, istruzioni del progetto e guida per Daniele nominano solo file e funzioni che esistono, e la guida d'uso spiega ogni opzione della funzione d'ingresso.**
   - Lo stesso controllo boccia il README di prima: delle 15 opzioni e delle 20 funzioni pubbliche dello Step 2 ne nominava una sola.
   - La guida per Daniele mette accanto a ogni parte del documento di progetto il codice e i test che la coprono. Dice anche che cosa è diventata ciascuna delle sue osservazioni del 2025. Come per le cellule, non contiene domande.
3. **La revisione indipendente ha trovato una quarantina di imprecisioni nei documenti: le ho corrette tutte, salvo una frase di dettaglio messa in coda; ciò che va cambiato nel codice è in coda per il prossimo lavoro sul codice.**
   - Due avevano sostanza. Il tempo previsto non basta neanche per il cervello di topo intero con i pezzi predefiniti (è la seconda decisione qui sotto). E una simulazione interrotta prima di un aggiornamento di R non si può riprendere dopo: è voluto, ma la guida diceva il contrario.
   - Le altre erano descrizioni troppo semplici del modello o delle opzioni. Per esempio mancava la coda dello spostamento che porta alcuni trascritti a decine di µm dalla cellula.
4. **Il 5 novembre un aggiornamento di R sul server ha fatto fallire due controlli che confrontano i file di oggi con quelli salvati prima, anche se i dati sono identici.**
   - Le cellule dell'esempio del README sono le stesse di fine ottobre. Lo sono anche quelle di una prova del livello cellulare. Nei file cambia un solo byte: quello che registra la versione di R.
   - Finché non correggiamo quella prova, la batteria dei test delle cellule segnala un errore che non c'è.
5. **Tutte le prove passano in una copia nuova del repository presa da GitHub, salvo quei due confronti, e il ramo di sviluppo ora contiene lo Step 2.**
   - Le prove della griglia passano tutte. Quelle delle cellule passano con i soliti errori noti, più il confronto del punto 4.
   - Il ramo principale, fermo al 2025, non è stato toccato, come avevi scelto.

**Cosa devi decidere.**

1. **Chiudiamo lo step anche se i due confronti con i file vecchi sono falliti? Consiglio di sì.**
   - I dati sono identici, l'ho controllato oggetto per oggetto. A cambiare è solo la versione di R scritta nell'intestazione dei file.
   - Alternativa: correggere subito le due prove perché confrontino il contenuto e non il file. Tocca il codice dei test, quindi serve uno step a parte prima di chiudere questo.
2. **Il cervello di topo intero con i pezzi predefiniti supera il tempo previsto: lo aggiungiamo al limite che hai già accettato per il tumore del colon e l'intestino? Consiglio di sì.**
   - Impiega 335 secondi per mm² di tessuto contro i 300 del bilancio, circa un decimo in più. L'11 novembre il limite è stato scritto su prove con pezzi più grandi, e lì il cervello restava sotto.
   - Alternativa: fare l'ottimizzazione del livello cellulare e dei pezzi prima di passare ai metodi, invece che dopo. Guadagno e tempo di lavoro non sono ancora stimati.

**Comunque**

- La guida per Daniele è pronta: puoi mandargliela.
- Quando torneremo sul codice, correggeremo per prime le prove che confrontano i file salvati, perché ogni step che rifà la batteria delle cellule le incontrerà. Insieme: una cartella di ripresa senza il suo file di controllo verrà rifiutata; una prova della griglia segnalerà gli errori anche con lo stato d'uscita; la documentazione dentro il codice della funzione d'ingresso dirà ciò che dice già la guida.
- Il ramo principale lo aggiorneremo quando il simulatore avrà una prima versione usabile per i metodi.

Appena decidi questi 2 punti, lo step è chiuso.

**Sui tempi di calcolo.** Li ho misurati io sul server. Una produzione della guida richiede circa 3 minuti. Le prove complete in una copia nuova richiedono da 40 minuti a un'ora di orologio, con una decina di processi alla volta: l'ora è quella dell'ultimo giro, rallentato da altri lavori sul server. Le ho rifatte tre volte, perché le correzioni della revisione cambiavano i documenti provati. La stima era di circa un'ora di orologio per un solo giro.
