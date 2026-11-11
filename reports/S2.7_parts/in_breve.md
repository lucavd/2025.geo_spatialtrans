## In breve: la simulazione a pezzi dentro il simulatore dà lo stesso vetrino della simulazione in un pezzo solo, nel tempo e nella memoria del bilancio?

**A cosa serviva.** Un vetrino vero copre circa 45 mm². Il simulatore di oggi lo regge in un pezzo solo solo per tre tessuti su cinque, e la simulazione a pezzi provata la volta scorsa sbagliava alle giunture: cellule simulate due volte, coda lunga dei trascritti ripiegata dentro ogni pezzo. Qui abbiamo scritto dentro il simulatore una simulazione a pezzi in cui ogni cellula si simula una volta sola e i trascritti che escono da un pezzo passano al pezzo dove arrivano. Poi abbiamo controllato che dia lo stesso risultato del pezzo solo e del simulatore di oggi, e quanto costa.

**Cosa abbiamo trovato.**

1. **Il risultato non dipende più da come si taglia il vetrino: a pezzi e in un pezzo solo i file sono identici, bit per bit.**
   - Lo abbiamo verificato su quadrati di 4 mm² dei sei tessuti, sulla sezione intera del cuore e sulle sezioni del tumore del colon e del cervello, con pezzi da mezzo millimetro e da un millimetro e con un numero diverso di processi.
   - Le giunture non esistono più come difetto: fra i pezzi di una sezione passano da 5 a 52 milioni di trascritti e nessuno si perde.
2. **Il simulatore a pezzi riproduce quello di oggi: UMI e geni per quadratino coincidono entro il 2 % in tutti i sei tessuti.**
   - Nel cervello, con nove ripetizioni, il nuovo dava l'1,6 % in più. Con 27 ripetizioni lo scarto scende allo 0,4 %: era il caso.
   - La prova vede una differenza quando c'è: abbassando apposta la cattura del 5 %, la vede in tutti e sei i tessuti.
3. **Sulle sezioni intere il vetro lontano dal tessuto ora è giusto.**
   - Nell'intestino, nel linfonodo e nel cuore la sezione a pezzi coincide con quella in un pezzo solo: nel tessuto entro l'1 %. Sul vetro a più di mezzo millimetro dal tessuto, nell'intestino e nel linfonodo, entro il 2 %. Prima, a pezzi, quel vetro aveva la metà dei trascritti.
   - Ora si simulano anche il tumore del colon e il cervello, che in un pezzo solo non stavano nella memoria del server.
4. **La memoria sta nel bilancio con pezzi da mezzo millimetro; il tempo no, nel tumore e nell'intestino.**
   - Con pezzi da un millimetro un pezzo del tumore usa 12,5 GB e uno del cervello 9,6, sopra gli 8 GB del bilancio. Con pezzi da mezzo millimetro scendono sotto i 4 GB e il risultato è lo stesso. Ho reso il mezzo millimetro il valore predefinito.
   - Per mm² di tessuto il tumore costa 6 minuti e mezzo e l'intestino 5 e mezzo, contro i 5 del bilancio; linfonodo, cervello e cuore stanno sotto. Una sezione intera del tumore richiede 3 ore e mezza di processore, mezz'ora con 8 processi. Da un sesto a due terzi del tempo va nel calcolo esatto di quanta parte di ogni quadratino copre ogni cellula.
5. **Il tetto ai trascritti generati per le strisce funziona, e abbiamo trovato e corretto un errore del simulatore di oggi.**
   - Con il tetto fisso che hai scelto, sulla sezione del linfonodo si generano 265 milioni di trascritti invece di 496. Il taglio tocca circa un quadratino su 100 000 e meno di 6 UMI su un milione.
   - Sulle aree molto grandi il simulatore di oggi poteva attribuire alcuni trascritti al gene sbagliato. Nelle sezioni intere simulate la volta scorsa in un pezzo solo era successo a circa l'1,3 % dei trascritti nell'intestino e nel linfonodo. Il numero totale di trascritti non cambiava, e le conclusioni di allora usavano solo quello. Ora è corretto, e la cellula del linfonodo che bloccava il calcolo dei triangoli si gestisce da sola.

**Cosa devi decidere.**

1. **Accettiamo che tumore del colon e intestino costino più di 5 minuti per mm² di tessuto? Consiglio di sì, come limite dichiarato.**
   - Una sezione intera resta fattibile in mezz'ora con 8 processi, e la memoria è nel bilancio.
   - Il guadagno possibile sta nel calcolo esatto delle coperture e nella posizione dei trascritti. È lavoro di ottimizzazione da fare insieme a quella del livello delle cellule, già in coda.
   - Alternativa: un passo di ottimizzazione subito, prima di chiudere lo Step 2. Costa un altro step e rimanda la revisione d'insieme.

**Comunque**: quando chiuderemo lo Step 2 sceglieremo quale delle due funzioni, quella a pezzi o quella di oggi, diventa l'ingresso del pacchetto. In quell'occasione scriveremo anche nella guida che il riquadro dell'area vale come bordo aperto sia per la coda lunga sia per la distanza dal bordo del tessuto. La prova del riquadro della finestra nell'intestino, rimasta indietro, si fa lì.

Appena decidi questo punto, lo step è chiuso.

**Sui tempi di calcolo.** Li ho misurati io sul server: 26 ore-processo per le prove pre-registrate, contro le 21 stimate, più 9 per le verifiche successive. In tutto circa 4 ore di orologio, con al massimo 16 processi e la temperatura al massimo di 82 °C. Una prima esecuzione di circa un'ora e mezza è stata cancellata e rifatta, dopo aver tolto una lentezza del codice che non cambiava i risultati.
