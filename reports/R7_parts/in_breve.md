### In breve: quanto sono troppo stretti gli intervalli dei valori veri presi da 5 aree per tessuto, e cambia qualche conclusione?

**A cosa serviva.** Quasi tutti i valori veri con cui misuriamo il simulatore (densità e aree dei nuclei, trascritti per quadratino, tipi cellulari) vengono da 5 aree di 1 mm² per tessuto, e la loro incertezza è calcolata fra queste 5 aree. Le aree sono poche e in parte si sovrappongono. Volevamo sapere di quanto quegli intervalli sono troppo stretti, correggerli e vedere se qualche conclusione passata cambia, prima di confrontare il simulatore completo con i tessuti veri.

**Cosa abbiamo trovato.**

1. **Gli intervalli che usiamo contengono il valore vero del tessuto fra 76 e 84 volte su 100, non 95.**
   - La causa principale è il numero di aree: anche con 5 aree indipendenti si arriva solo a 84 su 100.
   - Le zone comuni peggiorano le cose dove sono più grandi: 76 su 100 nello stroma del tumore, 77 nel linfonodo.
   - L'intervallo fra la più bassa e la più alta delle 5 aree, che usiamo come margine, regge meglio: fra 89 e 94 su 100.
2. **Per tornare a 95 su 100 gli intervalli vanno allargati di circa 1,6-1,9 volte nel cuore, nel tumore del colon e nell'intestino, e di 2,6-3,2 volte nello stroma e nel linfonodo.**
   - Nei primi tre basta la geometria delle zone comuni.
   - Nello stroma e nel linfonodo il tessuto ha strutture di 270-480 µm, per esempio i follicoli del linfonodo. Così due aree che si toccano si somigliano molto più di quanto dice la sola parte in comune, e il fattore viene dalla struttura misurata dentro le aree.
   - Nel cervello di topo il fattore è 1,8 per i trascritti per quadratino e 3,9 per le misure sui nuclei. Quest'ultimo è incerto: la densità dei nuclei cambia su distanze più lunghe delle aree, e la stima va oltre i dati.
   - La correzione regge con valori a code pesanti. Non basta del tutto quando i valori delle aree sono molto asimmetrici, per esempio frazioni vicine a zero o tipi rari (24 intervalli su circa 240 controllabili): per quelli abbiamo calcolato un allargamento in più.
3. **Contare due volte i nuclei delle zone comuni non sposta i valori veri.** Lo spostamento è al più del 7 %, in una frazione molto piccola (i quadratini vuoti nel linfonodo), e sempre meno di metà dell'intervallo corretto. Vanno corretti gli intervalli, non i valori.
4. **Con gli intervalli corretti un verdetto passato su dieci cambia, ma nessuna scelta del simulatore oggi in uso è smentita.**
   - Su circa 7 900 verdetti di 29 step ne cambiano 709. Circa metà diventano non decisi; poco meno di metà diventano confermati, perché si allargano anche i margini veri. Nessun verdetto passa da confermato a smentito; 22 passano da non decisi a smentiti.
   - Quattro regole del livello cellulare ora finiscono in pareggio o con un punto di scarto, e per come erano scritte sceglierebbero l'altra variante. Riguardano come si assegna la taglia ai nuclei, la distanza minima fra nuclei vicini, la forma della variabilità delle taglie e le quote dei tipi. Anche la verifica finale delle cellule passa da confermata a non decisa (163 conferme, ne servivano 167). Ognuno dei quattro cambi di regola dipende da una sola riga o da un pareggio, e torna indietro con piccole variazioni della correzione. Quindi, con l'incertezza onesta, i dati non distinguono quelle varianti.
   - Il limite dichiarato «troppi nuclei a contatto nel tumore del colon e nel linfonodo» non regge più nel linfonodo, in nessuna delle letture corrette. Nel tumore del colon diventa non deciso o compatibile, secondo la lettura. Nello stroma la risposta dipende da come si calcola l'incertezza dei conteggi fatti a mano.
5. **Nel confronto finale del simulatore con i tessuti veri, gli intervalli dei trascritti per quadratino diventano molto più larghi nello stroma e nel linfonodo.** Per esempio, nello stroma i trascritti per quadratino da 8 µm passano da 260-350 a 144-431.

**Cosa devi decidere.**

1. **Adottiamo gli intervalli corretti come valori veri da qui in avanti, compreso il confronto finale? Consiglio di sì.**
   - Con quelli di oggi il simulatore risulterebbe «fuori dal vero» troppo spesso, per un difetto della misura e non suo. I fattori vengono dai dati, con regole scritte prima del calcolo e ricontrollate da quattro revisori.
   - Per i conteggi fatti a mano nelle finestre consiglio di prendere, fra il vecchio intervallo e quello corretto, il più largo: nello stroma quello corretto verrebbe più stretto di quanto permettono 60 nuclei contati.
   - Nel cervello, per le misure sui nuclei, consiglio il fattore grande, come dice la regola, segnato come incerto. Nel confronto finale conta poco, perché lì contano i trascritti per quadratino, che nel cervello cambiano del 3 %.
   - Alternativa: usare solo la correzione geometrica. È più semplice, ma nello stroma, nel linfonodo e nel cervello gli intervalli conterrebbero il valore vero solo fra 74 e 90 volte su 100.
   - Alternativa: tenere gli intervalli di oggi, che contengono il valore vero fra 76 e 84 volte su 100.
2. **Teniamo le varianti oggi in uso nel livello cellulare, scrivendo che i dati non le distinguono dalle alternative? Consiglio di sì.**
   - Le regole cambiano solo per pareggi o per un punto, e tornano indietro con piccole variazioni. Nessuna variante risulta migliore.
   - Alternativa: applicare le regole come sono scritte. Cambierebbero quattro scelte del modello che genera le cellule, e servirebbe uno step per ritarare e risimulare il livello cellulare, per un guadagno che i dati non mostrano.
   - Alternativa: raccogliere i dati che servono per distinguerle davvero, cioè almeno 10 aree disgiunte per tessuto e altre finestre contate a mano. La segmentazione costa poco; costano rifare le misure dei passi successivi e le tue ore di conteggio, che non ho stimato.
3. **Riscriviamo il limite dichiarato sui contatti fra nuclei come lo vedono ora i dati? Consiglio di sì.**
   - Diventerebbe: nel linfonodo nessun eccesso dimostrato; nel tumore del colon non deciso; nello stroma e per le coppie di nuclei a contatto il limite resta.
   - Alternativa: lasciarlo com'è fino a nuove finestre contate a mano.

**Comunque, nei prossimi step:**

- quando confronteremo il simulatore completo con i tessuti veri useremo la tabella dei valori veri corretti, e ogni nuovo intervallo calcolato su poche aree o finestre avrà la sua correzione;
- quando rivedremo le tarature del livello cellulare terremo conto che le validazioni lasciando fuori un'area condividevano fino al 60 % dei nuclei con le aree di stima, nel linfonodo e nello stroma;
- quando torneremo sulle strisce e sul fondo del vetro misureremo lo stesso problema sugli intervalli calcolati su tessere vicine della stessa sezione.

Appena decidi questi 3 punti, lo step è chiuso.

**Sui tempi di calcolo.** Tutto il calcolo di questo step è girato sul server con al massimo 6 processi. La parte più lunga, la struttura del tessuto con il bootstrap rifatta dopo la revisione, ha richiesto circa 1 ora e 15 minuti con 6 processi (dall'avvio del processo sul server alla notifica di fine del job). Per l'alternativa delle aree nuove ho solo il tempo misurato di Cellpose, circa 1,2 secondi per un'immagine di 1920 × 1216 pixel; non ho misurato quanto costa rifare i passi successivi, né le ore di conteggio a mano.
