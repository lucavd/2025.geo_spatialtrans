## In breve: il simulatore regge una sezione intera, e la simulazione a pezzi di oggi sbaglia alle giunture?

**A cosa serviva.** Finora il simulatore era stato misurato su aree di 1 mm²; un vetrino vero ne copre 45. Qui abbiamo misurato come crescono tempo e memoria con l'area, quali errori fa una simulazione a pezzi fatta con le funzioni di oggi (ogni pezzo simulato da solo, con un bordo di tessuto attorno), e quanto costa una sezione intera dei cinque vetrini veri. Serve a decidere che cosa sistemare prima di chiudere lo Step 2.

**Cosa abbiamo trovato.**

1. **In un pezzo solo il simulatore cresce in proporzione all'area e sta nei 5 minuti per mm², ma per le sezioni più cariche la memoria del server non basta.**
   - Da 1 a 9 mm² il tempo per mm² resta lo stesso: da 1 minuto per mm² nel cuore a 4 minuti e mezzo nel tumore del colon.
   - Abbiamo simulato in un pezzo solo la sezione intera dell'intestino (107 GB, quasi 2 ore), del linfonodo (87 GB) e del cuore (26 GB). Il tumore del colon e il cervello supererebbero i 251 GB del server.
   - Su una sezione intera si generano più trascritti del previsto: per fare le strisce il simulatore genera tanti trascritti quanti ne servono alla striscia più forte di tutto il vetrino, e su tutta la sezione del linfonodo sono 7,6 volte quelli che restano, contro 3-6 volte sui quadrati di 9 mm².
2. **A pezzi, il bordo di tessuto attorno a ogni pezzo costa troppo e lascia due difetti.**
   - Con i pezzi da mezzo millimetro del documento di progetto si simula 4,5 volte l'area utile: il tempo per mm² va da quasi 2 a più di 4 volte il bilancio in cinque tessuti su sei (non nel cuore). Con pezzi da 1 mm una sezione intera costa da 2 a 7 ore-processo.
   - Le cellule del bordo vengono simulate due volte con due profili diversi: nel blocco accanto alla giuntura fanno da un sesto a due quinti dei trascritti, e la somiglianza fra blocchi vicini attraverso la giuntura cala da un terzo a più della metà.
   - Le strisce cambiano da un pezzo all'altro, per costruzione.
3. **Con il bordo, i trascritti medi vicino alle giunture restano giusti dove la misura è abbastanza precisa; senza bordo calano dal 6 al 13 %.**
   - Nel linfonodo restano entro il 2 % in 8 fasce di distanza su 10; negli altri tessuti tre semi non bastano per dirlo. Un calo del 5 % nello stroma sembrava reale, ma era il caso: con nove semi in un pezzo solo sparisce.
   - Senza bordo il calo accanto alla giuntura si vede chiaramente in quattro tessuti su sei.
4. **Sul vetro di una sezione intera, la simulazione a pezzi sposta la coda lunga dei trascritti.**
   - Nell'intestino e nel linfonodo, dove la coda arriva a 2 mm, il vetro a più di mezzo millimetro dal tessuto ha la metà dei trascritti; quello vicino al tessuto fino alla metà in più. Nel cuore, dove la coda è corta, a pezzi e in un pezzo solo coincidono, salvo un angolo del vetrino.
   - Il motivo: ogni pezzo ripiega la coda dentro il suo riquadro, anche sul bordo del vetrino, invece di lasciarla arrivare dove deve.
5. **Il livello delle cellule non è il collo di bottiglia, ma su una sezione intera cresce più che in proporzione.**
   - Sui quadrati costa un secondo più 0,7 millesimi di secondo per cellula. Sulla sezione intera dell'intestino invece richiede 36 minuti, quasi tutti nel ritaglio delle celle sul bordo a gradini della maschera del tessuto. Nel linfonodo arriva a 10,5 GB, sopra il bilancio di 8.
   - Fra le 631 000 cellule del linfonodo, una ha uno spigolo di meno di un decimillesimo di micron su cui il calcolo dei triangoli di GEOS si blocca e ferma tutta la simulazione. Per questa prova l'ho corretta a mano.

**Cosa devi decidere.** Niente: la regola scritta prima del calcolo, con la tua delega, sceglie lo step successivo.

**Comunque**: il prossimo step scriverà la simulazione a pezzi dentro il pacchetto, con ogni cellula simulata una volta sola, i trascritti che escono da un pezzo passati al pezzo dove arrivano (senza bordo simulato due volte), le strisce e la coda lunga calcolate sull'intera area, pezzi multipli di 8 µm, e una sovra-generazione che non cresca con l'area. Lì renderemo robusto anche il calcolo dei triangoli, e misureremo le giunture con più semi. L'ottimizzazione del ritaglio nel livello delle cellule e il vetro vicino al tessuto del linfonodo, che nel simulato è circa tre volte quello vero, restano in coda per quando li riprenderemo.

Non c'è nulla da decidere: se non hai dubbi, lo step è chiuso.

**Sui tempi di calcolo.** Misurati da me sul server: 52 ore-processo in tutto (crescita con l'area 9, giunture 16, sezioni intere 25, verifiche successive 2), contro le 44 stimate all'apertura. Circa 7 ore e tre quarti di orologio, con fino a 25 processi e la temperatura sotto gli 85 °C. In parte la macchina era condivisa con altri due lavori, quindi i tempi delle giunture e delle sezioni sono a macchina carica; i tempi della crescita con l'area sono stati misurati con al più 6 processi.
