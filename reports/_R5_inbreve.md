## In breve: di quanto si spostano i trascritti fra il tessuto e i quadratini di Visium HD, e quanto è diverso davvero il nucleo dal citoplasma?

**A cosa serviva.** Il simulatore mette i trascritti nelle cellule e poi li deve spostare come fa la tecnologia vera, prima di contarli nei quadratini. Il design prevedeva una sfocatura «a campana» di larghezza da misurare (il commento di Daniele). L'abbiamo misurata dove si vede meglio: ai bordi del tessuto, dove si passa al vetro nudo, in 314 finestre di bordo delle cinque sezioni vere. Poi abbiamo controllato se la stessa sfocatura spiega quello che si vede dentro il tessuto (la differenza fra nucleo e citoplasma) e nei buchi vuoti del tessuto.

**Cosa abbiamo trovato.**

1. **Ai bordi i trascritti si allargano di 3-5 µm, e una parte va molto più lontano.**
   - Larghezza: 3.1 µm nel cervello, 3.4 nel linfonodo, 3.5 nel tumore del colon, 5.2 nell'intestino (oltre i 4 µm che avevamo fissato come limite). Nel cuore 3.2, ma con 5 sole finestre utilizzabili: incerto.
   - In tutte le sezioni c'è una coda: dal 6 % (tumore) al 15 % (cervello) dei trascritti finisce a decine di µm. Una campana sola non basta.
   - Nessuna direzione privilegiata fra righe e colonne della griglia.
   - La larghezza dipende da come si traccia il bordo nell'immagine: con un bordo tracciato in modo diverso cresce di 0.6-2 µm.
2. **Dentro il tessuto la sfocatura è molto più stretta, sotto 1 µm quasi ovunque.**
   - La differenza osservata fra nucleo e citoplasma si spiega con una sfocatura di 0.6-0.8 µm nel tumore, 0.7-2.5 µm nel cervello, circa mezzo µm o meno nello stroma e nel linfonodo. La procedura ritrova i valori giusti su dati costruiti apposta.
   - Con la larghezza dei bordi, la differenza osservata nel linfonodo sarebbe impossibile: nessun rapporto vero la produce.
   - Viceversa, imporre una sfocatura stretta ai bordi peggiora molto l'adattamento: i bordi sono davvero sfumati su 3-5 µm. Sono due scale diverse, e non sappiamo ancora perché (bordo del tessuto sottile o sfrangiato, bordo dell'immagine impreciso, trasporto sul vetro libero).
3. **Prima della diffusione il nucleo è un po' più diverso dal citoplasma di quanto si osserva.**
   - Nel tumore il nucleo ha 0.7-0.8 volte la densità del citoplasma (osservato circa 0.85); nel cervello circa 1.5 volte (osservato 1.4).
   - Nell'intestino e nel cuore il modello «nucleo e citoplasma con due densità» non riproduce i dati; nello stroma e nel linfonodo la sfocatura è troppo piccola per misurare la correzione.
   - Con la sfocatura dei bordi i valori sarebbero assurdi (nel tumore un nucleo con un sesto della densità del citoplasma, nel cuore 27 volte): solo nel cervello le due strade concordano.
4. **I buchi vuoti vicino al bordo hanno i trascritti previsti dalla sfocatura dei bordi** (entro un quarto in 9 fasce su 10, fino a 8 µm dal bordo del buco). Più all'interno, nel tumore del colon c'è più RNA del previsto (materiale nelle ghiandole), nel cuore meno.
5. **Nell'intestino l'immagine è spostata di 2.4 µm rispetto alla griglia.** Lo stesso spostamento si ritrova dentro il tessuto: allineando le maschere dei nuclei la differenza nucleo/citoplasma cresce, spostandole al contrario si inverte. Nelle altre sezioni lo spostamento è sotto 1 µm o incerto. Conseguenza: nelle impronte dei nuclei dell'intestino che abbiamo usato per i tipi cellulari solo un terzo dei quadratini sta al posto giusto.
6. **Il revisore numerico ha trovato un difetto grave nell'adattamento, e tre controlli scritti prima non sono passati come scritti.**
   - Il difetto: molte finestre restavano bloccate a un limite che non avevamo scritto. È corretto, e tutto è stato rifatto due volte con risultati identici. I numeri sono cambiati (intestino da 4.9 a 5.2 µm; coda del cuore da 0.50 a 0.38), le conclusioni no.
   - I tre controlli: i test su dati costruiti (98 prove su 105), dove falliscono i casi in cui il dato non contiene l'informazione; la lettura della matrice, che contiene geni in più; la controprova dello spostamento della maschera, dove era sbagliata la mia attesa.

**Cosa devi decidere.**

1. **Come sposta i trascritti il simulatore? Consiglio due livelli: una sfocatura piccola per tutti (circa 0.5-1 µm, per tessuto) e, per i trascritti che finiscono fuori dal tessuto (vetro, buchi), l'alone misurato ai bordi (3-5 µm più la coda).**
   - Così i dati simulati hanno sia la differenza nucleo/citoplasma dei dati veri sia l'alone sul vetro e nei buchi, che serve per mettere alla prova i metodi che cercano lo spazio vuoto.
   - Alternativa: una sola campana con la larghezza dei bordi, come nel design. Conseguenza: la differenza fra nucleo e citoplasma sparisce, e i metodi che assegnano i trascritti alle cellule avrebbero dati più difficili dei veri.
   - Alternativa: una sola campana stretta. Conseguenza: niente alone fuori dal tessuto, bordi più netti dei veri.
   - Il come (quali trascritti «escono» e con che regola) si scrive nel design dello Step 2.
2. **Chiudiamo lo step anche se tre controlli non sono passati come scritti? Consiglio di sì.**
   - In tutti e tre il motivo sta nel controllo, non nei dati; il metodo usato è quello registrato, e il revisore ha confermato a identità la lettura dei dati, le tabelle e il modello.
   - Alternativa: rifare i test con geometrie in cui lo spostamento si possa misurare. Conseguenza: un altro giorno di calcolo; i numeri sui dati veri non cambiano.
3. **Che cosa usiamo per il cuore? Consiglio i valori mediani delle altre quattro sezioni, con la riserva scritta nel preset.**
   - Nel cuore ci sono solo 5 bordi utilizzabili, la coda non si separa dalla campana, lo spostamento misurato ai bordi non vale dentro, e il modello dei nuclei non regge.
   - Alternativa: usare le stime del cuore così come sono. Conseguenza: preset con incertezza grande e una coda di peso incerto (0.15-0.62).

**Comunque**: la sfocatura entra nei preset per tessuto, uguale in tutte le direzioni e con la coda; lo spostamento di 2.4 µm dell'intestino va in BACKLOG come riserva sulle impronte dei nuclei di R6; i buchi vuoti si possono riprendere con la sfocatura dei bordi come riferimento per i primi 8 µm; resta da capire perché i bordi siano sfumati più dell'interno.

Appena decidi questi 3 punti, lo step è chiuso.
