### In breve: il simulatore cattura i trascritti come un vetrino vero, con le strisce di righe e colonne, il fondo sul vetro e la coda lunga dello spostamento?

**A cosa serviva.** Questo step aggiunge al simulatore il passo che trasforma i trascritti spostati nei conteggi di un vetrino: righe e colonne della griglia catturano in modo diverso (le «strisce»), sul vetro c'è un fondo, e una piccola parte dei trascritti finisce a centinaia di micron di distanza. Abbiamo controllato che le strisce simulate siano indistinguibili da quelle vere con lo stesso stimatore, che il vetro intorno alle cinque sezioni intere sia descritto bene, e quanto costa in tempo e memoria.

**Cosa abbiamo trovato.**

1. **Le strisce simulate sono come quelle vere in tutti e sei i tessuti.**
   - Differenza fra righe pari e dispari, ampiezza dei fattori di riga e colonna, somiglianza fra righe vicine: nessuno scarto rilevante; dove l'esito resta incerto è per la precisione, non per uno scarto.
   - Lo stimatore non trova strisce quando le spegniamo: il confronto ha potere.
2. **Fondo e coda lunga, rifatti insieme come avevi scelto, descrivono il vetro di quattro sezioni su cinque a ogni distanza; nel tumore del colon solo fino a un quarto di millimetro.**
   - In intestino, linfonodo e cuore il fondo costante sparisce: tutto il vetro si spiega con una coda di trascritti che arriva a 1-2 mm (nel cuore a un decimo di millimetro).
   - Nel tumore del colon il vetro vero ha 7 decimi dei trascritti previsti fra un quarto e mezzo millimetro dal tessuto, e meno della metà oltre; lì fondo e coda si separano male.
   - La coda raddoppia i trascritti nei lumi profondi del tumore del colon, che passano da circa metà del vero a 7-9 decimi.
3. **La parte lenta delle strisce non entra nel simulatore.**
   - Sul vetro la troviamo solo in due sezioni su cinque, e il revisore ha mostrato che può essere un errore del modello del vetro più che una vera variazione della cattura. Resta un limite dichiarato.
4. **I quadratini vuoti veri non vengono dalle strisce.**
   - Nell'intestino, nel tumore e nel suo stroma il simulato ha da un decimo a quattro decimi dei quadratini vuoti veri.
5. **La memoria supera gli 8 GB del progetto in tumore, stroma e cervello, ma il rimedio c'è.**
   - Per avere conteggi interi e una verità esatta si generano da 3 a 6 volte i trascritti; il picco è nella funzione che li sposta (fino a 15-20 GB nel cervello).
   - Lavorando a blocchi più piccoli il cervello sta in 6 GB e meno di 3 minuti. Lo imposterà la funzione che metterà insieme tutti i passi.
   - Trovato e corretto anche un difetto dei generatori casuali: semi consecutivi davano estrazioni legate fra loro. Le funzioni nuove ne sono esenti; quelle dei passi precedenti vanno corrette.

**Cosa devi decidere.**

1. **Chiudiamo lo step anche se un controllo automatico è fallito come scritto? Consiglio di sì.**
   - Su 120 simulazioni, una ha catturato un numero di trascritti a 3.01 errori casuali dall'atteso, contro una soglia scritta di 3; sulle 120 è quanto ci si aspetta dal caso.
   - Alternativa: rifare il controllo con la soglia per confronti multipli che avevamo scritto in testa ai controlli; passerebbe, ma sarebbe cambiare la regola dopo il risultato.
2. **Accettiamo come limite che, nel tumore del colon, il vetro oltre un quarto di millimetro dal tessuto abbia nel simulato da una volta e mezza a più del doppio dei trascritti veri? Consiglio di sì.**
   - Riguarda il vetro lontano di un solo tessuto; il tessuto e il vetro vicino sono descritti bene.
   - Alternativa: uno step dedicato alla forma della coda (per esempio una coda con due lunghezze o un fondo che decresce con la distanza). Conseguenza: un pezzo in più nel modello, per un effetto piccolo e misurato su una sola sezione.

**Comunque**: la correzione dei generatori casuali dei profili e dello spostamento va in uno step dedicato; l'origine dei quadratini vuoti veri si studia con la revisione della densità dei preset; la funzione che metterà insieme tutti i passi sceglierà i blocchi per stare negli 8 GB; prima di quella funzione resta da fare lo step sulle aree di riferimento sovrapposte, come avevi deciso.

Appena decidi questi 2 punti, lo step è chiuso.

**Sui tempi di calcolo.** Misurati da me sul server condiviso (carico medio fino a 740 per il lavoro di un'altra utente), un processo per area: da 24 s (cuore) a 280 s (tumore del colon) per un'area di 1 mm² con tutta la catena; con i blocchi più piccoli il cervello impiega 145 s. Il pannello (120 simulazioni) ha girato con 6 processi in circa un'ora e mezza.
