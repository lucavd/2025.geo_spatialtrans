## In breve: quanti trascritti cadono nei quadratini di Visium HD nei tessuti veri, e come si dividono fra nucleo, citoplasma e cellule?

**A cosa serviva.** Prima di scrivere la parte del simulatore che mette i trascritti sulla griglia di Visium HD servono valori misurati: quanti trascritti cadono in un quadratino da 2 µm e da 8 µm, quanti quadratini restano vuoti, quanti trascritti ha una cellula, se il nucleo è più o meno denso del citoplasma, se serve un passo che azzera valori a caso («dropout»). La vecchia pipeline usava valori senza fonte: 8 000 trascritti per cellula, un dropout che azzera circa metà dei valori, un tetto a 50 000. Li abbiamo misurati nelle 30 aree di riferimento (5 per tessuto), con i nuclei veri e il tipo cellulare che avevamo dato loro dall'espressione.

**Cosa abbiamo trovato.**

1. **In cinque tessuti su sei le cellule vere hanno da 7 a 100 volte meno trascritti di quanto assumeva la vecchia pipeline, e variano fra loro circa tre volte di più.**
   - Trascritti per cellula (mediana): 83 nel linfonodo, 273 nel cuore, 510 nell'intestino di topo, 800 nello stroma del tumore, 1 125 nel tumore del colon, 6 700 nel cervello di topo.
   - Nessuna cellula supera 50 000 (la più ricca, nel cervello, ne ha 28 000): il tetto della vecchia pipeline non toglieva nulla di vero.
   - Le aree di riferimento non somigliano alla media della sezione intera: nel linfonodo hanno 2.4 volte i trascritti per quadratino della media della sezione, nello stroma 0.6 volte. I valori dei preset vanno presi dalle aree del tessuto, non dai riassunti che 10x pubblica per ogni sezione.
2. **Non serve un passo di dropout: gli zeri si spiegano con le quantità diverse di trascritti fra cellule.**
   - Dentro un tipo cellulare, ogni gene ha tanti zeri quanti ne prevede il modello di conteggio che il simulatore usa già, purché si tenga conto che alcune cellule hanno più trascritti di altre: così nel tumore del colon, nello stroma, nel cervello e nel cuore.
   - Nell'intestino restano zeri in più, ma vengono da sottotipi dentro i tipi grossolani (cellule di Paneth, enterociti in punti diversi del villo, plasmacellule fra le cellule immunitarie), non da perdite tecniche.
   - Nel linfonodo le cellule hanno troppo pochi trascritti per giudicare.
   - Un dropout artificiale del 30 % aggiunto ai dati veri la misura lo vede sempre.
3. **Il nucleo non ha la stessa densità di trascritti del citoplasma in quattro tessuti su sei, e la direzione dipende dal tessuto.**
   - Nel cervello il nucleo ha circa il 40 % di trascritti in più per quadratino rispetto alla corona di 8 µm attorno, nel cuore quasi l'80 % in più; nel tumore del colon il 15 % in meno, nel linfonodo un quinto in meno.
   - Nell'intestino la densità è uniforme; nello stroma la differenza sta a cavallo di un decimo, e vale la densità uniforme come avevamo stabilito per i casi incerti.
   - La prima misura, fatta sui quadratini raggruppati, era disturbata dai compartimenti del tessuto (dove ci sono più nuclei ci sono anche più trascritti); misurata cellula per cellula il disturbo sparisce e le direzioni restano le stesse.
4. **I trascritti per cellula cambiano con il tipo cellulare.** Nel cervello i neuroni ne hanno quasi 4 volte gli oligodendrociti; differenze oltre un quarto anche nell'intestino, nello stroma e nel cuore. Nel tumore del colon e nel linfonodo la misura sta a cavallo della soglia di un quarto (primo punto da decidere).
5. **A 2 µm la griglia di cattura alterna righe e colonne «forti» e «deboli».**
   - Nell'intestino e nel tumore del colon le colonne dispari raccolgono 1.7 volte i trascritti delle pari, nel cervello 1.2, nel linfonodo 1.1; nel cuore l'effetto è minimo e rovesciato.
   - È l'effetto «a strisce» già descritto da chi ha sviluppato Bin2cell. Il reticolo che avevamo visto in un lume del tumore del colon è questo: un altro reticolo, diagonale, non c'è.
   - A 8 µm l'effetto sparisce, perché ogni quadratino da 8 µm contiene metà righe pari e metà dispari.

**Cosa devi decidere.**

1. **Diamo a ogni tipo cellulare il suo numero di trascritti in tutti e sei i tessuti? Consiglio di sì.**
   - In quattro tessuti la differenza fra tipi supera chiaramente un quarto. Nel tumore del colon e nel linfonodo la misura sta a cavallo, e la regola che avevamo scritto prima del calcolo dà allora un solo valore per tessuto.
   - Ma nel tumore del colon le cellule epiteliali tumorali hanno meno di un terzo dei trascritti dell'epitelio normale (373 contro 1 298): restano fuori dal confronto solo perché sono il 4.7 % dei nuclei invece del 5 %.
   - Usare il valore per tipo ovunque non costa nulla: dove i tipi si somigliano, i valori sono vicini.
   - Alternativa: seguire la regola come scritta, un solo valore nel tumore del colon e nel linfonodo. Conseguenza: nel tumore del colon simulato le cellule tumorali avranno circa tre volte i loro trascritti veri.
2. **Il simulatore deve riprodurre a 2 µm l'alternanza di righe e colonne? Consiglio di sì, come strato tecnico che si può spegnere.**
   - I metodi da testare lavorano sempre più a 2 µm, e alcuni (come Bin2cell) correggono proprio queste strisce: senza l'effetto i dati simulati sarebbero più puliti dei veri e quella correzione non si potrebbe mettere alla prova.
   - A 8 µm non cambia nulla.
   - Alternativa: no. Conseguenza: i dati simulati a 2 µm saranno ideali, utili solo per metodi che lavorano a 8 µm o più.

**Comunque**, quando scriveremo il design della griglia: i trascritti per quadratino, per cellula e per tipo misurati qui diventano i riferimenti dei sei tessuti; spariscono il tetto, il dropout e le soglie a memoria della vecchia pipeline; la forza delle strisce si misura sezione per sezione, anche nella parte che non alterna in modo regolare; la densità del nucleo va tarata insieme alla diffusione dei trascritti, che misureremo nello step sulla diffusione, perché togliendo i quadratini al bordo del nucleo la differenza cresce (nel cuore il rapporto passa da 2 a 4); nell'intestino i profili di espressione dovranno avere sottotipi. Gli zeri nel linfonodo restano da giudicare quando avremo cellule con più trascritti.

Appena decidi questi 2 punti, lo step è chiuso.
