
### In breve: i metodi di analisi si comportano sui dati simulati come sui dati veri, e che cosa va sistemato per primo nel simulatore?

**A cosa serviva.** Il simulatore deve servire a mettere alla prova i metodi di analisi: un metodo dovrebbe comportarsi sui dati simulati come su quelli veri. Abbiamo fatto girare i due metodi più usati su Visium HD, bin2cell (ricostruisce le cellule allargando i nuclei) e RCTD (dice quali tipi cellulari ci sono in ogni quadrato di 8 µm), sulle stesse 30 aree dei sei tessuti, simulate e vere. Il margine accettabile era quanto ogni misura varia fra le aree vere dello stesso tessuto.

**Cosa abbiamo trovato.**

1. **Su metà dei confronti, 18 su 36, i metodi si comportano sul simulato in modo diverso che sul vero, soprattutto nell'intestino di topo, nel tumore del colon, nello stroma del tumore e nel linfonodo.**
   - Le cellule ricostruite da bin2cell sono più «sporche»: nel tumore il 15 % dei trascritti dei geni marcatori viene da altri tipi cellulari, contro meno del 5 % nel vero.
   - RCTD trova più miscugli di due tipi: nel tumore quasi un quadrato su cinque, contro uno su venti nel vero.
   - I tipi trovati da RCTD sono meno ordinati nello spazio: nello stroma due quadrati vicini hanno lo stesso tipo nel 65 % dei casi, contro l'89 % nel vero.
   - Il numero di trascritti per cellula invece torna in tutti e sei i tessuti.
2. **La causa principale è che il simulatore mette i tipi cellulari a caso uno accanto all'altro, mentre nei tessuti veri si raggruppano.**
   - Nell'intestino simulato quasi la metà dei vicini di una cellula è di un altro tipo, nel vero un quinto; nel tumore un terzo contro uno su sedici.
   - Questo spiega quasi tutti i trascritti estranei in più dell'intestino e due terzi di quelli del tumore.
   - Effetto visibile: nell'intestino simulato RCTD vede lo stroma in un quadrato su dodici, nel vero in uno su quattro, anche se il simulatore ha la stessa quota di cellule stromali dei nuclei veri (una su cinque). Sparse fra le cellule epiteliali, nei quadrati misti spariscono.
3. **Il simulatore ha già il rimedio: la semina a compartimenti, che oggi non è quella predefinita. Con questa i vicini si mescolano come nel vero, e i metodi si avvicinano al vero in quattro tessuti.**
   - Si chiude da un terzo a quasi tutta la differenza su miscugli, ordine nello spazio e accordo con il nucleo, e nell'intestino e nel tumore circa metà di quella sulle cellule «sporche»; 5 delle 17 differenze spariscono.
   - Non basta per tutto: l'intestino resta un po' più «sporco» del vero (anche per la densità troppo alta delle cellule del preset), e nel cervello di topo non cambia quasi niente.
4. **La mia ipotesi di partenza era sbagliata: i profili di espressione presi dai nuclei veri con la loro corona non sono la causa.** Con profili presi dal solo nucleo le cellule simulate restano sporche come prima, tranne nel cervello.
5. **Restano tre cose non spiegate e un difetto del codice.**
   - Nel linfonodo RCTD chiama plasmacellule il 30 % dei quadrati simulati contro il 13 % dei veri, e la disposizione dei tipi non lo spiega.
   - Nel cuore, dove ci aspettavamo che il simulatore sbagliasse perché le cellule vere sono allungate, nessuna di queste misure vede la differenza.
   - Le cinque aree simulate di un tessuto si somigliano fra loro molto più di quanto si somiglino le cinque aree vere.
   - Con la semina a compartimenti una delle cinque aree del cervello si ferma con un errore di geometria.

**Cosa devi decidere.**

1. **Facciamo diventare predefinita la semina a compartimenti, in uno step dedicato? Consiglio di sì.**
   - È la scelta che più avvicina il comportamento dei metodi al vero, e il codice c'è già, controllato quando abbiamo costruito le cellule.
   - Lo step dedicato correggerebbe l'errore di geometria del cervello, rifarebbe i controlli delle cellule e dei quadratini con la nuova semina, e rifarebbe queste misure dei metodi con due semi come verifica.
   - Alternativa: lasciare la semina uniforme come predefinita e dire nella guida d'uso di attivare i compartimenti per provare i metodi. Chi usa il simulatore senza leggerlo troverebbe metodi più in difficoltà che sui dati veri.
   - Alternativa: farlo insieme alla revisione della densità delle cellule nei preset. Uno step solo, ma più lungo, e sarebbe più difficile capire quale cambiamento ha fatto cosa.

**Comunque**

- Quando cambieremo la semina o la densità, rifaremo queste misure come controllo: gli strumenti sono pronti e girano in poche ore.
- Le plasmacellule in eccesso del linfonodo, le aree simulate troppo simili fra loro e una misura che veda le cellule allungate del cuore restano in coda, ciascuna con quello che serve per affrontarla.
- Una lezione di metodo per i prossimi step: RCTD dà risultati leggermente diversi a ogni esecuzione se non si fissa il seme; ora lo fissiamo.

Appena decidi questo punto, lo step è chiuso.

**Sui tempi di calcolo.** Li ho misurati io sul server, con 6 processi accanto al lavoro di tesi. Le 30 aree vere con i due metodi hanno richiesto circa un'ora e mezza; le 85 simulazioni del pannello e della controprova con i due metodi circa 4 ore, dopo che una prima corsa è stata fermata per il seme di RCTD; le 30 aree con la semina a compartimenti circa 3 ore con 5 processi. Una simulazione di un'area richiede da 2 a 7 minuti, RCTD da 4 a 20 minuti per area (mediana 11). Rifare le misure con due semi nello step proposto costerebbe circa 4 ore con 6 processi.

