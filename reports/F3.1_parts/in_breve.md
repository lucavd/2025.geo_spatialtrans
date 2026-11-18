
### In breve: i metodi di analisi si comportano sui dati simulati come sui dati veri, e che cosa va sistemato per primo nel simulatore?

**A cosa serviva.** Il simulatore deve servire a mettere alla prova i metodi di analisi: un metodo dovrebbe comportarsi sui dati simulati come su quelli veri. Abbiamo fatto girare i due metodi più usati su Visium HD, bin2cell (ricostruisce le cellule allargando i nuclei) e RCTD (dice quali tipi cellulari ci sono in ogni quadrato di 8 µm), sulle stesse 30 aree dei sei tessuti, simulate e vere. Il margine accettabile era quanto ogni misura varia fra le aree vere dello stesso tessuto. Dopo la tua prima lettura abbiamo indagato, qui dentro, i quattro punti rimasti aperti.

**Cosa abbiamo trovato.**

1. **Su metà dei confronti, 18 su 36, i metodi si comportano sul simulato in modo diverso che sul vero, soprattutto nell'intestino di topo, nel tumore del colon, nello stroma del tumore e nel linfonodo.**
   - Le cellule ricostruite da bin2cell sono più «sporche»: nel tumore il 15 % dei trascritti dei geni marcatori viene da altri tipi cellulari, contro meno del 5 % nel vero.
   - RCTD trova più miscugli di due tipi (nel tumore quasi un quadrato su cinque, contro uno su venti) e tipi meno ordinati nello spazio (nello stroma due quadrati vicini hanno lo stesso tipo nel 65 % dei casi, contro l'89 %).
   - Il numero di trascritti per cellula invece torna in tutti e sei i tessuti.
2. **La causa principale resta la disposizione: il simulatore mette i tipi cellulari a caso uno accanto all'altro, mentre nei tessuti veri si raggruppano.**
   - Nell'intestino simulato quasi la metà dei vicini di una cellula è di un altro tipo, nel vero un quinto; nel tumore un terzo contro uno su sedici.
   - La semina a compartimenti, che il simulatore ha già, mescola i vicini come nel vero e avvicina i metodi al vero in quattro tessuti: 5 delle 17 differenze spariscono.
3. **Le indagini hanno trovato una seconda causa: il simulatore fa finire nei nuclei troppi trascritti delle cellule vicine.**
   - Nel linfonodo i geni degli anticorpi, tipici di linfociti B e plasmacellule, nei nuclei dei linfociti T simulati sono quasi cinque volte più che nei veri. Nel simulato un quinto di quello che si raccoglie sopra un nucleo T viene da linfociti B o plasmacellule; nel vero al massimo un ventesimo.
   - È questo, e non i profili dei linfociti B (che sono giusti), a far chiamare a RCTD plasmacellule il 30 % dei quadrati simulati contro il 13 % dei veri: tre quarti dell'eccesso viene da quadrati di altre cellule.
   - Lo stesso succede nell'intestino, nel tumore e nello stroma: i nuclei veri sono puliti quanto i profili di partenza, quelli simulati no. Anche con i profili presi dal solo nucleo l'eccesso resta, quindi nasce dove il simulatore mette e sposta i trascritti.
   - La semina a compartimenti ne toglie circa metà nell'intestino e nel tumore, meno nello stroma e nel linfonodo. Nell'intestino e nel linfonodo i vicini sono già mescolati come nel vero e l'eccesso resta: è una causa distinta dalla disposizione.
   - Il sospetto, da verificare, cade su quanto lontano il simulatore sposta i trascritti. Quei valori furono tarati su quantità che non distinguono il tipo di cellula; i geni marcatori, usati qui, sono la prima misura che lo vede.
4. **Il cuore: una misura giusta esiste già.** Le misure dei metodi non vedevano che le cellule vere del cuore sono allungate e quelle simulate no. Una misura dell'orientamento dei trascritti, la stessa usata nei controlli dei conteggi, lo vede chiaramente nel cuore e non trova differenze negli altri cinque tessuti: è quella da usare come controllo negativo.
5. **Le aree simulate si somigliano troppo, e sappiamo perché.** Fra le cinque aree simulate la composizione dei tipi varia da 30 a migliaia di volte meno che fra le vere, la densità delle cellule da 5 a 150 volte meno, perché il simulatore usa una sola composizione e una sola densità per tessuto. Nel linfonodo, per esempio, follicoli e paracorteccia simulati hanno la stessa quota di linfociti B (49 %), mentre nei veri va dal 33 % al 64 %. Le misure dei metodi variano di meno negli stessi tessuti.
6. **L'errore di geometria del cervello è capito e ha una correzione provata.** Un buco di una zona toccava il suo bordo in un solo punto, un caso che la funzione di triangolazione non sa trattare. Un restringimento di un milionesimo di micron, usato solo quando tutto il resto fallisce, lo risolve; un'altra area di controllo esce identica con e senza correzione.
7. **Nel cervello l'eccesso di trascritti estranei sta nei profili, non nello spostamento,** e la composizione dei tipi chiamati da RCTD (troppi neuroni) resta da guardare.

**Cosa devi decidere.**

1. **La prova della semina a compartimenti (lo step che hai già chiesto) deve misurare anche quanto sono sporchi i nuclei, per tenere separate le due cause? Consiglio di sì.**
   - Costa poco: le misure sono pronte e girano in pochi minuti sulle stesse simulazioni.
   - Lo step correggerebbe anche l'errore di geometria con la correzione provata qui.
   - Alternativa: rivedere nello stesso step anche quanto lontano si spostano i trascritti. Sconsiglio: due cambiamenti insieme non permettono di capire quale ha fatto cosa.
   - Alternativa: fare prima lo step sullo spostamento dei trascritti e poi la prova dei compartimenti. Sconsiglio: la disposizione resta la causa più grande, e lo spostamento va prima misurato bene tipo per tipo.

**Comunque**

- Lo spostamento dei trascritti verso i nuclei vicini diventa una voce a sé in coda, con la misura per verificarlo già scritta.
- La variabilità fra aree (una sola composizione per tessuto) e la composizione del cervello restano in coda: ora sappiamo da dove vengono, ma decidere come farle variare è una scelta a parte.
- Per il cuore useremo la misura dell'orientamento nelle prossime prove dei metodi.

Appena approvi, lo step è chiuso.

**Sui tempi di calcolo.** Li ho misurati io sul server. La prima parte dello step: le 30 aree vere con i due metodi circa un'ora e mezza; le 85 simulazioni con i due metodi circa 4 ore con 6 processi; le 30 aree con i compartimenti circa 3 ore con 5 processi. Le indagini di oggi hanno riusato quei risultati: la prova della geometria (tre simulazioni del cervello) circa 18 minuti, le altre letture da 1 a 2 minuti ciascuna; la durata della prima tornata, su linfonodo, cuore e variabilità, non l'ho registrata. La prova dei compartimenti con due semi costerebbe circa 4 ore con 6 processi.

