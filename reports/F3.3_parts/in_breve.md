# In breve: il simulatore di oggi regge fuori dalle aree su cui è stato tarato?

**A cosa serviva.** Per ogni tessuto il simulatore è tarato su 5 aree di una sola sezione, le stesse su cui finora lo abbiamo verificato. Qui lo abbiamo misurato, senza toccarlo, su 38 aree vere mai usate: aree nuove delle stesse sezioni, e aree di altri pazienti e di altre sezioni. Abbiamo usato gli stessi due metodi di analisi della prima prova, bin2cell e RCTD. Il margine è quello della prima prova: quanto varia ogni misura fra le aree vere di taratura dello stesso tessuto. Volevamo distinguere due difetti. Il primo: il simulatore va bene solo sulle aree di taratura, cioè è sovradattato. Il secondo: sbaglia ovunque allo stesso modo perché gli manca la variabilità.

**Cosa abbiamo trovato.**

1. **Dentro la stessa sezione non vediamo sovradattamento, ma con così poche aree vedremmo solo scarti grandi.**
   - Sulle 15 aree nuove di intestino, tumore del colon, linfonodo, cervello e cuore, i valori veri restano dove stavano le aree di taratura, e il simulato non peggiora.
   - Con così poche aree, però, uno spostamento più piccolo di circa un margine non lo vedremmo.
   - I difetti trovati nella prima prova restano anche sulle aree nuove: 10 su 12. I due che mancano sono nel cervello, dove c'è una sola area nuova.
2. **Gli altri pazienti e le altre sezioni sono molto diversi dalla sezione di taratura, e lì il simulato peggiora di circa tre margini.**
   - La differenza più grande è la profondità. Nel tumore del colon i trascritti per cellula sono circa 1 100 nel paziente di taratura (il simulato ne ha altrettanti), 640 nel paziente 1, 2 000 nel paziente 5 e 270 nel tumore da 6,5 mm. Nella seconda sezione di cervello sono un terzo della prima.
   - Una delle conferme della prima prova era costruita dalla taratura: i trascritti per cellula tornavano perché il fattore di profondità era tarato su quelle aree. Nello stroma la conferma cade in 2 sezioni nuove su 3.
   - Con 3 pazienti non sappiamo separare quanto dipende dal paziente e quanto dalla sezione. Se aggiungiamo il tumore da 6,5 mm, la variabilità fra sezioni non è trascurabile rispetto a quella fra aree per 4 delle 12 combinazioni di misura e tessuto (trascritti per cellula e quadrati con un solo tipo nel tumore; quadrati con un solo tipo e vicini dello stesso tipo nello stroma).
3. **Nello stroma degli altri pazienti i tipi cellulari sono più mescolati che nel paziente di taratura, cioè più simili al simulato.**
   - Due quadrati vicini hanno lo stesso tipo nell'89 % dei casi nel paziente di taratura, nel 72-75 % nei pazienti 1 e 5, nel 65 % nel simulato.
   - Una parte di quello che chiamavamo difetto dello stroma è quindi propria del paziente di taratura. La semina a compartimenti rende i tipi più ordinati: potrebbe avvicinare il simulato a quel paziente e allontanarlo dagli altri.
4. **I controlli funzionano.**
   - Nel cuore, dove ci aspettiamo che il simulatore fallisca perché le cellule vere sono allungate, la differenza di orientamento resta anche sulle aree nuove (0,11, contro 0,10 dentro campione).
   - Un falso «fuori campione», fatto con le sole aree di taratura, non mostra alcuno spostamento.
   - Come previsto, solo 18 valori veri su 228 cadono dentro la variabilità delle aree simulate: il simulatore ha una sola composizione per tessuto.

**Cosa devi decidere.**

1. **Il prossimo passo resta la prova della semina a compartimenti? Consiglio di sì, giudicandola anche sullo stroma degli altri pazienti.**
   - Nelle stesse sezioni i difetti di disposizione dei tipi si confermano sulle aree nuove, e non c'è sovradattamento visibile: la disposizione resta la causa più grande che si può correggere senza dati nuovi.
   - Il punto 3 dice dove la prova può fallire, e va guardato.
   - Alternativa: prima lo step sulla variabilità fra sezioni e pazienti (profondità e composizione che cambiano da sezione a sezione). È lo scarto più grande fuori campione. Ma con 3 pazienti bisogna prima progettare cosa far variare e di quanto, e il compartimento resterebbe da provare.
2. **Un cambiamento del simulatore si giudica fuori campione rispetto al simulatore di oggi? Consiglio di sì.**
   - Oggi la regola dice: nessuna misura confermata sui dati di taratura può diventare smentita in più di una sezione. Il simulatore di oggi la viola già, per i trascritti per cellula nello stroma, a causa della profondità. Così com'è scritta, la regola boccerebbe ogni cambiamento, compresi i compartimenti, per un difetto che il cambiamento non tocca.
   - Proposta: un cambiamento si boccia se rende smentita fuori campione una misura che oggi lì regge, oppure se peggiora oltre il suo margine una misura già smentita.
   - Alternativa: lasciare la regola com'è. Allora prima bisogna far variare la profondità per sezione, altrimenti nessun cambiamento può passare.

**Comunque**

- Quando progetteremo la variabilità fra sezioni e pazienti partiremo da queste misure: profondità fino a 7 volte diversa fra sezioni dello stesso tumore, stroma più mescolato negli altri pazienti.
- I due tessuti nuovi, rene umano e fegato di ratto, restano sigillati per quando rivedremo lo spostamento dei trascritti. Gli altri gruppi di aree di verifica sono stati letti una volta: ne restano due letture ciascuno.
- Nel linfonodo il modo di contare i geni marcatori di due tipi cambia il confronto fra aree nuove e aree di taratura. Quando quella misura entrerà in una decisione, riporteremo i due conteggi.

Appena decidi questi 2 punti, lo step è chiuso.

**Sui tempi di calcolo.** Li ho misurati io sul server. Le 38 aree di verifica e le 10 aree di taratura del tumore rifatte con i geni comuni hanno richiesto circa 2 ore con 4 processi; la prova di riproducibilità circa un'ora con 2. La prova dei compartimenti costerebbe circa 9 ore con 6 processi o poco più di 2 con 24 (stima di F3.2, dai tempi della prima prova dei metodi). Le misure sulle aree vere di verifica non si rifanno: servono solo le simulazioni nuove e i due metodi su di esse.
