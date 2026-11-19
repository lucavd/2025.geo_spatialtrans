## In breve: come teniamo da parte dati veri mai usati per verificare il simulatore, e con quali regole lo cambiamo?

**A cosa serviva.** Ogni numero del simulatore viene da una sola sezione per tessuto, e quasi sempre dalle stesse aree su cui poi lo abbiamo verificato: una verifica così conferma per costruzione. Questo step non ha toccato il simulatore. Ha contato da dove viene ogni valore, ha messo da parte dati veri mai usati, e ha scritto le regole per accettare un cambiamento senza adattare il simulatore a sei sezioni.

**Cosa abbiamo trovato.**

1. **380 dei 478 valori che il simulatore usa vengono dalle stesse 30 aree su cui li abbiamo verificati, una sezione per tessuto.**
   - Gli altri vengono dai bordi del tessuto o dalla sezione intera (lo spostamento dei trascritti e le strisce della griglia) oppure sono impostazioni.
   - Nessuno script dei passi precedenti legge aree che non abbiamo censito.
2. **Abbiamo 48 aree vere di verifica in otto tessuti, ma sono poche proprio dove servirebbero di più.**
   - Altri pazienti dello stesso tumore del colon (i pazienti 1 e 5): 7 aree di tumore e 4 di stroma; un secondo tumore del colon, di paziente non noto: 5 di tumore e 3 di stroma; una seconda sezione di cervello di topo: 4 aree.
   - Due tessuti mai visti, il rene umano e il fegato di ratto: 10 aree, per verificare le parti tecniche comuni a tutti i tessuti, come lo spostamento dei trascritti.
   - Intestino di topo, linfonodo e cuore non hanno una seconda sezione pubblica. Restano aree nuove della stessa sezione: 4 nell'intestino, 5 nel cuore, 2 follicoli nel linfonodo (in aree di 500 µm, perché lontano dalle aree già usate c'è soprattutto grasso), nessuna di paracorticale. Lo stroma del tumore originale non ne ha nessuna.
   - Tutti i file scaricati sono integri, e cambiare versione del programma di 10x tra i pazienti sposta i conteggi di circa l'1 %: trascurabile.
3. **Le aree le ho scelte io, ma non come previsto: due regole automatiche sui marcatori hanno fallito la loro prova sulle aree già approvate, e la terza funziona solo sulle sezioni che conosce.**
   - Le aree tenute sono quelle su cui concordano la regola (dove funziona) e due letture a vista indipendenti dell'immagine: la mia e quella di un secondo lettore che non sapeva nulla delle altre. I due lettori concordano su 9 aree su 10.
   - Nessuno dei due è un patologo: fermano gli errori grossolani (grasso, muscolo, mucosa normale mescolata al tumore), non quelli sottili.
4. **Con l'obiettivo scritto come intervallo, oggi il simulatore risulta fuori quasi ovunque.**
   - La regola dice che i dati veri devono cadere dentro la variabilità delle aree simulate. Sui dati già usati solo 13 aree vere su 240 ci cadono: le aree simulate di un tessuto si somigliano molto più delle vere, perché il simulatore usa una sola composizione e una sola densità per tessuto.
   - La semina a compartimenti più che raddoppia la variabilità fra aree simulate, ma resta circa un terzo di quella vera.

**Cosa devi decidere.**

1. **Adottiamo le regole contro il sovradattamento? Consiglio di sì.**
   - In sintesi: i dati di verifica non entrano mai nella taratura e si leggono solo dopo aver scritto le previsioni; un cambiamento si accetta se ha un meccanismo documentato, al più due parametri nuovi, comuni ai tessuti quando sono tecnici, se migliora almeno 4 tessuti su 6 senza peggiorarne nessuno e se regge sui dati di verifica; se non regge non si ritocca; al più tre indagini a posteriori per step; ogni insieme di dati di verifica giudica al più tre cambiamenti, poi diventa dato di taratura.
   - Alternativa: regole più strette, per esempio una sola lettura dei dati di verifica per insieme. Più sicura, ma con così pochi dati veri li consumeremmo dopo il primo cambiamento.
2. **La prova della semina a compartimenti resta il prossimo passo, giudicata subito con queste regole, tranne la parte sulla variabilità fra aree? Consiglio di sì.**
   - La parte «i dati veri dentro la variabilità del simulato» oggi boccia tutto per un difetto noto. Applicarla adesso non direbbe nulla sui compartimenti.
   - La si applica dal passo che farà variare composizione e densità fra le aree, da mettere subito dopo i compartimenti.
   - Alternativa: fare prima la variabilità fra aree e poi i compartimenti. È più pulito per le regole, ma rimanda il limite più grande, la disposizione dei tipi.
3. **Chiudiamo lo step anche se le prove delle due prime regole sulle aree sono fallite come scritte? Consiglio di sì.**
   - Le prove hanno fatto il loro lavoro: hanno impedito di usare regole che sbagliavano. Nessun dato e nessun risultato dipende da quelle regole.
   - Alternativa: tenere aperto lo step e cercare una regola automatica che passi. Sarebbe il quarto tentativo sulle stesse 30 aree: è proprio il sovradattamento che vogliamo evitare.

**Comunque**

- Quando proveremo un cambiamento, le previsioni sui dati di verifica si scriveranno prima di leggerli, e ogni lettura si annoterà in un registro.
- Nei confronti con i pazienti 1 e 5 useremo solo i geni che le due versioni del programma di 10x hanno in comune: la revisione indipendente ha trovato che l'1 % di differenza è la somma di due effetti opposti, geni in più in una versione e qualche trascritto in meno sugli altri.
- Per intestino, linfonodo e cuore scriveremo «verificato solo dentro la sezione»; il secondo tumore conta come altro paziente solo se un giorno sapremo chi è.
- Le aree poche (stroma originale, paracorticale del linfonodo, cervello) restano in coda insieme al modo di averne di più.

Appena decidi questi 3 punti, lo step è chiuso.

**Sui tempi di calcolo.** Li ho misurati io sul server: scaricare e controllare i 138 GB dei dati di verifica ha richiesto circa un'ora e un quarto; tutte le altre prove di questo step da pochi minuti a un quarto d'ora ciascuna. La prova dei compartimenti con i dati di verifica costerebbe, secondo la mia stima dai tempi di F3.1 (circa 3,5 simulazioni con i due metodi per processore all'ora), circa 9 ore con 6 processi o poco più di 2 ore con 24, per 48 aree, due semi e due modi di semina; in più circa due ore e mezza per far girare i due metodi sulle 48 aree vere (in F3.1 le 30 aree vere ne avevano richiesta una e mezza).

