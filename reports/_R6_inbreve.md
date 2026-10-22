# In breve

## In breve: di che tipo sono le cellule vere nei sei tessuti, e i tipi spiegano le taglie dei nuclei?

**A cosa serviva.** Il simulatore oggi piazza un solo tipo di cellula per tessuto, con una sola taglia di nucleo. Sapevamo che non basta: nel cervello di topo i nuclei hanno due taglie, nel tumore del colon la taglia ha una coda troppo lunga. Per dare ai preset più tipi servono riferimenti misurati sui tessuti veri: quanti nuclei di ogni tipo, quanto grandi, che forma. Qui abbiamo dato un tipo a ciascuno dei 240 000 nuclei veri delle 30 aree di riferimento (5 per tessuto): abbiamo sommato i trascritti dei quadratini di 2 µm che cadono sul nucleo e in una corona di 8 µm attorno (quanto Space Ranger espande i nuclei, misurato dai suoi dati), raggruppato i nuclei per somiglianza di espressione e dato un nome a ogni gruppo con i geni marcatori di un atlante pubblico (PanglaoDB), secondo il catalogo di tipi che hai approvato.

**Cosa abbiamo trovato.**

1. **Il criterio con cui avevamo deciso di dare i nomi ai gruppi era sbagliato; lo abbiamo corretto, e il confronto con Space Ranger dice che la correzione è giusta.**
   - Il criterio scritto prima chiamava «endoteliali» le cellule di Paneth dell'intestino (5 260 nuclei, riconoscibili dai loro geni) e lasciava senza nome tre gruppi di cardiomiociti: quando un tipo domina il tessuto, un punteggio relativo fra gruppi si azzera proprio su quel tipo.
   - Il criterio corretto (quanto i 50 geni più caratteristici di ogni gruppo sono arricchiti nei marcatori di ciascun tipo) dà nomi coerenti con i geni in tutti e sei i tessuti.
   - Nel cervello e nel cuore, dove Space Ranger ha una sua annotazione, l'accordo sui nuclei appaiati è 91 % e 79 % con il criterio corretto (88 % e 77 % con quello scritto prima). Spostando le impronte di 20 µm l'accordo scende al livello del caso: la classificazione dipende davvero dal nucleo e non dal tessuto attorno.
   - Per regola, i verdetti delle previsioni sono dati con il criterio scritto prima; i numeri che seguono sono quelli del criterio corretto, che proponiamo come riferimento.
2. **Nel cervello di topo i nuclei dei neuroni sono circa 3 volte più grandi di quelli della glia: la previsione è confermata con il criterio corretto (con quello scritto prima resta al limite), ma il divario è molto minore di quello che avevamo letto nei due picchi.**
   - Neuroni 111 µm², glia 34 µm² (astrociti 50, oligodendrociti 31, microglia 23).
   - I neuroni sono due nuclei su tre fra quelli che Cellpose trova; la glia è sottostimata perché Cellpose perde i nuclei piccoli, come sapevamo.
   - La quota di oligodendrociti va dall'8 % al 36 % fra le 5 aree: è la sostanza bianca, e spiega gran parte della variabilità fra le aree del cervello che avevamo visto.
3. **I tipi spiegano la taglia del nucleo nell'intestino, nel tumore e nello stroma del colon (il tipo più grande supera il più piccolo del 24-42 %), ma non nel linfonodo né nel cuore.**
   - Linfonodo: i nuclei dei linfociti B e T differiscono del 5 % (16,8 contro 16,0 µm²): per i linfociti basta una sola taglia.
   - Cuore: la differenza fra tipi è dell'11 %, mentre la mediana dei soli cardiomiociti va da 16 a 28 µm² fra le aree: la variabilità del cuore non viene dai tipi (ipotesi da verificare: orientamento delle fibre nella sezione).
   - Tumore del colon: la dispersione della taglia dentro il solo epitelio (0,56 sul logaritmo dell'area) è uguale a quella di tutte le cellule insieme (0,55): la coda alta non è una mescolanza di tipi, e un modello di taglia per tipo da solo non la toglierà.
4. **I nuclei dei cardiomiociti non sono più allungati degli altri nuclei del cuore: previsione smentita.**
   - Eccentricità mediana 0,82 contro 0,83. La forma del nucleo non dice nulla sulla forma allungata della cellula: per il cuore serve una misura della cellula, non del nucleo.
5. **Le quote dei tipi ipotizzate nel documento di progetto erano quasi tutte sbagliate; ora abbiamo quelle misurate.**
   - Intestino di topo: epitelio 68 % (ipotesi 88 %), stroma e muscolo liscio 21 %, immunitarie 10 %.
   - Tumore del colon: epitelio 71 % (ipotesi 70 %), fibroblasti 20 %. Il criterio di espressione per distinguere epitelio tumorale e normale non funziona: chiama «normale» l'epitelio delle aree che il patologo ha letto come tumore.
   - Stroma del tumore: fibroblasti 32 %, epitelio 41 %, immunitarie 19 %: le aree scelte come stroma contengono molte ghiandole tumorali.
   - Linfonodo: fra i nuclei con un tipo, B 61 % nei follicoli e 42 % nella paracorticale; ma il 47 % dei nuclei resta senza tipo (46 % non risolti, 1 % con meno di 20 trascritti). Con 84 trascritti per nucleo il linfonodo è al limite di questo metodo.
   - Cuore: cardiomiociti 76 % (ipotesi 85 %), fibroblasti 13 %, immunitarie e vascolari 11 %.

**Cosa devi decidere.**

1. **Adottiamo come riferimenti le etichette del criterio corretto, invece di quelle del criterio scritto prima? Consiglio di sì.**
   - Il criterio scritto prima dà riferimenti sbagliati per costruzione nell'intestino (le cellule di Paneth come endoteliali) e nel cuore (un cardiomiocita su quattro senza tipo); quello corretto concorda di più con Space Ranger nei due tessuti dove si può controllare.
   - Alternativa: tenere il criterio scritto prima. I riferimenti dell'intestino e del cuore sarebbero inutilizzabili per i preset.
2. **Nel colon, chiamiamo «tumorale» l'epitelio delle aree che il patologo ha letto come tumore e «normale» quello delle altre, invece di affidarci al criterio di espressione? Consiglio di sì.**
   - Il criterio di espressione (geni delle cripte contro geni degli enterociti) separa gruppi di taglia diversa, non tumore e mucosa sana; le 5 aree del tumore sono tutte tumore per il patologo.
   - Alternativa: cercare marcatori tumorali con fonte. Richiederebbe uno step dedicato, senza un candidato pronto.
3. **Per il linfonodo, teniamo i riferimenti di B e T come provvisori e un'unica taglia per i linfociti? Consiglio di sì.**
   - Quasi metà dei nuclei non ha tipo e i due tipi principali hanno la stessa taglia: per il simulatore il linfonodo può restare a un tipo di taglia, con due composizioni (follicolo, paracorticale) se vorrai.
   - Alternativa: uno step dedicato con più trascritti per nucleo (bin a 8 µm o la matrice per cellula di Space Ranger): circa 2 ore di calcolo sul server con 6 processi, nessuna annotazione manuale; nel caso migliore porterebbe sotto il 10 % i nuclei senza tipo, senza cambiare la conclusione sulla taglia.

**Comunque:**

- quando porteremo i tipi nei preset useremo i riferimenti per tipo di questo step: due taglie nel cervello (neuroni e glia) con una quota di oligodendrociti che cambia fra le aree; nel tumore del colon la coda alta resta dentro l'epitelio e andrà trattata nel modello di taglia, non con i tipi;
- la forma allungata delle cellule del cuore resta da misurare sulla cellula (membrane o citoplasma), non sul nucleo;
- i riferimenti per tipo entrano nella tabella dei riferimenti biologici con la loro fonte, segnati provvisori per il linfonodo e per lo stroma del tumore (aree miste).

Appena decidi questi 3 punti, lo step è chiuso.

**Sui tempi di calcolo.** Stima dell'agente dai tempi misurati in questo step (aggregazione dei bin 2-7 minuti per dataset, clustering 1 minuto per tessuto, analisi con bootstrap circa 2 ore su un processo). L'alternativa del punto 3 ripete aggregazione e clustering sul solo linfonodo con bin a 8 µm: circa 2 ore su lesexp-server con 6 processi, entro i limiti di calcolo del progetto.
