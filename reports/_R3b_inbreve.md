# In breve

## In breve: nel linfonodo e nel cuore, da quale insieme di nuclei veri prendiamo posizioni, aree e vicini?

**A cosa serviva.** Nel linfonodo e nel cuore il simulatore prendeva i centri dei nuclei da Space Ranger, ma aree e vicini da StarDist: tre misure da due maschere diverse. Abbiamo confrontato con i tuoi conteggi manuali due maschere uniche: «i nuclei di Space Ranger con la forma di StarDist» e «solo StarDist». Volevamo anche capire perché nel linfonodo i nuclei simulati si toccano più di quanto vedi tu.

**Cosa abbiamo trovato.**

1. **La regola sceglie in tutti e due i tessuti i nuclei di Space Ranger con la forma di StarDist: ne conta quanti ne conti tu, mentre StarDist da solo ne perde troppi.**
   - Nelle finestre di linfonodo questa maschera trova 416 nuclei contro i tuoi 404; StarDist da solo 314, quasi uno su quattro in meno. Nel cuore 262 e 202 contro i tuoi 242.
   - Il conteggio giusto viene in parte da una compensazione. Dei 102 nuclei che restano con la forma di Space Ranger (quelli che StarDist non vede) 65 sono nuclei veri, spesso pallidi, e 37 non corrispondono a nessun tuo punto. Nel cuore 47 veri su 67.
2. **L'eccesso di contatti simulati nel linfonodo non viene dai centri dei nuclei: con i tuoi punti come centri resta lo stesso.**
   - Con i tuoi punti i nuclei simulati a contatto sono 0.72, con la maschera scelta 0.69, contro il tuo 0.55; tutti dentro il tuo margine, che arriva a 0.75.
   - Lo 0.41 di StarDist dipendeva dai nuclei che StarDist non vede. La spiegazione che avevamo dato quando abbiamo portato il nucleo ritagliato nel simulatore (centri di Space Ranger troppo fitti) era sbagliata.
   - Nel cuore tutte le versioni stanno nel margine (0.05 contro il tuo 0.02).
3. **I nuclei rimasti con la forma di Space Ranger hanno aree e contatti poco affidabili.**
   - Sono più di un nucleo su quattro in entrambi i tessuti, con la forma a scalini di 2 µm. Dove lo stesso nucleo è visto anche da StarDist, il poligono di Space Ranger è più grande di un terzo, con grandi oscillazioni.
   - Se li usiamo per la taglia, la variabilità della taglia dei nuclei del linfonodo cresce di più di un terzo; se la stimiamo sui soli nuclei con la forma di StarDist resta quella di oggi. I contatti simulati nel linfonodo passano da 0.69 a 0.74, al bordo del tuo margine.
   - Sui nuclei che hai giudicato nel linfonodo, la maschera segna 44 contatti su 59 contro i tuoi 32: quasi due contatti su tre passano da un poligono di Space Ranger. Il riferimento del simulatore per i contatti resta comunque il tuo conteggio.
4. **Con la nuova maschera i territori non cambiano; cambiano il rapporto fra nucleo e cellula nel linfonodo e il tetto del nucleo nel cuore.**
   - Area e forma dei territori cambiano meno del 2 %.
   - Nel linfonodo il nucleo occupa in mediana il 48 % della cellula invece del 61 %: il vecchio valore usava le aree dei poligoni di Space Ranger. Nel cuore il tetto al nucleo scende dal 18 % al 15 % della cellula.
5. **Un controllo di integrità è fallito per come l'avevo scritto, non per un difetto della maschera.**
   - Avevo chiesto che ogni nucleo della maschera avesse il centro nel tessuto valido. Sulla maschera salvata 59 nuclei su 145 000 hanno il centro fuori, come nelle maschere originali di Space Ranger (60) e di StarDist (107).
   - Il filtro del tessuto si applica dopo: nessuno di quei nuclei entra nei territori o nei conteggi. I due revisori indipendenti hanno rifatto la maschera e tutti i numeri e trovano gli stessi valori.

**Cosa devi decidere.**

1. **Adottiamo per linfonodo e cuore i nuclei di Space Ranger con la forma di StarDist, ma stimiamo la taglia dei nuclei solo su quelli con la forma di StarDist? Consiglio di sì.**
   - Posizioni, territori e vicini vengono da un'unica maschera, che conta i nuclei come te. La taglia si misura solo dove la forma è affidabile, e il modello della taglia resta quello di oggi.
   - Alternativa: la taglia da tutti i nuclei della maschera. Più coerente sulla carta, ma un nucleo su quattro porta un'area a scalini e circa un terzo di questi non è un nucleo; la variabilità della taglia nel linfonodo crescerebbe di più di un terzo.
   - Alternativa: solo StarDist. Anche i suoi contatti simulati stanno nel margine (0.41), ma nel linfonodo ha un nucleo su quattro in meno, quindi territori più grandi di circa un terzo.
2. **Accettiamo che il controllo di integrità fallito sia un errore di come l'avevo scritto, e non un difetto della maschera? Consiglio di sì.**
   - Le maschere originali da cui partiamo hanno lo stesso comportamento, e il filtro del tessuto esclude quei nuclei da ogni calcolo.
   - Alternativa: filtrare già la maschera salvata e ripetere i controlli. Non cambierebbe nessun numero.

**Comunque**, appena decidi aggiorneremo i riferimenti di linfonodo e cuore (territori, taglia, tetto, rapporto fra nucleo e cellula) e la maschera di questi due tessuti nei banchi. I nuclei con la forma di Space Ranger che non sono nuclei veri li riguarderemo più avanti, per esempio tenendo solo quelli con un segnale di ematossilina. Nel linfonodo i contatti simulati stanno nella parte alta del tuo margine anche con i tuoi punti come centri: per stringere il confronto servirebbero più finestre contate a mano, e lo teniamo da parte. La forma vera delle cellule del cuore passa a quando scriveremo la composizione dei tessuti.

Appena decidi questi 2 punti, lo step è chiuso.
