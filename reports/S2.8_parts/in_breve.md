
### In breve: il codice della griglia di Visium HD fa ciò che dice il suo documento di progetto, e quale funzione deve essere l'ingresso del simulatore?

**A cosa serviva.** Prima di scrivere la guida d'uso e di unire lo Step 2 al ramo di sviluppo, dovevamo rivedere tutto il codice della griglia di Visium HD contro il suo documento di progetto. Dovevamo anche scegliere quale delle due funzioni di oggi diventa l'ingresso del simulatore: quella che simula l'area in un pezzo solo, scritta prima, o quella a pezzi, che regge anche il vetrino intero. Il modello non doveva cambiare: nessun conteggio simulato doveva muoversi.

**Cosa abbiamo trovato.**

1. **L'ingresso del simulatore sarà la funzione a pezzi: anche su un'area di 1 mm² costa poco più di quella in un pezzo solo.**
   - Lasciata lavorare in un solo pezzo, su 1 mm² impiega da un 3 % in meno a un 23 % in più del tempo della funzione vecchia, nei sei tessuti, e usa al più 6,2 GB.
   - Con i pezzi da mezzo millimetro, il valore predefinito, impiega da 85 secondi (cuore) a 4 minuti e mezzo (tumore del colon) e meno di 3,5 GB: dentro il bilancio in tutti i tessuti. La volta scorsa sembrava fuori bilancio per due ragioni: il server era carico, e la funzione vecchia riceveva già fatta la parte più lunga del calcolo.
   - La funzione vecchia resta nel pacchetto come riferimento per le prove.
2. **Nessun conteggio è cambiato.**
   - Abbiamo rifatto le simulazioni delle volte scorse con gli stessi semi. I file sono identici bit per bit, con la funzione a pezzi (pezzi da mezzo millimetro, da un millimetro o un pezzo solo) e con la funzione vecchia, anche dopo le ultime correzioni.
   - La prova vede una differenza quando c'è: cambiando di poco il tetto ai trascritti generati, i file cambiano in tutti i sei tessuti.
3. **Il documento di progetto era rimasto indietro in molti punti, il codice in pochi.**
   - Su 187 affermazioni del documento, 124 coincidevano con il codice. 37 erano superate da decisioni prese dopo, soprattutto la simulazione a pezzi, che il documento non conosceva: le ho annotate tutte con la loro fonte, senza cambiare nessuna regola.
   - 19 affermazioni erano giuste ma senza un test, quasi tutte valori dei parametri. Ora i parametri di tutti i tessuti sono controllati alla precisione scritta nel documento.
   - Il codice si scostava dal documento in 5 punti. Quattro erano già noti e ora sono chiusi o restano in coda; il quinto era un valore del fondo sul vetro che una funzione restituisce ancora ma che il simulatore non usa più: ora è scritto chiaro e provato.
4. **Una simulazione interrotta, se la si riprende, non può più mescolare risultati vecchi e nuovi.**
   - La funzione a pezzi salva un file per pezzo e, rilanciata nella stessa cartella, riusa quelli già fatti. Prima li riusava anche se nel frattempo erano cambiati il codice o i profili di espressione.
   - Ora si ferma e dice che cosa è cambiato. La revisione indipendente ha trovato un caso che mi era sfuggito: se cambiavano le cellule, per esempio i loro tipi, la ripresa restituiva in silenzio i conteggi vecchi. Ho corretto e provato anche questo.
5. **Su un'area ritagliata, il lato della finestra conta solo su una striscia sottile di vetro.**
   - Nell'intestino, dare al simulatore la finestra vera del ritaglio invece del contorno del tessuto cambia meno del 2 % i conteggi del tessuto e del vetro vicino.
   - Cambia solo la striscia di vetro larga 5 µm lungo il lato della finestra: da 0,4 a 1 UMI per quadratino. Lo scriveremo nella guida.

**Cosa devi decidere.** Niente: la scelta dell'ingresso l'ha presa la regola che avevi approvato prima della misura.

**Comunque**

- Il prossimo lavoro è sui documenti: la guida d'uso con la funzione a pezzi come ingresso, il README, il materiale per Daniele e l'unione nel ramo di sviluppo. All'apertura ti chiederò se unire anche il ramo principale, fermo al 2025.
- Su aree piccole il pezzo solo è più veloce dei pezzi da mezzo millimetro (da 1,1 a 1,2 volte). Far scegliere da sé al simulatore la grandezza dei pezzi è un'ottimizzazione che faremo insieme a quella del livello delle cellule, già in coda.
- Quando torneremo sul codice dei passi della griglia faremo tre cose rimaste annotate: verificare la forma delle strisce, provare che la posizione dei trascritti non dipende dal gene, scrivere i test di alcuni argomenti interni.

Non c'è nulla da decidere: se non hai dubbi, lo step è chiuso.

**Sui tempi di calcolo.** Li ho misurati io sul server: 81 simulazioni in 3,9 ore-processo (la stima era 6), 44 minuti di orologio con 6 processi, temperatura al massimo 86 °C. A parte, le batterie dei test, rifatte tre volte dopo le correzioni.

