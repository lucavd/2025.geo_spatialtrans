### In breve: il simulatore completo produce conteggi come quelli dei sei tessuti veri?

**A cosa serviva.** Questo step mette insieme tutti i passi del simulatore di Visium HD in una sola funzione: dalle cellule ai conteggi nei quadratini da 2 e da 8 µm, con la verità di riferimento per mettere alla prova i metodi. Poi confronta le 30 aree simulate con le stesse 30 aree vere su sette misure: trascritti per quadratino, quadratini vuoti, geni rilevati (a 2 e a 8 µm), contrasto fra nucleo e citoplasma attorno. Come soglia usa i margini corretti che hai adottato dopo lo step sulle aree sovrapposte, e confronta sugli stessi quadratini dei due lati, come avevi scelto.

**Cosa abbiamo trovato.**

1. **Dopo una sola ritaratura della profondità, trascritti e geni per quadratino tornano sul vero in tutti e sei i tessuti, a 2 e a 8 µm.**
   - Prima della ritaratura l'intestino aveva 14 trascritti su 100 in meno del vero e il cervello 12 in meno. La causa era nota: il fattore di profondità era stato tarato senza i trascritti che lo spostamento porta sul vetro, e la verifica era rinviata a qui.
   - Il fattore nuovo cambia poco nello stroma (1 % in meno) e di più nell'intestino (17 % in più). Che i trascritti tornino giusti è in gran parte per costruzione, perché il fattore si ritara sulle stesse aree; i geni rilevati e i blocchi da 8 µm invece non sono tarati, e tornano anche loro.
   - Il confronto ha potere: senza fattore di profondità il simulato risulta fuori dal vero nel linfonodo, nel cervello e nel cuore.
2. **Restano tre scarti, tutti già visti prima o attesi.**
   - Nell'intestino e nel tumore del colon il simulato ha da un diciassettesimo a un tredicesimo dei quadratini vuoti veri (già visto quando abbiamo aggiunto le strisce).
   - Nel tumore del colon il contrasto fra nucleo e citoplasma è più alto del vero: 0,90 contro 0,82-0,88. Con il margine del ±10 % che usavamo prima sarebbe giusto.
   - Nel cuore i conteggi veri hanno una direzione, quella delle fibre (indice 0,25), mentre il simulato non ne ha (0,05) perché le cellule sono rotonde. Era il controllo negativo previsto: la misura vede il limite.
3. **Il simulatore completo sta nei 5 minuti e negli 8 GB per un'area di 1 mm², anche nel tumore del colon, il tessuto più carico.**
   - A macchina libera il tumore impiega da 2 minuti e mezzo a 4 minuti, con un picco da 6 a 7,2 GB. Con i blocchi di cellule di prima un'area del tumore arrivava a 8,9 GB: li ho dimezzati, senza perdere tempo.
   - Sotto 7 GB non si scende per un'area così carica: un'area di cattura intera andrà fatta a pezzi.
4. **Il confronto dipende molto da quali quadratini si chiamano tessuto.** Sugli stessi quadratini dei due lati, come hai scelto, l'intestino torna sul vero. Con i quadratini di Space Ranger, che contano più bordo con pochi trascritti, lo stesso simulato risulterebbe sopra il vero nei geni rilevati.
5. **Tre revisori indipendenti hanno ricalcolato tutto con codice proprio e trovato due difetti, corretti, che non cambiano nessun numero.** Riguardavano lo stato del generatore casuale dell'utente e un controllo degli argomenti. Due esiti del cervello prima della ritaratura stanno così vicini al margine che li decide il caso del ricampionamento; dopo la ritaratura sono confermati.

**Cosa devi decidere.**

1. **Teniamo come limiti dichiarati i pochi quadratini vuoti nell'intestino e nel tumore e il contrasto troppo alto nel tumore, da riprendere quando rivedremo la densità delle cellule nei preset? Consiglio di sì.**
   - Le due cose toccano gli stessi tessuti densi, e le cause possibili (buchi piccoli nel tessuto, densità delle cellule, nuclei simulati con meno trascritti del vero) sono tutte nella densità dei preset, che dovevamo già rivedere.
   - Alternativa: uno step subito su questi due scarti. Conseguenza: lavoro sul livello delle cellule prima di chiudere lo Step 2, senza un'ipotesi già pronta da provare.
2. **Teniamo come limite dichiarato che il cuore simulato non ha la direzione delle fibre? Consiglio di sì.**
   - Era già un limite scritto del livello delle cellule; ora sappiamo che si vede nei conteggi, con una misura che lo coglie.
   - Alternativa: cellule allungate e orientate nel cuore. Conseguenza: cambia il modo in cui il livello delle cellule disegna i territori, uno step grande.

**Comunque**: lo Step 2 ha la sua funzione completa e il fattore di profondità nuovo nei preset; la verifica su una sezione intera e un'area di cattura intera a pezzi andranno nello step sulle prestazioni; la correzione dei generatori casuali dei passi precedenti e la variabilità fuori dal nucleo restano negli step già previsti.

Appena decidi questi 2 punti, lo step è chiuso.

**Sui tempi di calcolo.** Misurati da me sul server, con al massimo 6 processi: il pannello principale (60 simulazioni) circa 18 minuti, la controprova (30) circa 10, il pannello dopo la ritaratura (60) circa 26 minuti, con il server carico per i revisori. Per un'area del tumore del colon da sola, a macchina libera, da 149 a 232 secondi.
