
Decisione di Luca: [F3.1: la semina a compartimenti diventerà predefinita dopo una prova; i punti non spiegati si indagano in F3.1 prima di chiuderlo (DEC-088)](../docs/decisions/DEC-088.md). Analisi a posteriori: non cambiano gli esiti di F3.1. Dati: le uscite già prodotte in `/mnt/micron/geo_spatialtrans/F3.1/main/` (vero, simulato U6, profili dal nucleo, compartimenti Q6). Un dubbio si dice **spiegato** se il criterio scritto qui sotto è soddisfatto; altrimenti resta rinviato con quello che si è imparato.

**PH-4, linfonodo: perché RCTD chiama plasmacellule tanti linfociti B simulati?** (BL-260)
- Domanda 1: i profili veri generati per i B simulati (prima di spostamento e cattura) assomigliano alle plasmacellule più dei profili dei nuclei B veri da cui vengono? Metodo: per ogni cellula, tipo di massima verosimiglianza multinomiale contro i profili medi per tipo del riferimento di F3.1 (stesso riferimento di RCTD), sui profili veri del simulato (`profiles.parquet`) e sui conteggi dei nuclei veri con la corona (R6) per i B; quota chiamata plasmacellula, per sottotipo.
- Domanda 2: la confusione nasce prima o dopo lo spostamento dei trascritti? Stessa classificazione sulle cellule ricostruite da bin2cell, sui due lati.
- Domanda 3: la composizione per area (follicoli f1, f2 contro paracorteccia p1-p3) è diversa fra simulato e vero? Tipi veri del simulato e tipi di R6 dei nuclei veri, per area.
- Spiegato se: la quota di B simulati chiamati plasmacellule è già alta sui profili veri generati (almeno la metà dell'eccesso visto con RCTD) e un sottotipo o un passo del generatore la spiega; oppure se compare solo dopo lo spostamento e la mescolanza dei vicini la spiega con i compartimenti.

**PH-5, cuore: una misura orientata vede le cellule allungate?** (BL-263)
- Misura: anisotropia del campo dei tipi di RCTD a 8 µm. Per ogni direzione θ (0, 22.5, ..., 157.5 gradi) e distanza 1-3 blocchi, quota delle coppie di blocchi non scartati con lo stesso primo tipo; anisotropia = massimo su θ meno minimo su θ della quota media a distanza 1-3. Stessa misura sul campo degli UMI per blocco (correlazione direzionale), già usato in S2.5.
- Spiegato (misura utile) se: nel cuore l'anisotropia vera supera quella simulata di più di 2 errori standard fra aree, e negli altri tessuti la differenza è minore che nel cuore (controprova: la misura non deve vedere anisotropia dove non c'è). Altrimenti: le misure di F3.1 e questa non vedono il cuore, e lo si scrive.

**PH-6, perché le aree simulate si somigliano più delle vere?** (BL-262)
- Misura: coefficiente di variazione fra le 5 aree, sui due lati, delle sei misure M-1...M-6 e delle quantità di base: nuclei o cellule per mm² di tessuto, composizione per tipo, UMI per quadratino.
- Spiegato se: la variabilità fra aree delle quantità di base del simulato (densità, composizione) è minore di quella del vero di almeno la metà, e la differenza nelle misure dei metodi segue quella delle quantità di base (stessi tessuti). Il preset usa una sola composizione e una sola densità per tessuto: se è così, la causa è questa.

**PH-7, errore di geometria del cervello con i compartimenti** (BL-259)
- Riprodurre l'errore in isolamento (area A5_r3, livello Q6, seme della tabella di F3.1), estrarre la geometria 168, descriverla (vertici, spigoli più corti, validità, area) e confrontarla con il caso di BL-241; provare in una copia di lavoro, non nel pacchetto, la correzione candidata.
- Spiegato se: l'errore si riproduce, la geometria che lo causa è descritta, e la correzione candidata lo toglie senza cambiare le altre aree (stessi conteggi su un'area di controllo). La correzione nel pacchetto va nello step della prova.



