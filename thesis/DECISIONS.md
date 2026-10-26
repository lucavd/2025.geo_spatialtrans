# DECISIONS — tesi Pierangela (thesis/expression-sources)

Una decisione vale solo se è scritta qui, con la data. Si aggiunge in fondo, non si cancella: se una decisione cambia, si scrive una riga nuova che la sostituisce e cita la precedente.

| Data | Decisione | Chi |
|---|---|---|
| 2026-10-14 | Tema, titolo e obiettivi come in PROTOCOLLO.md, sez. 1 | Luca |
| 2026-10-14 | Periodo 19/10–31/12/2026 a tempo pieno; tesi scritta e consegnata entro il 31/12/2026 | Luca |
| 2026-10-14 | Bracci E, P, M, G, C−, E′; braccio G opzionale, condizionato al gate G0 entro il 30/10 | Luca |
| 2026-10-14 | Tessuti: colon/tumore del colon, tessuto linfoide, cervello | Luca |
| 2026-10-14 | Metodi: RCTD, CARD, cell2location (GPU), NNLS; SPOTlight facoltativo | Luca |
| 2026-10-14 | Rumore binomiale negativo unico per tutti i bracci; pseudo-spot con composizione Dirichlet; nessun uso del cell layer | Luca |
| 2026-10-14 | Split per donatore o campione, mai per cellula | Luca |
| 2026-10-14 | Account Synthesize Bio personale di Pierangela; richiesta di permesso inviata da lei, Luca in copia | Luca |
| 2026-10-14 | Branch dedicato con push riservato; nessun merge durante la tesi | Luca |
| 2026-10-14 | Lavoro sul server con utente proprio (`pierangela`, non admin), clone proprio in `/mnt/micron/geo_spatialtrans/tesi_pp1209/repo`, token GitHub personale (classic, scope `public_repo`) | Luca |
| 2026-10-14 | Regola di taglio: SPOTlight → braccio G → tessuto cervello; mai E, P, M, C−, E′, colon, linfoide | Luca |
| 2026-10-14 | Fino alla VPN lavoro sul PC con clone parziale del branch; percorsi dati da `TESI_DATA`; dati solo da script con md5, mai copiati dal PC al server | Luca |
| 2026-10-14 | Riserva: senza VPN entro il 30/10, ADEMP approvato con pilota su CPU (RCTD, CARD, NNLS); repliche di cell2location fissate nei primi giorni sul server con emendamento datato prima di ogni risultato | Luca |
| 2026-10-26 | Dataset: cervello = PsychAD (10x 3′ v3, single-nucleus); colon = Pelka 2021, sottoinsieme 10x 3′ v2; linfoide = Wells 2025, linfonodo (10x 5′ v2). Un dataset per tessuto: training, held-out ed E′ sono donatori diversi dello stesso dataset. 5′ e single-nucleus accettati perché i confronti sono dentro il tessuto. Limite da dichiarare: il 5′ potrebbe svantaggiare GEM-1 per motivi tecnici | Pierangela, approvato da Luca |
| 2026-10-26 | Colon: lo split è per donatore, non per campione; il campione normale e quello tumorale di una persona finiscono sempre nello stesso split (evita pseudoreplicazione e leakage). Gli pseudo-spot si costruiscono dentro un singolo campione, mai mescolando sano e tumore. Precisa la riga del 2026-10-14 "split per donatore o campione" | Luca |
| 2026-10-26 | E′ è il riferimento della variabilità tra donatori ("buono quanto un altro donatore reale"), non il tetto: una fonte che supera E′ non indica leakage. Sostituisce il ruolo "controllo positivo: il tetto" di PROTOCOLLO §5.2 e la controprova 2 di §5.6 | Pierangela, approvato da Luca |
| 2026-10-26 | Tetto (noise ceiling): per ogni tessuto e tipo cellulare, le cellule del donatore held-out si dividono a caso in due metà (divisioni ripetute); l'affidabilità tra le metà si corregge con Spearman–Brown (2r/(1+r)) e la sua radice è la correlazione massima che una fonte qualunque può raggiungere con l'held-out. Stessa correlazione e stessi geni dell'Obiettivo 1; dettagli nell'ADEMP | Luca |
| 2026-10-26 | Controllo del leakage sulla provenienza dei dati: test automatico che i donatori dei vari split siano disgiunti; per GEM-1 solo donatori post-cutoff, con la fonte della data di cutoff annotata | Luca |
| 2026-10-26 | Metodi: RCTD, cell2location (GPU), spatialDWLS (pacchetto standalone `spatialDWLS`, drieslab, senza Giotto), MuSiC (seconda baseline), NNLS (baseline). CARD escluso perché sfrutta la struttura spaziale, che gli pseudo-spot per costruzione non hanno: con coordinate finte sarebbe penalizzato ingiustamente; non per prestazioni (tra i migliori in Li et al. 2023, Nat Commun, doi:10.1038/s41467-023-37168-7). SPOTlight tolto. Regola di taglio: braccio G → tessuto cervello. Sostituisce le righe "Metodi" e "Regola di taglio" del 2026-10-14 | Luca |
| 2026-10-26 | Braccio G: si valuta solo sui tessuti con donatori held-out post-cutoff nello stesso dataset usato dagli altri bracci (oggi cervello e linfoide; il colon va ricontrollato con il criterio di cutoff della riga sotto, perché Pelka risulta su CELLxGENE solo dal 2025-08-05); il tessuto senza donatori post-cutoff è dichiarato "non valutato". Se un dataset post-cutoff diverso ha almeno 10 donatori, si valuta di farlo diventare il dataset principale di quel tessuto. Decisione finale al gate G0 del 30/10 | Luca |
| 2026-10-26 | Benchmark indipendente per le letture (PROTOCOLLO §10): Spotless (Sang-aram et al. 2024, eLife, doi:10.7554/eLife.88431) | Luca |
| 2026-10-26 | Cutoff di GEM-1-sc: versione di luglio 2024 di CELLxGENE (più 10 dataset scPerturb), dal preprint GEM-1 (doi:10.1101/2025.09.08.674753, Methods). Un dataset è "post-cutoff" se non era in CELLxGENE a luglio 2024. Da confermare per il modello servito dall'API con la risposta di Synthesize Bio (gate G0) | Pierangela |
