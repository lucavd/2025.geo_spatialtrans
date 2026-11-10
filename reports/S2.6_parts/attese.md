
Commit dei cancelli: `ffc521d` («S2.6: perimetro, esplorazioni, strumenti, cancelli»). Ok di Luca sulle previsioni P-1...P-8, sulle soglie (2 % alle giunture e nel tessuto della sezione; un quarto per la crescita; bilancio del documento di progetto; 0.8-1.25 per il vetro), sulle regole R-1 e R-2, sulle controprove e sul calcolo di circa 44 ore-processo: 2026-11-10, «Ok, come proposte (Recommended)».

Semi (Mersenne-Twister, `sample.int(1e8, n)`): livelli dei quadrati da 20261209 (ordine: tessuto A1...A6, lato 1 000, 2 000, 3 000 µm, seme k = 0, 1); simulazioni in un pezzo solo della domanda 1 da 20261210 (stesso ordine); pezzi del disegno della domanda 1 da 20261212 (tessuto, pezzo); domanda 2 da 20261211 (tessuto, disegno ONE, P1000, P504, P1000H0, k = 0, 1, 2); livelli delle sezioni da 20261213; simulazioni delle sezioni da 20261214 (sezione, in un pezzo solo / a pezzi); dentro un'esecuzione a pezzi il seme di ogni pezzo è estratto dal seme dell'esecuzione. Bootstrap B = 2000, seme 20261215. Processi: al più 20 a un thread (limite termico, file `WORKERS` e guardia della temperatura come in S2.0), memoria prevista totale al più 180 GB. Uscite in `/mnt/micron/geo_spatialtrans/S2.6/` (un file per esecuzione e per pezzo, rilanciabile); tabelle in `results/S2.6/`.

# Gli strumenti

- `tools/S2.6_pieces.R`: `one` (simulazione in un pezzo solo; per i quadrati `frame_open = TRUE`, per le sezioni intere `FALSE`) e `pieces` (pezzi: area utile U di lato T, riquadro = U allargato del bordo H e tagliato al riquadro dell'area intera; cellule con il centro nel riquadro; territori e nuclei tagliati al riquadro; regioni intersecate; quattro ancore di 1e-4 µm agli angoli del riquadro, così la griglia del pezzo è allineata a quella dell'area intera e la riflessione della seconda coda avviene sui lati del riquadro; `frame_open = TRUE` in ogni pezzo; si tengono i quadratini con il centro in U, più il margine di vetro sui lati dell'area intera; seme del pezzo estratto dal seme dell'esecuzione). Un pezzo senza cellule non si simula (oggi `simulate_stripes()` richiede tessuto): i suoi quadratini restano a zero.
- `tools/S2.6_section_layer.R`: livello cellulare della sezione dalla maschera di Space Ranger a 8 µm (pixel di 8 µm), più quattro ancore agli angoli dell'area di cattura [0, 6 704] µm.
- `tools/S2.6_explore.R pc`: livello cellulare di un quadrato pieno (pixel di 2 µm). Argomenti di `place_cells()` come i livelli U6 di S1.6 in tutti i casi.
- `tools/S2.6_measures.R`: misure delle giunture (sotto).

# Definizioni delle misure

- **Tempo**: tempo di orologio della chiamata (`place_cells()`, `simulate_visium_hd()`, o di ogni pezzo); **memoria**: picco del processo (VmHWM). Per i pezzi la memoria è quella del processo figlio che simula il pezzo. Costo per mm² utile = somma dei tempi dei pezzi / area utile.
- **Giunture** (domanda 2): blocchi di 8 µm della griglia dell'area intera con il centro nel quadrato; distanza di un blocco = distanza del centro dalla giuntura interna più vicina; fasce 0-8, 8-16, 16-40, 40-100, 100-250 µm. Nel pezzo solo le stesse linee sono giunture virtuali. **J-1**: UMI medio per blocco nella fascia; rapporto (media dei 3 semi a pezzi) / (media dei 3 semi in un pezzo solo), errore standard col metodo delta dalle deviazioni standard fra semi. **J-2**: quota dei trascritti catturati nella fascia nati da cellule del bordo (cellule simulate anche in un altro pezzo, con un'altra realizzazione). **J-3**: coppie di blocchi adiacenti da lati opposti di una giuntura: correlazione di Pearson degli UMI e coseno medio dei vettori dei geni.
- **Sezione intera** (domanda 3): blocchi di 8 µm dell'area di cattura; tessuto = maschera di Space Ranger; vetro per distanza dal tessuto (16-64, 64-256, 256-512, oltre 512 µm, come S2.4 P-5, distanza fra centri dei blocchi); esclusi i 50 µm al bordo dell'area di cattura. Rapporto a pezzi / in un pezzo solo dell'UMI medio per blocco, intervallo dal bootstrap delle tessere di 500 µm (come S2.4 P-5), B = 2000. Lato vero (descrittivo): matrice di Space Ranger a 8 µm di tutti i blocchi, stesse fasce.

# Regole di decisione dei dati (scritte prima del calcolo)

**R-1, sezione intera in un pezzo solo.** Si simula in un pezzo solo una sezione il cui picco previsto (0.092 GB per milione di trascritti generati + 0.87 GB per mm² di tessuto, dal tumore a 1 mm² con e senza strisce, E-2 ed E-5; trascritti generati per mm² dai quadrati di 1 mm² di E-5) non supera 120 GB; una alla volta, con un controllo che ferma il processo oltre 200 GB. Previsto: intestino circa 100 GB, linfonodo circa 47, cuore circa 49 (si fanno); cervello circa 154, tumore circa 240 (non si fanno: per loro la sezione intera è solo a pezzi).

**R-2, lo step successivo** (decisione dei dati con la delega di Luca, [Dopo S2.5 si fa prima una prova di prestazione su aree grandi; ordine e perimetro delle sessioni successive li prepara l'agente (DEC-071)](../docs/decisions/DEC-071.md)). In ordine:

1. Se P-8 è smentita per il tempo in almeno una sezione, il livello cellulare va ottimizzato prima (id proposto S1.7).
2. Altrimenti, se fra P-3, P-4 (disegno del documento di progetto), P-5, P-6 e P-7 c'è almeno un esito smentito, lo step successivo è la simulazione a pezzi corretta dentro il pacchetto (S2.7), con le richieste che gli esiti indicano: niente bordo simulato due volte se P-3 o P-4 falliscono; seconda coda trattata sull'area intera se P-6 o P-7 falliscono; giunture corrette se P-5 fallisce.
3. In ogni caso S2.7 deve dare a ogni cellula una sola realizzazione e strisce dell'area intera: oggi non è così per costruzione (J-2, strisce per pezzo). Se nessun esito è smentito, questi due difetti si portano a Luca come scelta (limite dichiarato oppure S2.7).

# Cancelli (al massimo 5)

| id | controllo | tolleranza | comando |
|---|---|---|---|
| G-1 | **batteria degli strumenti** (33 asserzioni): un pezzo unico senza bordo = `simulate_visium_hd()` bit per bit (cuore, tumore: UMI per quadratino, matrice a 8 µm, origine); con bordo, ogni quadratino della griglia dell'area intera coperto una volta, UMI cuciti = somma dei pezzi, 8 µm = somma dei 16 quadratini, semi dei pezzi diversi, cellule doppie presenti con il bordo e assenti senza; ogni cellula con un solo pezzo proprietario, territori tagliati dentro il riquadro, riquadro del pezzo = riquadro di territori, regioni e ancore anche su un tessuto a disco; misure su un caso sintetico con un deficit noto e diverso ai due lati della giuntura. 8 mutanti (ancore tolte: riquadro della riflessione diverso dal riquadro del pezzo dove il tessuto non arriva ai lati, taglio saltato, colonna di quadratini persa, trascritti tutti del proprietario, stesso seme nei pezzi, blocchi sfasati nella cucitura, riflessione spenta nei pezzi, giunture spostate nelle misure), tutti rilevati (prima stesura: 6/8; corretti i due test deboli, sotto «Correzioni prima dei cancelli») | identità esatte; mutanti 8/8 | `R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla R/testing/test_S2.6.R`; `.venv/bin/python tools/S2.6_mutants.py` |
| G-2 | **invarianti di ogni esecuzione**: griglia del pezzo allineata (altrimenti lo strumento si ferma), copertura una volta, UMI cuciti = somma dei pezzi, 8 µm = somma dei 2 µm, nessun NA, conservazione di `simulate_visium_hd()` (generati = catturati + scartati + persi) in ogni pezzo | identità | `tools/S2.6_check.R` |
| G-3 | **riproducibilità**: un'esecuzione a pezzi (cuore, 504 µm, k = 0) rifatta con un numero di processi diverso e misure rifatte | md5 identici (UMI, matrice, tabelle delle misure) | `tools/S2.6_repro.sh` |
| G-4 | **integrità**: md5 dei file letti prima e dopo (maschere e matrici delle cinque sezioni, parametri dei preset); `git diff` vuoto su `R/` salvo `R/testing/test_S2.6.R`; batterie di S2.1...S2.5 invariate | identità | `tools/S2.6_md5.sh before / after`; batterie |
| G-5 | **revisione avversariale**: sotto-agenti con codice proprio ricalcolano (a) J-1, J-2, J-3 di due esecuzioni dai file grezzi; (b) le fasce del vetro di una sezione; (c) tempi e memoria dalle tabelle per esecuzione; (d) gli esiti dalle tabelle; (e) leggono gli strumenti contro questo file | quantità deterministiche entro 1e-9 relativo; ogni discrepanza è un dubbio | sotto-agenti, report |

Correzioni prima dei cancelli (dichiarate): (1) la prima stesura della batteria non vedeva due mutanti: sul disco di prova il tessuto toccava tutti i lati di ogni riquadro (le ancore non cambiavano nulla) e nel caso sintetico il deficit restava nella fascia 0-8 µm anche spostando le giunture; ora due dischi lontani e un deficit asimmetrico. (2) Un territorio tangente al lato del riquadro dava, tagliato, poligono più segmento (errore di GEOS «mixed-dimension» nei pezzi senza bordo): lo strumento tiene solo la parte poligonale. (3) Durante E-5 lo strumento è stato riscritto mentre girava (stroma non eseguito): da qui in poi gli strumenti non si toccano durante un'esecuzione.

# Previsioni

Esito come da protocollo: stima ± 2 errori standard tutta dentro l'intervallo, confermata; tutta fuori, smentita; a cavallo, indeterminata.

| id | quantità | tessuti | atteso | intervallo di rilevanza | fonte dell'intervallo | 2 errori standard tipici |
|---|---|---|---|---|---|---|
| P-1 | livello cellulare, tempo per mm² a 9 mm² / a 1 mm² (quadrati pieni) | A1-A6 | confermata 6/6 (E-2: 0.92 fra 4 e 1 mm²) | [0.8, 1.25] | da approvare: la precisione che serve per pianificare (un tempo stimato sbagliato di un quarto non cambia macchina né numero di processi) | pochi percento (2 semi) |
| P-2 | simulazione in un pezzo solo: (a) tempo per mm² a 9 / a 1 mm²; (b) memoria aggiunta per mm² fra 4 e 9 mm² / fra 1 e 4 mm² (crescita lineare, su cui si basano R-1 e le stime delle sezioni) | A1-A6 | attesa: (a) confermata 6/6; (b) confermata 6/6. Base debole, dichiarata: una sola coppia di aree misurata (tumore, 1 e 4 mm², E-2: tempo 1.01, memoria aggiunta 6.7 GB per mm²); nessuna misura oltre 4 mm² e nessuna per lo stroma; le altre attese sono estrapolazioni lineari | [0.8, 1.25] | come P-1 | pochi percento |
| P-3 | memoria di un pezzo del disegno del documento di progetto (504 µm, bordo 280 µm, tutti gli strati) | A1-A6 | attesa (estrapolazione lineare dai quadrati di 1 mm² di E-5; lo stroma è stimato dal rapporto stroma/tumore delle aree di S2.5, non misurato): confermata intestino (circa 6.3 GB), stroma (5.7), linfonodo (3.3), cuore (3.7); indeterminata cervello (7.7); **smentita tumore** (9.7) | [0, 8 GB] | bilancio del documento di progetto, sezione 12 | pochi percento (4 pezzi) |
| P-4 | costo per mm² utile: (a) disegno del documento di progetto; (b) un pezzo solo a 9 mm² | A1-A6 | attesa (come P-3; stroma stimato): (a) **smentita** intestino (circa 920 s), tumore (1 210), stroma (540), linfonodo (790), cervello (780); indeterminata cuore (250); (b) confermata 5/6, indeterminata tumore (270) | [0, 300 s] | bilancio del documento di progetto, sezione 12 (1 mm² in 5 minuti) | pochi percento |
| P-5 | J-1, a pezzi / pezzo solo, cinque fasce, pezzi di 1 000 e 504 µm con bordo | A1-A6 | confermata in tutte le fasce (E-6: 0.99-1.01); nel cuore possibile indeterminata (poche cellule grandi) | **[0.98, 1.02]** | da approvare: la precisione del lato simulato in S2.5 (1-3 %): uno scarto più piccolo non è distinguibile con le verifiche che abbiamo fatto del simulatore | circa 1 % (E-6) |
| P-6 | vetro della sezione intera, a pezzi / pezzo solo, quattro fasce | A1, A4 (coda lunga 1.6-1.8 mm) | oltre 512 µm **smentita per difetto** (i pezzi senza cellule restano a zero; la coda resta nel riquadro del pezzo); 256-512 µm indeterminata o smentita; 16-256 µm confermata o smentita per eccesso (coda ripiegata vicino al tessuto) | [0.8, 1.25] | Luca, 2026-10-31 (R5 CP-2), usato in S2.3 P-4 e S2.4 P-5 | bootstrap delle tessere |
| P-7 | UMI medio per blocco nel tessuto della sezione intera, a pezzi / pezzo solo | A1, A4, A6 | confermata o indeterminata (la coda ripiegata resta nel tessuto: + 1-1.5 % nell'intestino e nel linfonodo) | [0.98, 1.02] | come P-5 | bootstrap delle tessere |
| P-8 | livello cellulare della sezione intera: (a) tempo per mm² di tessuto; (b) picco di memoria | 5 sezioni | (a) confermata 5/5 (linfonodo 70 s per mm², E-4); (b) **smentita nel linfonodo** (11.1 GB, E-4), confermata nel cuore (1.0 GB), indeterminata o confermata altrove | (a) [0, 300 s]; (b) [0, 8 GB] | bilancio dello Step 1 (C11) e del documento di progetto, sezione 12 | un solo seme: nessun errore Monte Carlo, si dice |

Il cuore non è qui un controllo negativo del Voronoi isotropo: nessuna misura di questo step dipende dalla forma delle cellule.

# Controprove

- **CP-1, la misura delle giunture vede un difetto quando c'è?** Pezzi di 1 000 µm senza bordo: J-1 nella fascia 0-8 µm smentita per difetto 6/6 (E-6: 0.924 nel linfonodo). Se in un tessuto non lo è, lì P-5 non ha potere.
- **CP-2, il confronto del vetro non inventa differenze?** Cuore (coda di 113 µm): vetro a pezzi / pezzo solo confermato in tutte le fasce con vetro. Se è smentito, il confronto del vetro ha un errore suo.

# Misure (senza intervallo)

- J-2 e J-3 per tessuto e disegno; tempo per fase dello Step 1 (tassellazione) e della simulazione; esponente di crescita.
- Sezione intera: tempo totale di processo e di orologio, memoria massima per pezzo, pezzi senza cellule; picco della sezione in un pezzo solo contro la previsione di R-1.
- Descrittive contro la sezione vera (BL-205): UMI medio per blocco nel tessuto, vetro per fascia (simulato / vero), per le sezioni a pezzi e in un pezzo solo.

# Stima dei tempi (dell'agente, dalle esplorazioni)

Domanda 1: circa 9 ore-processo (simulazioni in un pezzo solo 7.8, pezzi del disegno 1.3, livelli cellulari 0.3). Domanda 2: circa 13 ore-processo. Domanda 3: circa 21 ore-processo (livelli cellulari 1, pezzo solo 3.7, a pezzi 16.5). In tutto **circa 44 ore-processo**, poco sopra le 20-40 della preparazione; con 16-20 processi circa 3-4 ore di orologio, limitate dalle simulazioni grandi in un pezzo solo (fino a 100 GB, una alla volta).
