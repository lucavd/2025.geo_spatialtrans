
Commit dei cancelli: `e03f6bd` (perimetro, stime di preparazione, regola R-1, cancelli e [Il consolidamento dello Step 2 si fa in due step: S2.8 per il codice, S2.8b per i documenti e l'unione nel ramo di sviluppo (DEC-078)](../docs/decisions/DEC-078.md)). Ok di Luca sulle previsioni P-1...P-3, sulla regola R-1 con la soglia di 1.5, sulla controprova e sul calcolo di circa 6 ore-processo: 2026-11-11, «Ok, come proposte (Recommended)».

Semi: quelli di P-1 di S2.7, k = 0, 1, 2 per area: funzione a pezzi dalla tabella `tasks_P1.json` di S2.7 (variante «new»); motore di S2.5 k = 0, 1 dai json di `S2.5/panel/main_r2/` e k = 2 da `tasks_P1.json` (variante «old»). Nessun seme nuovo: tutte le esecuzioni hanno un'uscita salvata in S2.7 o in S2.5 con cui confrontarsi (G-2). Unità indipendente per l'incertezza: il seme (3 per area e configurazione); stima ± 2 errori standard fra semi. Uscite in `/mnt/micron/geo_spatialtrans/S2.8/`, tabelle in `results/S2.8/`.

### Cancelli (al massimo 5)

| id | controllo | tolleranza | comando |
|---|---|---|---|
| G-1 | **batteria** `R/testing/test_S2.8.R`: (a) l'impronta contiene l'impronta del codice e quella dei profili; una ripresa con codice e argomenti uguali riusa i file, una ripresa dopo una modifica di una funzione dello Step 2 o con profili diversi si ferma con un errore chiaro; (b) `truth$layers$stripes$M` = `metadata$overgeneration` = M*, il massimo grezzo è nel campo nuovo, `bins$c` = min(`bins$c_raw`, M*); (c) `metadata$package_version` viene da `DESCRIPTION`; (d) i test che la mappa trova mancanti (dichiarati prima dell'esecuzione ufficiale). **Mutanti** (3): impronta senza il codice; impronta senza i profili; M lasciato grezzo. Tutti rilevati. **Batterie di S2.1...S2.7**: tutte PASS, con le asserzioni di S2.7 G-4 (21, 36, 39, 37, 40, 33) e di S2.7 G-1 (78); cambiano solo asserzioni che leggono i nomi rinominati, dichiarate | identità esatte; mutanti 3/3 | `R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla R/testing/test_S2.8.R`; batterie |
| G-2 | **i numeri non cambiano**: (a) E504, 18 esecuzioni: impronte md5 di `counts_2um`, `counts_8um`, `bins_8um`, `overlap`, `overlap_8um`, `origin`, `origin_8um`, `background`, `background_8um`, `arrived_2um`, `cells`, `genes` uguali a quelle di S2.7 (`p1/new/*.md5`); `bins` e `layers` uguali dopo aver tolto i campi rinominati; (b) E1008 ed EInf uguali a E504; (c) E0: file binari dei conteggi a 2 e 8 µm uguali a quelli salvati (S2.7 `p1/old` per k = 2, S2.5 `main_r2` per k = 0, 1) | md5 identici | `tools/S2.8_identity.R` |
| G-3 | **mappa completa**: ogni sezione 1-17 del documento di progetto e ogni argomento esportato di `R/05b*.R` ha una riga; ogni scarto ha una classe (difetto corretto con test, documento annotato con fonte, voce di BACKLOG); 0 righe senza classe; documento v0.6 committato; controllo OKF senza errori | esatta | `results/S2.8/S2.8_design_map.csv`; `python3 -B tools/okf.py check --bundle docs` |
| G-4 | **integrità**: md5 dei file letti prima e dopo (livelli delle sei aree di S1.6, griglie di S2.3, profili, uscite di S2.5 e S2.7 usate nei confronti); `git diff 2ddb679 -- R/` limitato a `R/05b5_simulate_visium_hd.R`, `R/05b6_simulate_visium_hd_pieces.R` e ai file dichiarati in questo file per i difetti trovati dalla mappa; `R/04b*.R` identici | identità | `tools/S2.8_md5.sh before / after` |
| G-5 | **revisione avversariale**: sotto-agenti con codice proprio ricalcolano (a) i confronti di G-2 dai file; (b) gli esiti di P-1...P-3 e la regola R-1 dai json; (c) le descrittive di BL-237; (d) leggono le modifiche del codice contro questo file; (e) ricontrollano un campione della mappa contro il documento e il codice | quantità deterministiche identiche; ogni discrepanza è un dubbio | sotto-agenti, report |

### Previsioni

Intervallo di rilevanza per tutte: il bilancio del documento di progetto, sezione 12 («un'area di 1 mm² in meno di 5 minuti e 8 GB»), già usato in S2.7 (P-3). Uno scarto oltre il bilancio cambia ciò che la guida dirà su come simulare le aree piccole.

| id | quantità | tessuti | atteso | intervallo di rilevanza | fonte dell'intervallo | 2 errori standard tipici |
|---|---|---|---|---|---|---|
| P-1 | tempo di orologio della funzione a pezzi con il lato predefinito (E504) su un'area di 1 mm², media dei 3 semi | A1-A6 (6 giudizi) | confermata nello stroma e nel cuore; tumore del colon smentita (S2.7 sotto carico: 385 s); intestino, linfonodo, cervello confermate o indeterminate (S2.7 sotto carico: 281-292 s) | [0, 300 s] | bilancio, design sezione 12 | circa 10-20 s (S2.7: da 278 a 291 s fra 9 semi nell'intestino) |
| P-2 | picco di memoria di E504 (massimo fra processo principale e pezzi, massimo sui semi) | A1-A6 (6 giudizi) | confermata 6/6 (S2.7: al più 3.9 GB) | [0, 8 GB] | bilancio, design sezione 12 | piccolo (picchi di S2.7 stabili fra semi) |
| P-3 | picco di memoria della funzione a pezzi in un pezzo solo (EInf) su un'area di 1 mm² | A1-A6 (6 giudizi) | smentita nel tumore e nel cervello (stima da E-1: circa 9-10 GB), confermata negli altri quattro | [0, 8 GB] | bilancio, design sezione 12 | piccolo |

Il cuore qui non è un controllo negativo del Voronoi isotropo: nessuna misura dipende dalla forma delle cellule.

### Controprova

- **CP-1, la mappa legge davvero il documento contro il codice?** Il documento v0.5 è stato scritto prima di S2.6 e S2.7: la sezione 1 non nomina `simulate_visium_hd_pieces()`, la sezione 12 prescrive ancora tessere di 500 µm con un bordo di 280 µm, la sezione 15 si ferma a S2.5. La mappa deve segnalare come indietro almeno queste tre sezioni. Se non le segnala, non distingue un documento vecchio da uno giusto e G-3 non vale nulla.
- **CP-2, il confronto delle impronte vede una differenza quando c'è?** Le esecuzioni E504 con k = 0 rifatte con un tetto diverso della sovra-generazione (`ceiling_q = 1e-4` invece di 1e-5) devono dare un md5 di `counts_2um` diverso da S2.7 in 6 tessuti su 6.

### Misure (senza intervallo)

- Tempi e picchi di E0, E1008, EInf; rapporti della regola R-1; tempo per fase della funzione a pezzi.
- **Descrittiva di BL-237** (BL-248): area A1_r1 (intestino; il tessuto comincia a 5.48 µm dal lato sinistro della finestra), funzione a pezzi E504, `frame_open = TRUE`, stessi semi k = 0, 1, 2, con `frame` = finestra di R2 dell'area invece del riquadro di regioni e territori (predefinito, già in E504): rapporto fra le due versioni dell'UMI medio per quadratino di 2 µm in tutto il tessuto e in fasce di 0-8, 8-16, 16-32, 32-64 µm dal lato sinistro del tessuto; UMI dei quadratini di vetro fra il lato della finestra e il tessuto; valori veri degli stessi quadratini dove S2.5 li ha.

### Stima dei tempi (dell'agente, dai tempi misurati in S2.7)

E504 18 esecuzioni circa 1.2 ore-processo (S2.7: 108-385 s); E1008 ed EInf circa 1 ora ciascuna; E0 con la verità geometrica circa 1 ora (S2.7: 29-238 s, più la verità geometrica); CP-2 6 esecuzioni circa 0.4; BL-237 3 esecuzioni circa 0.3; batterie di S2.1...S2.8 circa 1. In tutto **circa 6 ore-processo**, circa 1 ora e mezza di orologio con 6 processi. La mappa del documento è lettura (sotto-agente), senza calcolo.

