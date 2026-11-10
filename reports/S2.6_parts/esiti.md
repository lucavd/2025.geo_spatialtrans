
Tabelle: `results/S2.6/S2.6_verdicts.csv` (una riga per giudizio), `S2.6_verdict_counts.csv`, `S2.6_q1_runs.csv`, `S2.6_q1_design_pieces.csv`, `S2.6_q2_measures.csv`, `S2.6_q2_pairs.csv`, `S2.6_q2_J1_J2.csv`, `S2.6_q2_J3.csv`, `S2.6_q3_step1.csv`, `S2.6_q3_runs.csv`, `S2.6_section_<A>.csv`, `S2.6_PH3_edge.csv`. Uscite grezze in `/mnt/micron/geo_spatialtrans/S2.6/` (`layers/`, `q1/`, `q2/`, `q3/`, `sections/`, `ph1/`, `ph2/`).

# Cancelli

- **G-1** superato: 33/33 asserzioni, mutanti 8/8 (`results/S2.6/S2.6_G1_tests.txt`, `S2.6_mutants.txt`; prima stesura 6/8, vedi «Correzioni prima dei cancelli»).
- **G-2** superato: 140/140 esecuzioni con tutti gli invarianti (`S2.6_G2_check.csv`).
- **G-3** superato: pezzi di 504 µm del cuore (k = 0) rifatti con 2 processi invece di 4: UMI, origine, matrice a 8 µm e tabelle delle misure identici (`S2.6_G3_repro.txt`).
- **G-4** superato: 69 file letti identici prima e dopo (`S2.6_md5_diff.txt` vuoto); `R/` invariato (`S2.6_gitdiff_R_after.txt` vuoto). Il livello della sezione del linfonodo è un'uscita dello step, non un dato letto: la sua correzione è nelle deviazioni.
- **G-5** superato: tre revisori, 279 confronti, nessuna discrepanza deterministica (sotto «Revisione avversariale»); due correzioni dello strumento delle sezioni, sotto «Correzioni dopo il calcolo».

# Previsioni (giudizi: confermate / indeterminate / smentite; non valutabili a parte)

| id | giudizi | conf. | indet. | smentite | in breve |
|---|---|---|---|---|---|
| P-1 livello cellulare, tempo per mm² 9/1 | 6 | 3 | 0 | 3 | smentita per difetto nello stroma (0.71), cervello (0.49), cuore (0.45): costo fisso di circa 1 s |
| P-2a simulazione, tempo per mm² 9/1 | 6 | 5 | 1 | 0 | tumore indeterminata (0.96 ± 0.22) |
| P-2b memoria aggiunta per mm², 4-9 / 1-4 | 6 | 0 | 6 | 0 | 2 semi, errore standard a 1 grado di libertà |
| P-3 memoria di un pezzo del disegno | 6 | 5 | 1 | 0 | tumore 8.9 ± 1.3 GB (GB = 1024³ byte; in GB decimali 9.5 ± 1.4, smentita) |
| P-4a costo per mm² utile, disegno | 6 | 1 | 0 | 5 | smentita 527-1 325 s; cuore confermata (272 ± 12 s) |
| P-4b costo per mm² utile, pezzo solo a 9 mm² | 6 | 6 | 0 | 0 | 59-278 s |
| P-5 giunture con bordo, cinque fasce, due disegni | 60 | 17 | 42 | 1 | confermate 8/10 nel linfonodo; smentita stroma, 504 µm, 0-8 µm (0.951 ± 0.024, spiegata da PH-2) |
| P-6 vetro della sezione, intestino e linfonodo | 8 | 2 | 1 | 5 | oltre 512 µm 0.55 e 0.49; 256-512 µm linfonodo 0.65; vicino all'intestino 1.30 e 1.51 |
| P-7 tessuto della sezione | 3 | 1 | 2 | 0 | intestino 1.016 [1.011, 1.022], linfonodo 1.015 [1.006, 1.024], cuore 1.001 |
| P-8a livello cellulare della sezione, tempo per mm² | 5 | 5 | 0 | 0 | 1.2-96.5 s per mm² |
| P-8b livello cellulare della sezione, memoria | 5 | 4 | 0 | 1 | linfonodo 10.5 GB |

# Controprove

- **CP-1** (pezzi senza bordo, fascia 0-8 µm): smentita per difetto in intestino (0.919 ± 0.037), stroma (0.929 ± 0.047), linfonodo (0.939 ± 0.019), cervello (0.873 ± 0.040); indeterminata nel tumore (0.939 ± 0.043) e nel cuore (0.914 ± 0.075). Regge in 4 tessuti su 6: nel tumore e nel cuore, con 3 semi, P-5 non ha potere per un difetto del 6-9 %.
- **CP-2** (cuore, vetro a pezzi / pezzo solo): confermata a 16-64 µm (1.026 [0.988, 1.061]) e a 64-256 µm (1.069 [1.021, 1.117]); **smentita a 256-512 µm** (2.26 [2.19, 4.0]: 61 UMI contro 27, 108 blocchi in 3 tessere, 100 all'angolo dell'area di cattura); oltre 512 µm non c'è vetro. La smentita non è un errore del confronto: PH-3 mostra che viene dalla riflessione della coda sul bordo dell'area di cattura nei pezzi.

# Regola R-2: lo step successivo

1. P-8a confermata 5/5: il livello cellulare non va ottimizzato prima.
2. Smentite fra le previsioni che contano: P-4a (5 tessuti: il bordo simulato due volte costa troppo) e P-6 (5 giudizi: la seconda coda nei pezzi va trattata sull'area intera); la smentita di P-5 è spiegata dal caso (PH-2). **Lo step successivo è S2.7, simulazione a pezzi corretta dentro il pacchetto**, con: nessun bordo simulato due volte (trascritti in uscita passati al pezzo di arrivo), seconda coda sull'area intera, una sola realizzazione per cellula, strisce dell'area intera (punto 3 di R-2), pezzi multipli di 8 µm. Registrato in [S2.6: lo step successivo è la simulazione a pezzi corretta nel pacchetto (S2.7), con i requisiti dati dalla prova (DEC-073)](../docs/decisions/DEC-073.md).

# Misure (senza intervallo)

- **J-2, cellule simulate due volte**: nella fascia 0-8 µm dalla giuntura da 16 % (linfonodo) a 39 % (cuore) dei trascritti catturati; su tutta l'area dall'1.2 al 2.3 % con pezzi di 1 000 µm, dal 2.7 al 6.1 % con pezzi di 504 µm (`S2.6_q2_J1_J2.csv`).
- **J-3, somiglianza attraverso la giuntura** (correlazione degli UMI fra blocchi adiacenti, media di 3 semi): pezzo solo / pezzi di 1 000 µm con bordo / senza bordo: intestino 0.41 / 0.21 / −0.09; tumore 0.44 / 0.18 / −0.05; stroma 0.57 / 0.23 / −0.21; linfonodo 0.13 / 0.09 / −0.05; cervello 0.75 / 0.48 / −0.24; cuore 0.66 / 0.34 / −0.19 (`S2.6_q2_J3.csv`). Con il bordo la correlazione scende da un terzo a più della metà.
- **Livello cellulare**: sui quadrati 0.8-1.1 s fissi più 0.62-0.74 ms per cellula (retta sui 6 punti per tessuto; nel linfonodo la retta non ha parte fissa, −0.7 s); sulle sezioni intestino 2 172 s (tassellazione 1 876 s, 86 %), linfonodo 1 569 s (1 101 s), tumore 736 s, cervello 62 s, cuore 53 s (`S2.6_q3_step1.csv`).
- **Sovra-generazione delle strisce** (massimo dei fattori / media sul tessuto): quadrati di 9 mm² da 3.4 (cuore) a 5.6 (stroma); sezione intera del linfonodo 7.6. Cresce con la griglia, perché dipende dalla striscia più forte di tutta l'area: sulla sezione del linfonodo 496 milioni di trascritti generati contro 295 previsti da R-1.
- **Sezioni intere** (`S2.6_q3_runs.csv`): in un pezzo solo intestino 6 830 s e 107.4 GB, linfonodo 6 159 s e 87.4 GB, cuore 2 642 s e 26.3 GB (R-1 prevedeva 100, 47 e 49 GB: la regola sbaglia fino a quasi il doppio nei due sensi); a pezzi di 1 000 µm (49 pezzi, nessuno vuoto) da 7 289 s di processo (cuore) a 25 128 s (tumore), 12-86 minuti di orologio, picco del pezzo fino a 20.8 GB (tumore; processo figlio, comprende le pagine del livello dell'area intera).
- **Descrittive contro le sezioni vere** (BL-205): UMI medio per blocco nel tessuto, simulato in un pezzo solo / vero: intestino 1.12, linfonodo 2.32, cuore 0.82; a pezzi / vero: tumore 1.44, cervello 1.41. Sono i rapporti fra le aree di riferimento e la sezione intera misurati in R4 (1.08, 2.39, 0.80, 1.38, 1.26): la sezione simulata con un solo preset ha la profondità delle aree di riferimento. Vetro vicino (16-64 µm), simulato / vero: intestino 1.69, linfonodo 6.9, cuore 0.49 (in un pezzo solo); rapportato al tessuto, il linfonodo simulato ha circa tre volte il vetro vicino vero (BL-244).
