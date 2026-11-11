
Tabelle in `results/S2.7/`: `S2.7_G1_tests.txt`, `S2.7_mutants.txt`, `S2.7_G2_identity.csv`, `S2.7_G4_test_S2.*.txt`, `S2.7_md5_before.txt`, `S2.7_md5_after.txt`, `S2.7_P1.csv`, `S2.7_section_<A>.csv`, `S2.7_P3_P4.csv`, `PH1/S2.7_P1.csv`, `S2.7_PH2_md5.txt`, `S2.7_PH_A1*.csv`. Uscite grezze in `/mnt/micron/geo_spatialtrans/S2.7/` (`G2/`, `p1/`, `sections/`, `p1ph/`, `sections_ph/`). Codice delle esecuzioni: commit `d5e0cfa`.

### Cancelli

- **G-1** superato: batteria 78/78 sul codice finale; mutanti 8/8 (prima stesura 6/8, correzione C-1).
- **G-2** superato: 210 confronti di impronte, 210 identici (sei quadrati di 4 mm², pezzo solo contro pezzi di 504 e 1 008 µm con 4 processi; sezione intera del cuore, pezzo solo contro pezzi di 1 008 µm con 8 e 3 processi e di 504 µm con 8). Il revisore dei numeri ha confrontato byte per byte 75 coppie di file: nessuna differenza.
- **G-3** superato: in 108 esecuzioni della funzione nuova tutti i controlli scritti dalla funzione sono veri (ogni cellula in un solo pezzo, ogni quadratino una volta, conservazione per cellula e totale, UMI = catturati + fondo, nessun NA, somma a 8 µm). Due di questi controlli sono veri per costruzione (vedi la revisione, R-04); quello che rileva i trascritti persi fra pezzi è la conservazione per cellula (mutante M2). Un ripiego della triangolazione in tutto: la cellula 624286 della sezione del linfonodo.
- **G-4** superato: 192 file letti identici prima e dopo; `R/` cambiato solo nei file dichiarati; batterie di S2.1...S2.6 invariate (21, 36, 39, 37, 40, 33 asserzioni, nessun FAIL).
- **G-5** superato: due revisori, nessuna discrepanza numerica, nessun difetto che cambi i risultati (sotto «Revisione avversariale»).

### Previsioni

| id | giudizi | conf. | indet. | smentite | in breve |
|---|---|---|---|---|---|
| P-1 stessa legge del motore di S2.5 (1 mm²) | 24 | 21 | 3 | 0 | cervello UMI e geni 1.016 [1.003, 1.030], 1.016, 1.014: tutte e tre per eccesso; con 27 semi 1.004 [0.995, 1.014], confermate (PH-1) |
| P-2 sezioni a pezzi / pezzo solo di S2.6 | 14 | 13 | 1 | 0 | tessuto 0.999, 1.001, 1.002; vetro oltre 512 µm 1.018 (intestino), 1.013 (linfonodo); cuore 256-512 µm 1.26 [1.15, 4.0] su 108 blocchi all'angolo (34 contro 27 UMI) |
| P-3a picco per pezzo, pezzi di 1 008 µm | 5 | 3 | 0 | 2 | tumore 12.5 GB, cervello 9.6 GB (fase B); con pezzi di 504 µm 3.8 e 3.4 GB (PH-2) |
| P-3b tempo di processo per mm² di tessuto | 5 | 3 | 0 | 2 | intestino 330 s, tumore 388 s; linfonodo 233, cervello 263, cuore 72 |
| P-4 quota troncata dal tetto | 5 | 5 | 0 | 0 | da 4.8e-7 a 5.3e-6; da 142 a 267 quadratini sopra il tetto per sezione |

P-3 è un solo seme: la stima è l'esito, senza errore Monte Carlo (dichiarato).

### Controprova

- **CP-1** superata: UMI per quadratino delle uscite dei pezzi assottigliate a 0.95, smentita per difetto in 6/6 (da 0.947 a 0.965; cervello 0.965 [0.953, 0.978]).

### Misure (senza intervallo)

- Trascritti passati fra pezzi nelle sezioni: da 4.6 milioni (cuore) a 52 milioni (cervello); generati: linfonodo 265 milioni contro 496 della sezione in un pezzo solo di S2.6.
- Tempo per fase (sezioni, pezzi di 1 008 µm, 8 processi): verità geometrica da 1 390 a 2 750 s (dal 17 % del tempo nel tumore al 69 % nel cuore); orologio da 7 a 32 minuti; processo principale fino a 34 GB (tumore, cucitura). Tetto M*: 4.13, 4.70, 4.04, 3.81, 2.98 (intestino, tumore, linfonodo, cervello, cuore).
- Quadrati di 4 mm² in un pezzo solo (G-2): da 249 s (cuore) a 1 382 s (tumore), contro 238-1 190 s del motore di S2.5 in S2.6; a pezzi, tempo di processo da 1.1 a 1.7 volte il pezzo solo (504 e 1 008 µm).
- Non fatte: descrittive di BL-237 (riquadro della finestra in A1_r1) e descrittive delle sezioni contro il vero oltre a quelle nelle tabelle delle sezioni (colonne `ratio_pieces_real`): vedi «Deviazioni».
