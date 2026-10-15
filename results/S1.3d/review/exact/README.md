# S1.3d — revisione avversariale «esatta» (revisore indipendente)

Data: 2026-10-15. Repo `~/2025.geo_spatialtrans`, branch `step1-cell-layer`, commit `1431d75` (sola lettura: nessun file del run modificato).
Script e log qui (`results/S1.3d/review/exact/`); dati pesanti in `/mnt/micron/geo_spatialtrans/S1.3d/review_exact/`.

## Esito in breve
1. **Verità dell'arbitro confermata con un metodo esatto diverso** (integrazione per strisce verticali in `fractions.Fraction`, `rx.py`; non importa
   `exactgeom.py` né `S1.3d_arbiter.py`): 5 358 pezzi O2 (N, H, S), 260 coppie O1, 519 aree di cella: **uguaglianza al double più vicino in 100 %**
   (differenza relativa 0); nessun pezzo mancante o spurio nella verità sulle 42 celle ricalcolate per intero (4 938 pezzi a 2 e 8 µm).
2. **Tutte le 565 violazioni registrate** (coppie e pezzi; V0, F2, F3, F5, anche S) hanno verità esatta corretta e sono **reali** (errore > τ).
   **F1: 0 violazioni** registrate in qualsiasi insieme e operazione (nulla da ricalcolare; confermato anche dai CSV aggregati: CP-3 su H 0/300, 0/600, 0/300).
3. **F5**: polyclip **tronca** (codice `interface.cpp`: `(cInt)((x - x0)/eps)`), la documentazione dice «rounding to the nearest integer».
   Ma **il difetto topologico non è il troncamento**: è Clipper **6.4.0**, una sola riga (`AddLocalMinPoly`) corretta in 6.4.2. Su 3 187 809 pezzi
   (11 casi): 6.4.0 troncando 43 errori > 1e-3 µm², 6.4.0 arrotondando 29, 6.4.0 + la sola guardia del 6.4.2 → 0 (in entrambi i modi), 6.4.2 → 0.
   Il troncamento da solo produce schegge oltre τ (≤ 4e-6 µm²) in U5_cocirc/CP1P con ogni versione di Clipper.
4. **Le 3 cadute di F5 su N (e le 18 su H, le 69 su S) non sono errori di Clipper**: l'uscita intera di Clipper ha area esatta entro τ
   (0.13–0.69 τ); l'errore nasce in `pc_area()` del wrapper (shoelace in virgola mobile su coordinate assolute ~1e3 µm: errori quantizzati a 2^-33…2^-36 µm²).
5. **F4**: 14 celle invalide (13 in N + 1 in H), conteggio confermato da una scansione indipendente di tutte le 18 745 076 celle F4.
   Diagnosi confermata e precisata: l'anello della cella passa **due volte** per un vertice dove **due anelli della regione si toccano**
   (guscio–buco o buco–buco), 3 punti distinti (A5 r1, r2, r4). `st_make_valid` → POLYGON con buco, valido, **stessi vertici, area esatta identica** (Δ = 0).
6. **V0 in O1**: in tutte le 15 unioni sbagliate di N l'ingresso è valido in esatto, la verità è chiusa fra limiti esatti (larghezza ≤ 3.9e-13 µm²)
   e il deficit di `st_union` è **l'area esatta di una sola cella** (scarto ≤ 2.4e-9 µm²): errore dell'operazione. Le 20 cadute di V0 con
   copertura valida sono reali (22 coppie: esatto = 0, V0 = 17–1 933 µm², 3.1e5–5.1e6 τ). Su 167 coppie violate, 162 restituiscono l'area **intera** di una delle due celle.
7. **D-3** confermato: in R `paste()` di un double intero differisce da `sprintf("%d")` solo per multipli di 1e5 (11 valori in [0, 2e6]);
   fuori da S max n = 31 070, max indice di bin 508 → nessun caso toccato; in S 70 + 14 pezzi esatti delle celle 1e5…4e5 → 140 + 28 violazioni fittizie.
**Nessun rilievo cambia la decisione (F1 per O1, O2, O3)**; due cambiano la *lettura* di F5 (punti 3–4).

## Metodo
- `rx.py`: aree esatte per **strisce verticali**: fra ascisse di rottura consecutive (vertici, incroci fra lati, incroci con y0/y1 del bin) la
  lunghezza della sezione verticale è lineare in x, quindi l'integrale = larghezza × lunghezza nel punto medio, esattamente. Sezione = regola
  pari-dispari. Coordinate = `Fraction(double)` dai parquet di `geom/`. Confronto: `float(Fraction)` (arrotondamento corretto) == valore della verità.
  `test_rx.py`: 17/17 PASS (casi analitici, buchi, incavo, simmetria, accordo con shapely ≤ 1e-13).
- Validità esatta degli anelli (`v0_union.py`): test di orientazione in Fraction su tutte le coppie di lati non adiacenti.
- Clipper: sorgenti 6.4.0 (da `polyclip_1.10-7.tar.gz`, CRAN) e 6.4.2 (da `pyclipper-1.4.0`), compilati senza R (`clipdrv.cpp`, g++ 13.3),
  più la variante 6.4.0 con la **sola** riga 1930 del 6.4.2; pyclipper 1.4.0 e pyclipr 0.1.8 (Clipper2 2.0.1) in un venv proprio
  (`/mnt/micron/.../review_exact/venv_clip`). Area dell'uscita intera calcolata in esatto (`__int128`).
- Campioni a seed fisso (`select_sample.py`, seed 20261015; `select_sample_extra.py`, seed 20261016).

## Tabella dei risultati
| # | Quantità | Riportato | Ricalcolato (esatto, indipendente) | Esito |
|---|---|---|---|---|
| 1a | 300 pezzi O2 2 µm, 10 nulli N (A1×2, A2×2, A3×2, A4×2, A5, A6) | truth/*_o2_s2 | 300/300 uguali al double | ✓ |
| 1b | pezzi di 30 celle intere (N), 2 e 8 µm | truth | 2 792/2 792 uguali; 0 mancanti, 0 spuri; Σ pezzi = area cella (Fraction) 60/60 | ✓ |
| 1c | H (3 nulli) e S: pezzi + celle intere | truth | 2 161 + 105 uguali | ✓ |
| 1d | 200 coppie O1 casuali (20 nulli N, 3–4 per archetipo) + 10 di prova + 40 H + 10 S | truth/*_o1_pairs | 260/260 uguali, flag `zero` 260/260, simmetria esatta | ✓ |
| 1e | aree di cella O1 | truth/*_o1_cells | 519/519 uguali | ✓ |
| 1f | violazioni F1 | 0 | 0 righe in `res/*.rds` (tutti gli insiemi) | ✓ |
| 1g | violazioni registrate (coppie V0, pezzi F2/F3/F5, S) | 565 | verità 565/565 uguale; 565/565 errori reali | ✓ |
| 2a | caso CP1T_A1_r3_CSR_16_d08, cella 2206, bin (151,109) | 0.9130726 vs 1.653643 | esatto 0.9130726; polyclip 1.653643 riprodotto bit per bit; 6.4.0 standalone 1.653643; 6.4.2 / Clipper2 ≈ 0.91307 | ✓ |
| 2a | F5 errori > 1e-3 µm² fra le violazioni | 145 (CP1T) | 145/145 riprodotti dal 6.4.0 standalone (|Δ| < 1e-6); 0/145 con 6.4.0+guardia, 6.4.2, pyclipper, Clipper2 | ✓ |
| 2a | «polyclip tronca» | ipotesi | sì: sorgente + sonda (1.7e-6 → 1e-6) + stessi risultati con interi pre-troncati | ✓ |
| 2a | il troncamento causa il difetto? | ipotesi | **no**: con arrotondamento 6.4.0 dà ancora 29 errori grandi su 3.19 M pezzi; con la guardia 0 | ✗ (ipotesi respinta) |
| 2b | celle F4 invalide | 13 (N) | 13 in N + 1 in H (A5_r4_RSArule_63); scansione completa: 14/18 745 076; S: 0 | ✓ (+1 H) |
| 2b | st_make_valid su F4 | da verificare | 14/14 valide, POLYGON con buco, Δarea esatta = 0 (τ_i 2e-4–9e-4) | ✓ |
| 2c | V0: unioni sbagliate in N | 15 | 15/15: deficit = area esatta di una cella (≤ 2.4e-9 µm²); ingressi validi; coppie complete | ✓ |
| 2c | V0 N copertura valida ma silenziosa | 20 | 20 reali: 4 con unione (3 + coppia, 1 sola unione), 16 solo coppie; 22 coppie, esatto 0, V0 17–1 933 µm² (3.1e5–5.1e6 τ), 20/22 = cella intera | ✓ |
| 2d | F2/F3 H A4_p2_CSR_62 8601 (380,469) | 1.42e-17 vs 2.08e-14 | esatto 1.4190e-17; estensione 0.0039 g × 0.0072 g; cella arrotondata a g → 2.0810829e-14 (= candidato a 5e-7 rel.); 1.07 τ | ✓ scheggia sub-g |
| 2d | 3 F5 su N | violazioni | pezzi 3.5–14 g, larghezza media 1.1–3.1 g; uscita intera di Clipper entro τ (0.13–0.40); errore = `pc_area` in virgola mobile | rilettura |
| 2e | D-3: valori toccati da paste() | multipli di 1e5 | 11 valori in [0, 2e6], tutti multipli di 1e5; fuori da S nessuno; S 140 + 28 | ✓ |

## Rilievi
- **R-1 (confermato)** Verità dell'arbitro corretta sui campioni e su tutte le violazioni (tabella righe 1a–1g).
- **R-2 (discrepanza di lettura, F5)** Il difetto topologico di F5 è di **Clipper 6.4.0** (riga 1930, `AddLocalMinPoly`: manca la guardia
  `prevE->Top.Y < Pt.Y && e->Top.Y < Pt.Y` del 6.4.2), non del port R né del troncamento. Le dev_notes lo lasciavano indistinto. File: `f5/f5_versions.csv`, `f5/broad_summary.csv`, `src/clipper_640_vs_642.diff`.
- **R-3 (bug documentale di polyclip)** polyclip tronca verso zero invece di arrotondare come documentato: spostamento fino a √2·eps invece di eps/√2;
  rapporto errore/τ massimo con Clipper corretto 0.684 (troncando) contro 0.502 (arrotondando) su 2.68 M pezzi; 26 pezzi oltre τ (U5_cocirc 25, CP1P 1) con ogni versione.
- **R-4 (metodo, wrapper F5)** `pc_area()` calcola lo shoelace su coordinate assolute: errore di cancellazione ~2^-34 µm² (5.8e-11) che supera τ per pezzi
  minuscoli. Tutte le violazioni F5 di N (3), H (18 righe) e S (69) sono di questo tipo (uscita intera entro τ). Correzione: sottrarre un'origine locale o usare l'area intera.
  Non cambia la decisione (F5 resta non eleggibile in O2 per R-2 su CP1T, e non è ≥ 2× più veloce in O1).
- **R-5 (confermato, F4)** Le 14 celle invalide: un vertice ripetuto due volte nell'anello, che è il punto di contatto fra due anelli della regione
  (A5 r1: buchi 182/183; r2: guscio/buco 2; r4 parte 2: guscio/buco 2) — non «la regione che tocca sé stessa» con un solo anello. `st_make_valid` le ripara senza cambiare area né vertici.
- **R-6 (confermato, V0)** Gli errori di V0 sono dell'operazione: unione = cella intera persa (cfr. GEOS #965); intersezione di coppia = area intera di una cella (162/167) o parte grande (5/167).
- **R-7 (metodo, tolleranza dei pezzi sub-g)** Per un pezzo esatto sub-g τ = g·P_esatto (1.9e-14) è più stretto dello spostamento ammesso dalla griglia;
  il candidato F2/F3 è il pezzo esatto della cella arrotondata (sub-g, perimetro 0.74 g). Scheggia sotto la risoluzione: la violazione è un artefatto della definizione di τ, non un errore di F2/F3. Non tocca F1.
- **R-8 (confermato, D-3)** Solo multipli di 1e5; i percorsi `paste()` rimasti (`o3_aggregate`, `tau3`) sono innocui: chiavi intere lato verità e `as.integer("1e+05") = 100000`.

## Limiti dichiarati
- Le unioni/differenze simmetriche sbagliate ricalcolate in esatto: le 15 di N (non le 6 di U2, 4 di H, 1 di CP1P). Per le altre violazioni di unione ho verificato solo che la verità di coppie e celle è corretta.
- La controprova estesa di Clipper usa solo i pezzi della verità (non prova i pezzi spuri) su 11 casi.
- I valori dei candidati non violanti non sono salvati nei `res` (solo i riassunti): F1 è verificato tramite assenza di violazioni e `max_ratio`, non pezzo per pezzo.

## File
`rx.py`, `test_rx.py`(+log), `select_sample*.py`, `review_sample.py`(+log), `dump_viol.R`(+log), `review_viol.py`(+log), `f5_cases.R`(+log), `clipdrv.cpp`,
`f5_clipper_versions.py`(+log), `f5_broad.py`(+log), `sliver.py`(+log), `f4_check.py`, `f4_makevalid.R`, `f4_makevalid_exact.py`, `f4_scan.py`(+log),
`v0_union.py`(+log), `v0_pairs.py`(+log), `d3_paste.R`(+log). Dati: `review_exact/{sample_all.csv, viol_all.parquet, viol_recheck_all.csv, f5/, f4/, v0/, v0_union_summary.csv, v0_pairs.csv, sliver.csv, build/, src/}`.
