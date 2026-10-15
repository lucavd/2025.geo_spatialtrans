# S1.4b — revisione avversariale del CODICE (revisore: code)

Commit rivisto: `1c13c4c` (branch `step1-cell-layer`). Nessun file del repo o dei grezzi modificato; tutto il codice e i
risultati di questa revisione sono in questa cartella, intermedi in `/tmp/code_adv`, `/tmp/code_vc`, `/tmp/code_tgt`.
File rivisti: `tools/S1.4b_laguerre.R`, `tools/S1.4b_bench.R`, `tools/S1.4b_arbiter.py`, `tools/S1.4b_cp1.py`,
`tools/S1.4b_analysis.R`, `R/testing/test_S1.4b.R`, `tools/S1.4b_mutants.R`, contro pre-registrazione e dev notes.

## Sintesi
- Il diagramma di potenza è **corretto**. Non ho trovato controesempi per `lag_neighbors()` nei casi non degeneri.
  Un arbitro esatto (certificato) sui 30 insiemi reali trova 0 celle fuori tolleranza. Il FAIL di C-L2 è un
  **artefatto della statistica** di C-L2, non del codice (CODE-05).
- Tutti i numeri B, CP-1 e il punteggio 4 vs 4 si **riproducono a tolleranza 0** (sezione «Ricalcolo»).
- La decisione «(c)» poggia su un **pareggio** risolto da una regola aggiunta nelle dev notes (CODE-08). Si **ribalta** se
  per V-c si usa l'operatore di osservazione della pre-registrazione, cioè il Voronoi dei centri dei nuclei simulati
  (CODE-07): V-c 4 → 3 PASS e {a,b} vince.
- Il troncamento di V-a/V-b dipende dall'algoritmo: la passata unica tronca **tutti** i nuclei in conflitto, cioè circa 2×
  il minimo teorico. B-4 resta FAIL in A1–A5 per il vincolo stesso; in A6 diventa PASS con un algoritmo iterativo. Il
  conteggio PASS (4) non cambia (CODE-06).

## Rilievi

| id | gravità | cambia un verdetto? |
|---|---|---|
| CODE-01 | bug (minore) | no |
| CODE-02 | osservazione | no |
| CODE-03 | osservazione | no |
| CODE-04 | osservazione (verifica) | no |
| CODE-05 | difetto di metodo | **sì** (C-L2 e l'eleggibilità di V-a; la decisione no, a operatore invariato) |
| CODE-06 | difetto di metodo | **sì** (B-4 A6 di V-a/V-b dipende dall'algoritmo; punteggio e decisione no) |
| CODE-07 | difetto di metodo | **sì** (punteggio V-c 4 → 3; decisione (c) → (b), oppure (a) con CODE-05) |
| CODE-08 | difetto di metodo | no da solo; rende fragile la decisione |
| CODE-09 | osservazione | no |
| CODE-10 | osservazione | no |
| CODE-11 | osservazione | no |
| CODE-12 | osservazione | no |

**CODE-01 — `lag_cells()` con n = 1 restituisce la cella in coordinate centrate.** Quando un punto non nascosto non ha
vicini (solo n = 1), `cells[[i]] <- P` salta la traslazione `+ c(mx, my)`. Caso C8d (un punto in (3,3), regione
10×10): area 64 invece di 100, errore relativo 0.36 (`adv_R_summary.csv`). Nel banco non ha effetto (n ≥ 524).

**CODE-02 — qhull fallisce sui casi degeneri.** Con tutti i punti collineari (C3a, n = 20) e con i 5 punti sollevati
coplanari (C8c: quadrato più centro a pesi equilibrati) `convhulln(..., "Qt n")` termina con QH6154 «initial simplex is
flat». Il fallimento è esplicito (nessun risultato sbagliato) e non c'è ripiego. Nei dati del banco non si presenta.
Proposta per S1.4c: ripiego su tutte le coppie per n piccolo, oppure `QJ`.

**CODE-03 — generatori quasi coincidenti a pari peso.** qhull fonde i punti a distanza ≤ 1e-10 µm. Uno risulta nascosto
e l'altro prende l'unione delle due celle: errore d'area 41.7 µm² (offset 1e-10) e 24.5 µm² (offset 1e-12). Fino a
1e-8 µm le celle coincidono con l'arbitro. Con pesi diversi (+0.5) il nascosto è sempre corretto. La distanza minima
fra generatori reali nei 30 insiemi è **1.40 µm** (`cl2_zstat_exact.csv`), quindi il caso è irrilevante.

**CODE-04 — casi avversari: nessun controesempio per `power_diagram()` / `lag_neighbors()`.** Ho usato 17 configurazioni
non degeneri × 2 modi (`adv_cases.R`): pesi eterogenei fino a 1e4 µm², pesi negativi, fino a 170/200 nascosti,
grappoli di nascosti multipli, reticoli cocircolari, collineari parziali, regione multipla con buchi e generatori fuori,
coordinate traslate di 1e4 e 1e5 µm, celle che non contengono il generatore, RSA con dischi disgiunti (teorema).
Gli arbitri sono miei (`adv_arbiter.py`):
(i) intersezione di **tutti** gli n−1 semipiani per ogni cella;
(ii) rasterizzazione a passo 0.1 µm con argmin della potenza per forza bruta.
Risultati (`adv_arbiter.csv`):
- modo esatto: errore relativo d'area ≤ 4.7e-11 (1e4 µm) e ≤ 3.5e-12 (1e5 µm), soglia 1e-9;
- modo su griglia: |ΔA| ≤ 5.5e-5 µm², cioè entro perimetro × 1e-6;
- celle vuote identiche;
- **0 pixel assegnati male, scoperti o doppi** (oltre 1e-4 µm dal bordo);
- teorema: 0 nuclei fuori dalla cella su 400.

**CODE-05 — C-L2: la soglia sulle aree assume la normalità. Il FAIL 22/30 è un artefatto.** I conteggi per cella sono
binomiali con media ~2–60 (coda superiore pesante). P(|z| > 4) esatta, sommata sulle celle dei 30 insiemi, dà
**E = 72.5**; l'approssimazione normale del codice dà **15.28**; le osservate sono **68** (`cl2_zstat.py`,
`cl2_zstat_exact.csv`).
Controprova esatta, non campionaria (`exact_cells_arbiter.py`, `exact_cells_arbiter.csv`): ogni cella è ricostruita dai
semipiani dei 40 vicini più tutti i violatori trovati, con un certificato ai vertici contro tutti i generatori. Nei 30
insiemi:
- max |ΔA|/perimetro = 4.5e-7;
- 0 celle con differenza simmetrica > perimetro × 1e-6;
- 0/0 celle vuote;
- Σ aree relativa ≤ 3.6e-16.

Il caso più estremo, A1 fold 3 (8 celle |z| > 4, p < 0.002 in PH-B), dipende dal seme dei punti: con altri 3 semi si
ottengono 3, 1 e 3 celle. Il conteggio per potenza minima coincide cella per cella con il conteggio nei poligoni di R
(`inspect_cl2.py`). Il codice è quindi corretto e V-a sarebbe eleggibile. Da sola, però, l'eleggibilità di V-a **non
cambia** la decisione: V-a e V-b hanno metriche B identiche, quindi resta 4 vs 4 e il pareggio va a (c).

**CODE-06 — `nonoverlap_scale()` (una passata) è corretta ma conservativa; B-4 di A6 dipende dall'algoritmo.**
- Correttezza: s_i ≤ f_ij per entrambi i membri della coppia, quindi il vincolo è soddisfatto. `near_pairs()` è
  completo: un ricalcolo con tutte le coppie da KD-tree coincide con `s` del banco entro 1e-12 in 10/10 insiemi
  (`nonoverlap_opt.csv`, colonna `onepass_equals_bench`).
- Conservatività: la frazione troncata è **esattamente** la frazione di nuclei con almeno un conflitto (10/10).
- Minimo esatto (MILP di vertex cover con vincoli forzati), A2/A4, s = 1: **A2 0.229–0.285 (media 0.261), A4
  0.227–0.242 (media 0.233)**, contro una passata 0.559 e 0.523.

Sensibilità su tutti i 30 fold × 10 semi con le stesse regole (`nonoverlap_sensitivity.py`,
`nonoverlap_sens_archetypes.csv`). La passata unica riproduce esattamente i B-4 del CSV. Frazione troncata B-4:

| algoritmo | A1 | A2 | A3 | A4 | A5 | A6 |
|---|---|---|---|---|---|---|
| una passata | 0.380 | 0.557 | 0.357 | 0.525 | 0.185 | 0.0875 |
| greedy («riduci il più grande») | 0.193 | 0.287 | 0.184 | 0.271 | 0.120 | **0.0488 (PASS)** |
| LP max Σ R·s | 0.205 | 0.320 | 0.203 | 0.280 | 0.131 | 0.0518 |

- B-4 FAIL in A1–A5 dipende dal vincolo, non dall'algoritmo: il minimo teorico in A2/A4 è > 0.05.
- In A6 B-4 dipende dall'algoritmo (0.0875 → 0.0488).
- P2 resta FAIL: in A2 il minimo è 0.261 > 0.5 × 0.327.
- Il punteggio {B-1, B-2, B-3, B-7} di V-a/V-b vale **4 con tutti e tre gli algoritmi**. Cambiano alcune celle senza
  PASS guadagnati o persi: A2 B-3 FAIL → WARN; A4 B-1 e B-3 FAIL → WARN (LP); A3 B-2 FAIL → WARN (LP).

**CODE-07 — l'operatore di B-7 (e di B-1/B-3) per V-c non è quello dei reali; con l'operatore della pre-registrazione la
decisione si ribalta.**
- Reali: cerchio di pari area **centrato nel centroide del nucleo**, confrontato con il **Voronoi dei centroidi dei
  nuclei**. Il generatore è il centroide del nucleo (`S1.4_bench.R`, `d_gen`).
- V-b/V-a: centro = generatore, quindi stesso operatore.
- V-c (`S1.4b_bench.R` r. 96): cerchio centrato nel **POI**, confrontato con il Voronoi dei **generatori**
  (`center_dist`). Quel Voronoi non è osservabile dai nuclei simulati. Il POI dista dal generatore (mediana) 1.7, 1.8,
  2.7, 0.6, 5.8 e **7.0 µm** in A1–A6.
- Con il clamp (k32·r·s ≤ 0.85·d) si ha r < d, quindi `b7` ≡ 0 per costruzione (P4 lo dichiara). La definizione è nelle
  dev notes, quindi non è un bug rispetto a esse. Contraddice però l'«operatore di osservazione comune» della
  pre-registrazione («rispetto al Voronoi dei centroidi»).
- Ricalcolo con il Voronoi dei POI ritagliato alla stessa finestra (`vc_poi.R`, `verdicts_recompute.R`,
  `recomputed_archetypes.csv`):
  - b7_sim di V-c = 0.000–0.016 (B-7 resta FAIL ovunque);
  - A6 B-1 WARN → FAIL;
  - A6 B-3 PASS → WARN (fold 1: d_sp 0.131 → 0.156; mediana del metro 0.151);
  - punteggio V-c **4 → 3**.
- Con la regola implementata, passo 2: {a,b} 4 > 3. Passo 3: V-a non è eleggibile (C-L2), quindi la decisione diventa
  **(b)**. Con CODE-05 (C-L2 ricondotto all'artefatto) e CP-1 che non smentisce (a), diventa **(a)**.
- Lo stesso operatore era già in S1.4 (C-0 lo eredita).

**CODE-08 — la decisione poggia su un pareggio e su una regola non nella pre-registrazione approvata.** La regola
«a parità resta (c)» è stata aggiunta nelle dev notes con «da confermare con Luca nel report». È implementata come
scritta (`if (s_ab > s_c) "ab" else "c"`). Il margine è 0 PASS e un solo fold (CODE-07) lo ribalta. Il punto va portato
a Luca prima di dichiarare la decisione.

**CODE-09 — C-0 verifica davvero l'identità, ma non copre B-7.** `identical()` bit a bit su (ks_nc, ks_an, d_sp,
clamp) × 330 corse. `o$seeds` è una matrice numerica 10×4, quindi il confronto riga = seme è significativo e
fallirebbe se le righe fossero sfasate. Il mio ricalcolo coincide con i grezzi a tolleranza 0 in 30/30 fold per V-c e
per V-a/V-b. C-0 non riguarda la nuova metrica B-7 di V-c, che è identicamente 0 (CODE-07).

**CODE-10 — CP-1: implementazione corretta; il placebo è debole; il controllo R3 è una riproduzione.**
- `lab2gen` ricostruito con un dizionario indipendente: stessi pixel. `interior` allineato (assert su `idx`).
  Permutazioni riprodotte dallo stesso flusso RNG. Sollevamento √(W − w) → |p − c|² + W − w: segno corretto.
- Ricalcolo senza sollevamento (`cp1_indep.py`: k vicini euclidei con certificato d_K² ≥ pow_own + W, forza bruta sui
  pixel non certificati), 660 condizioni:
  - max |Δm1| = 2.7e-5 e max |Δm2| = 1.0e-4, dovuti ai pareggi esatti sulla griglia dei pixel, più frequenti con
    λ = 0; il KD-tree li conta «fuori» in modo arbitrario;
  - «batte» 30/30 ROI, cioè 6/6 archetipi, come nel CSV;
  - medie d'archetipo identiche alla tabella (`cp1_archetype.py`).
- Il controllo «λ = 0 = R3» usa lo stesso `lab2gen`, la stessa convenzione pixel (c·upp, r·upp) e lo stesso KD-tree di
  `R3_containment.py`: verifica la riproduzione, non è un controllo indipendente.
- Placebo debole: **anche il Voronoi (pesi nulli) batte tutte e 20 le permutazioni in 29/30 ROI** (eccezione A6 r2).
  Pesi permutati aggiungono rumore, quindi «battere le permutazioni» è quasi garantito. Il contrasto informativo è
  Laguerre vs Voronoi: 29/30 ROI; in A6 r5 m2 peggiora del 2.5 %.

**CODE-11 — analisi: regola e parità implementate come scritte; due dettagli.**
(i) Ramo irraggiungibile ma sbagliato: con `elig_a & refuted & !elig_b` il codice sceglie «(a) Laguerre» benché (a) sia
smentita.
(ii) P1 usa il rapporto delle medie (1 − mean m2_lag / mean m2_vor). Le dev notes («media sui ROI») sono ambigue. La
media dei rapporti cambia solo A5 (0.656 → 0.755); P1 resta FAIL (A1–A4 ≤ 0.23).
Ricalcolati e identici anche: regola d'archetipo (≥ 4 P → PASS; ≥ 2 F → FAIL), mediana sui 10 semi, metro, B-4 su
tutte le celle, eleggibilità (a: C-L1, C-L1 mut, C-L2, C-L3, C3, C-NO, C-mut; b: C-NO; c: C-0).

**CODE-12 — test e mutanti.**
- `test_S1.4b.R` T6 calcola `pdh` senza usarlo.
- Mancano test per n = 1 (CODE-01), collineari (CODE-02) e quasi coincidenti (CODE-03).
- Nessun mutante colpisce la selezione dell'inviluppo inferiore (`normals[,3] < 0`) né la gestione dei nascosti; il
  difetto sarebbe comunque visibile a C-L2 e C-L3.
- I mutanti non verificano la calibrazione della statistica di C-L2 (CODE-05).

## Ricalcolo (confronto con `S1.4b_verdicts.csv`)

| voce | esito |
|---|---|
| V-c e V-a/V-b, B-1/B-2/B-3/B-7 per fold e per archetipo (stringhe fold incluse); B-4 0.3803/0.5572/0.3569/0.5247/0.1847/0.0875 (V-a/b) e 0.1636/0.3267/0.1786/0.4782/0.0679/0.0030 (V-c) | **identici** |
| punteggio PASS {B-1, B-2, B-3, B-7} | **4 vs 4, identico** |
| metriche per seme, mio ricalcolo vs grezzi `loo` | uguali a tolleranza 0 in 30/30 fold |
| CP-1: m1/m2 vor, lag, perm per archetipo; 5/5 ROI in 6/6; P1 | **identici** |
| C-L2: 22/30 e 68 celle con \|z\| > 4 | riprodotti; E esatta 72.5 (non 15.28); arbitro esatto 30/30 corretti |
| decisione | **(c) riprodotta**; diventa (b) o (a) con l'operatore di CODE-07 |

## File di questa revisione
`adv_cases.R`, `adv_arbiter.py`, `adv_arbiter.csv`, `adv_R_summary.csv`, `inspect_c2.py` (CODE-01..04) ·
`cl2_zstat.py`, `cl2_zstat_exact.csv`, `exact_cells_arbiter.py`, `exact_cells_arbiter.csv`, `inspect_cl2.py` (CODE-05) ·
`nonoverlap_opt.py`, `nonoverlap_opt.csv`, `nonoverlap_classify.py`, `nonoverlap_alt_verdicts_s1.csv`, `export_targets.R`,
`nonoverlap_sensitivity.py`, `nonoverlap_sens_aggregate.py`, `nonoverlap_sens_raw.csv`, `nonoverlap_sens_folds.csv`,
`nonoverlap_sens_archetypes.csv` (CODE-06) ·
`vc_poi.R`, `verdicts_recompute.R`, `recomputed_folds.csv`, `recomputed_archetypes.csv`, `inspect2.R` (CODE-07..09, ricalcolo) ·
`cp1_indep.py`, `cp1_indep_raw.csv`, `cp1_indep_roi.csv`, `cp1_compare.py`, `cp1_archetype.py` (CODE-10, CODE-11).

## Limiti
- L'ottimo esatto di non sovrapposizione (MILP) è calcolato solo per A2/A4 con s = 1, come chiesto. La sensibilità
  dell'algoritmo (greedy, LP) copre invece 30 fold × 10 semi.
- L'insieme di valutazione del Voronoi dei POI (CODE-07) usa `interior` del Voronoi dei generatori.
- I casi avversari sono sintetici, con n ≤ 400.

## Nota di verifica (CODE-01)
`grep C8d_n1 adv_R_summary.csv` →
```
/home/user/2025.geo_spatialtrans/results/S1.4b/review/code/adv_R_summary.csv:52:"C8d_n1","snap",1,"",0,0,0.36,NA
/home/user/2025.geo_spatialtrans/results/S1.4b/review/code/adv_R_summary.csv:53:"C8d_n1","exact",1,"",0,0,0.36,NA
/tmp/code_adv/R_summary.csv:52:"C8d_n1","snap",1,"",0,0,0.36,NA
/tmp/code_adv/R_summary.csv:53:"C8d_n1","exact",1,"",0,0,0.36,NA
53:add("C8d_n1", 3, 3, 1, st_polygon(list(sq(0, 0, 10, 10))))
```
