# S1.3: revisione del codice (RA-code)

**Revisore**: sotto-agente "revisore del codice", 2026-10-09. **Oggetto**: `R/04b3_tessellate_voronoi.R` al commit `37e68ad` (branch `step1-cell-layer`).
**Vincoli rispettati**: si è scritto solo in `results/S1.3/review/code/` e in `/mnt/micron/geo_spatialtrans/S1.3/{review_code_data,review_venv}`. Pacchetto, risultati, report e git non sono stati toccati.
Ogni numero di questo file viene da uno dei file elencati sotto. Esiti dettagliati: `findings.csv`, 15 rilievi (2 bug, 5 discrepanze, 3 di metodo, 5 conferme).

## Metodo
1. **Lettura del codice** contro la pre-registrazione (D-S1.3.1–3, deviazioni 1–3, definizioni operative), l'addendum C-10 e il design §5.3. Casi limite eseguiti in `edge_cases.R`: 3 punti collineari, 4 cocircolari, generatore sul bordo e su un vertice, isola MULTIPOLYGON senza generatore, due quadrati che si toccano in un punto, smussatura di territori multiparte, `st_make_valid` nella smussatura. Inoltre `.tv_fragments` è stata tracciata su tutti i 40 ROI reali (lunghezza condivisa, unione multiparte, esito dell'aggancio), in una copia isolata della funzione.
2. **Reimplementazione indipendente in Python** (`reimpl_voronoi.py`; venv con numpy 2.5.3, scipy 1.18.1, shapely 2.2.0 / GEOS 3.14.1):
   - tile con **Qhull** (`scipy.spatial.Voronoi`), regione per regione, più 8 generatori fittizi a distanza > diagonale della regione: tutte le celle sono limitate e, dentro la regione, uguali alle celle vere;
   - area di tile ∩ regione **senza GEOS**: Sutherland–Hodgman di ogni anello della regione (buchi compresi, con segno) contro i semipiani delle bisettrici, in coordinate locali al generatore;
   - pezzi del ritaglio con shapely; la regola **D-S1.3.2 applicata da noi**, con lunghezze condivise stimate campionando i bordi dell'orfano (passo 0.002 µm, spostamento di 1e-5 µm verso l'esterno). Due varianti: "simultanea" (contro i soli pezzi principali, catene nelle passate successive; `reimpl_voronoi.py`) e "sequenziale" nell'ordine del codice R (`check_fragments_seq.py`).
   - Ingressi esportati da R (`export_inputs.R`), con coordinate a 17 cifre e geometrie in WKB esadecimale (esatto):
     (a) I1–I5 seed 1, con `load_case()` + `seed_case()` presi dal testo di `R/testing/test_S1.3.R`;
     (b) A1 r1, A2 r1, A3 r1, A5 r1 (cellpose_rgb) e A4 f1, A6 r1 (spaceranger): centroidi da `R3/real/*_gen.parquet`, regione `one_region(r3_roi_window())`;
     (c) A3 r3 CSR replica 3, punti da `gen_points()`.
3. **Mutanti e check**: lettura di `tools/S1.3_mutants.R` e riconteggio dei FAIL per check in `S1.3_test_geos_{M1..M6,none}.csv`.
4. **Arbitro** (`arbiter_check.R`, `check_arbiter_c2.py`, `check_fragments_seq.py`): `hp_cell` su esagono, quadrato, quadrato perturbato di 1e-9, celle sul bordo di un buco, un controesempio costruito alla regola di arresto, e tutte le celle con tile nel quadrato di A3 r1 e A5 r1. Le **26/26 righe** di `S1.3_c10a_arbitration.csv` sono ricalcolate con Qhull (punti nulli rigenerati con `gen_points`, replica 1).
5. **C-2, A3 r3 CSR rep 3** (`check_arbiter_c2.py`, `exact_overlap_766_1987.py`): predicati e intersezione con shapely/GEOS 3.14.1; intersezione 766 ∩ 1987 in **aritmetica razionale**; copertura con 2e5 punti.

## Comandi per rifare (sul server)
```bash
cd ~/2025.geo_spatialtrans
export R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent
Rscript --vanilla results/S1.3/review/code/export_inputs.R        # ~1 min -> review_code_data/*_{cen,reg,tv,info}.csv
Rscript --vanilla results/S1.3/review/code/arbiter_check.R        # ~1 min -> arbiter_known_cases.csv, review_code_data/arbiter_systematic.csv, arb_null_*.csv
PY=/mnt/micron/geo_spatialtrans/S1.3/review_venv/bin/python        # creato con: ~/.local/bin/uv venv $V --python 3.12; uv pip install numpy scipy shapely pyarrow pandas pillow
$PY results/S1.3/review/code/reimpl_voronoi.py I1 I2 I3 I4 I5 A1_r1 A2_r1 A3_r1 A4_f1 A5_r1 A6_r1 A3_r3_CSR03   # ~20 s
$PY results/S1.3/review/code/check_arbiter_c2.py                   # arbiter_recompute_qhull.csv + caso 766/1987
$PY results/S1.3/review/code/check_fragments_seq.py A1_r1 A2_r1 A4_f1 A5_r1 A6_r1 I1 I2 I3 I4 I5 A3_r3_CSR03
Rscript --vanilla results/S1.3/review/code/edge_cases.R           # ~5 min (30 core) -> edge_cases.csv, fragments_trace_40roi.csv
$PY results/S1.3/review/code/exact_overlap_766_1987.py; $PY results/S1.3/review/code/zero_edges.py
```
`reimpl_summary.csv`, `reimpl_orphans_simultaneous.csv` e `reimpl_regions.csv` concatenano i file per caso `review_code_data/review_<caso>_*.csv`. I log sono in `logs/`.

## Esiti
**Compito 2 (reimplementazione), (i) celle non coinvolte nei frammenti.** 62 121 celle in 12 casi. Errore relativo massimo dell'area R rispetto all'area senza GEOS: **5.3e-13** (A3 r3 CSR03); Σ per regione: ≤ 3.8e-15. → RA-code-01, confermato.
**(ii) Celle coinvolte.** Nell'ordine del codice (emulazione sequenziale) 620/622 orfani hanno lo stesso destinatario e le aree coincidono a 5.6e-14 in 10/11 casi. Le 2 eccezioni sono contatti puntiformi in A6 r1 (max rel 2.4e-3). Con la variante simultanea cambiano 12/622 destinatari, con aree fino al 68 % diverse: la regola dipende dall'ordine (RA-code-03). Contatto puntiforme e auto-riassegnazione al donatore: **bug** RA-code-04. Agganci inefficaci (39/40): RA-code-05.
**Compito 1 (codice vs pre-registrazione).** Casi limite corretti (RA-code-14). Divergenze: pezzi interi scartati dalla smussatura (RA-code-06), `st_make_valid` non contato (RA-code-07), confini non annodati (RA-code-11), testo del design §5.3 non aggiornato (RA-code-15).
**Compito 3 (mutanti).** M1–M6 implementano il difetto descritto e i conteggi tornano (RA-code-13). C-8 passa per costruzione; nessun test sul destinatario degli orfani; M4 rilevato solo da C-7; M1 non verificato con C-5 (RA-code-12).
**Compito 4 (arbitro).** La regola di arresto di `hp_cell` può sbagliare: controesempio +3.9 %, bordo di un buco fino a +11.6 %, 2 celle reali in A5 r1 (**bug** dello strumento, RA-code-08). Sulle 26 righe di C-10a l'errore non compare: aree GEOS = Qhull in 26/26, deldir sbaglia > 1e-9 in 16/26. Sui lati invece Qhull dà ragione a deldir in 10/26 righe (vertici degeneri in A4), quindi la frase del report va corretta (RA-code-09).
**Compito 5 (C-2).** Intersezione esatta 766 ∩ 1987 = **1.19e-14 µm²**, mentre GEOS 3.12.1 e 3.14.1 riportano 340.427 µm²; 0/5000 punti della 766 cadono nella 1987. La diagnosi del report è confermata (RA-code-10). La causa è un vertice condiviso che differisce all'ultima cifra (RA-code-11).

## Limiti dichiarati
- La reimplementazione copre la tassellazione a cs = 0. La smussatura è verificata solo con controlli mirati (RA-code-06/07), non reimplementata.
- I ritagli per ottenere i *pezzi* usano shapely, cioè GEOS 3.14.1, una versione diversa da quella di sf (3.12.1). Le aree di controllo (Sutherland–Hodgman) non dipendono da GEOS.
- Le lunghezze condivise sono stimate per campionamento (passo 0.002 µm): 9 orfani hanno scarto fra primo e secondo confinante < 0.01 µm, e, nell'ordine sequenziale, nessuno di loro cambia destinatario rispetto a R, salvo i 2 contatti puntiformi di A6 r1 (RA-code-04).
- I ROI reali reimplementati sono 6 di 40, più un nullo. Il tracciamento degli agganci e dei contatti puntiformi copre tutti i 40 ROI reali; l'arbitro sistematico copre A3 r1 e A5 r1.
