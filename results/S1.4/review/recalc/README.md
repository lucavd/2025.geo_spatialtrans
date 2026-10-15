# S1.4 — revisione avversariale «recalc» (ricalcolo indipendente)

Revisore: sotto-agente `recalc`, 2026-10-15, commit `ec2a874` (branch `step1-cell-layer`). Non committato.
Script e output: `/mnt/micron/geo_spatialtrans/S1.4/review_recalc/` (`recalc.R`, `verdicts.R`, `robust.R`, `seed_rank.R`,
`export_geom.R`, `geom_check.py`, `poi_check.py`; CSV/JSON/log con lo stesso nome).
Codice proprio: OLS con `lm()`, KS con `ecdf()`, Spearman come Pearson dei ranghi, simulazione del cerchio riscritta
dalla specifica (D-S1.4.1/D-S1.4.2: `set.seed(seme); z <- rnorm(n)`, target, k32, clamp `0.85·d_POI`); geometria in
Python/shapely 2.1.2 (GEOS `maximum_inscribed_circle`, `STRtree`). Nessuna funzione di `tools/S1.4_*.R`;
`derive_cell_geometry()` non è stata chiamata: la riscrittura riproduce esattamente i suoi output (punto 2).

## Esito sintetico
Tutti i numeri chiave ricalcolati coincidono con quelli riportati; nessun verdetto cambia. Rilievi elencati sotto.

| # | Quantità | Riportato | Ricalcolato | Esito |
|---|---|---|---|---|
| 1 | Calibrazione (ρ, β, σ, A_ref, n, β_perm, mediana N/C, ecc), 8 righe (6 primari + StarDist A4/A6) | `S1.4_calibration.csv` | max \|Δ\| ≤ 1e-15; n identici (33203 … 2610) | conferma |
| 1b | se_beta | es. A1 0.0042498 | Δ rel ≈ 1/(2n) (1.5e-5 A1; 2.5e-4 A5): denominatore diverso di 1 | conferma (irrilevante) |
| 2 | LOO: 160 righe (6 archetipi × 5 fold × M0–M3, + StarDist A4/A6) | `S1.4_loo_folds.csv` | KS(N/C), KS(A_N), clamp, metro (mediana, max): Δ = 0; \|ΔSpearman\|: max 2e-8; ρ, β, σ: ≤ 1e-15; semi = 20261015+1000i+10k (40/40) | conferma |
| 3 | B-1, B-2, B-3 (stringhe di fold e verdetto), B-4, CP-1, CP-3, P1–P5 | `S1.4_verdicts.csv` | identici (18/18 B-1..B-3, 6/6 B-4, CP-1 6/6 + «basta» M1: A2, A5; M3: A3, A6; P1 6/6, P2 6/6, P3 0/6, P4 FAIL, P5 5/6) | conferma |
| 4 | PH-1 circ_poi_085 A1 r1–r5 | 0.17988 0.17074 0.17491 0.19221 0.17795 | d_POI salvata: identici; d da GEOS MIC indipendente: 0.17988 0.17074 0.17491 0.19239 0.17795 | conferma |
| 4b | PH-1 circ_poi_085 A4 f1 f2 p1 p2 p3 | 0.37681 0.41632 0.47602 0.50099 0.46981 | GEOS MIC: 0.37687 0.41621 0.47598 0.50103 0.46972 (\|Δ\| ≤ 1.1e-4) | conferma |
| 5 | C8 [q05, q95], spec, ×2, ×0.5 — tutti e 6 gli archetipi | `S1.4_c8.csv` | identici (es. A1 [2.8104, 7.7138] 5/5/5; A5 [6.9535, 26.8783] 5/4/1) | conferma |
| 6 | BL-060 appaiamento (shapely within e intersects), A4 f1 f2 p1 p2 p3, A6 r1 r3 | A4 f1 n_terr_0/1/2p = 6892/16706/109; A6 r1 272/325/77 | identici in tutti i 7 set; vettore `matched` concorde 100 %; a_n identica; frazione appaiata delle interne A4 0.708–0.759, A6 0.472–0.770 (dev notes «71–76 %, 47–77 %») | conferma |
| 7 | Aree territori real/*.rds vs R3 `*_cells.parquet` (tutti i 40 set) | (atteso: coincidenti) | **non coincidono** per 0–305 celle/set (max rel 2.32); Σ aree identica; vedi R-1 | rilievo |

## Rilievi
- **R-1 (rilievo) territori ≠ R3.** In 36/40 set (0 in A3 r1, A3 r3, A4 p3 spaceranger e stardist_he) alcune celle hanno area diversa da R3 (rel > 1e-6): soprattutto non interne,
  ma anche interne in entrambi (fino a 17/set, A6 r5; Δ rel fino a 0.77 in A1 r2) e celle dell'insieme di valutazione (fino a 21/set, A1 r1 = 0.28 %).
  Le interne del pacchetto sono 0–16 in più di R3 per set (mai meno; A1 r1 +13, BL-071 riporta «+12»). Σ aree identica
  (A1 r1 798307.54051 µm² in entrambi) → riassegnazione dei frammenti D-S1.3.2 (BL-071, aperta), dichiarata nella pre-registrazione
  («territori S1.3»). Non è un bug di S1.4, ma il report non deve dire che i territori sono quelli di R3. File `recalc_r3_area*.csv`.
- **R-2 (deviazione) sensibilità al seme non riportata.** Le dev notes dicono «10 semi aggiuntivi per M2, solo riportata»: i valori sono in
  `loo/*.rds$seeds` ma non sono usati da `S1.4_analysis.R` né citati in `S1.4.qmd`. Ricalcolo dei verdetti d'archetipo sui 10 semi:
  A5 B-1 FAIL (primario) → WARN in 9/10 semi, FAIL 1/10; A6 B-2 PASS 5/10, WARN 5/10; A6 B-3 WARN 7/10, PASS 3/10; tutti gli altri 10/10 uguali al primario.
  Il seme primario non è anomalo (rango medio fra 11: KS(N/C) 6.7, KS(A_N) 6.97, atteso 6, SE 0.58). P2/P3 non cambiano.
- **R-3 (rilievo) CP-2 sottopotenziato in A5/A6.** Nullo per permutazione (200 permutazioni): sd(β_perm) 0.0246 (A5), 0.0200 (A6);
  P(\|β\| > 0.02 \| nullo) = 0.435 (A5), 0.320 (A6); 0 in A1–A4. Il FAIL di A5 (0.0224) è compatibile col caso.
- **R-4 (rilievo) POI di polylabelr subottimale nei territori con buchi.** Su 12 set controllati, 10 territori (A4 f1 1, A6 r1 8, A6 r3 1), tutti con buchi
  della maschera e tutti non interni (0 nell'insieme di valutazione): d_POI salvata più piccola del raggio inscritto massimo GEOS/shapely
  fino a 5.87 µm (A6 r1 cella 144: 20.24 vs 26.06). Effetto solo sul clamp di celle di bordo (≤ 8/674 in A6 r1). Provvisorio: 12/40 set.
- **R-5 (rilievo) A_ref = Σ A_T/n su tutte le celle = area finestra/n**, più grande della media delle interne (A1 94.4 vs 78.5; A6 940 vs 855 µm²).
  La stima e la simulazione sono invarianti (stesso A_ref nei due passi), quindi i verdetti no; ma S1.5 deve definire la densità allo stesso modo
  (n/area della finestra con lumi) perché (ρ, β) siano trasferibili.
- **R-6 (rilievo) celle reali con A_N ≥ A_T**: 730/186 983 nell'insieme di valutazione (0.39 %; fino a 1.9 % in A2 r4; rapporto max 1.96), irraggiungibili
  dal simulatore; effetto sulle KS ≤ frazione, nessun verdetto cambia. «N/C» = A_N/A_T (definito nelle dev notes; non è nucleo/citoplasma).
- **R-7 (rilievo) riproducibilità**: `real/*.rds$cells$idx` è un vettore ALTREP di arrow; letto senza arrow nel percorso di librerie torna a lunghezza zero
  (avviso «cannot unserialize ALTVEC»).
- PH-1: le colonne `*_085` usano il raggio circoscritto (k32), le `*_1` il raggio equivalente: definizioni diverse fra colonne (minore).
