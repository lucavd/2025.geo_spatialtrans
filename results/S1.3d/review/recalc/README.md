# S1.3d — revisione avversariale «ricalcolo» (BL-084)

Revisore: sotto-agente «ricalcolo», 2026-10-15 (sessione S1.3d). Repo `step1-cell-layer` @ `1431d75`, run ufficiale in `/mnt/micron/geo_spatialtrans/S1.3d/`
(sola lettura). Script e log in questa cartella; dati pesanti in `/mnt/micron/geo_spatialtrans/S1.3d/review_recalc/`.
Nessun file del run ufficiale è stato modificato; nessun commit.

## Esito (3 righe)
1. **Eleggibilità e decisione riprodotte indipendentemente**: 20/20 righe di `S1.3d_eligibility.csv` identiche su tutte le colonne (casi, S, silenziosi,
   eccezioni, copertura, celle invalide, CP-1, eleggibile, garantito, t_S, rss_S); decisione per regola **F1 per O1, O2, O3**, come dichiarato.
2. **CP-3 su H confermato due volte**: dai res 0/300 (O1), 0/600 (O2), 0/300 (O3); con un ricalcolo indipendente di F1 dalla verità (codice proprio) lo stesso.
3. **Nessun difetto che cambi il verdetto.** Rilievi: F5-O1 su S rilanciato e ancora in corso dentro la cartella del run ufficiale (non dichiarato; D-2 non
   corrisponde allo stato finale); `S1.3d_verdicts.csv` mappa le previsioni F4/F5 sull'eleggibilità invece che sul testo pre-registrato; difetti minori di documentazione.

## Metodo
| script | cosa fa | log |
|---|---|---|
| `recalc_from_res.R` | Legge i grezzi (`res/*_o1.rds`, `*_o12.rds`, `S_cperf_*.rds`, `coverage/*.json`, `geom/*_meta.json`, `checks/*`, `logs/perf_*.log`, `b1/verdicts_*.csv`), **ri-deriva** i flag (τ_U da `meta$P`/`Preg`, violazioni unione/diff. simmetrica dagli errori grezzi, conteggi delle tabelle `viol` e loro τ ricalcolati), costruisce gli insiemi dai file di chiavi, applica la regola pre-registrata scritta da capo (gerarchia F1, F5, F4-solo-O1; C-6; F1 resta salvo ≥ 2× più veloce e memoria ≤ 1.5×; tempi = mediana delle 3 ripetizioni di preparazione + operazione dai log). Non usa `tools/S1.3d_analysis.R` né i CSV aggregati. | `recalc_from_res.log` |
| `recalc_validity.R` | Lettore parquet proprio; `st_is_valid` degli ingressi per **tutte** le 2 669 chiavi (N, H, U, CP1, S); validità F1/F3 (`geos_set_precision` g = 1e-6) anche su H (assente in `checks/validity_F1F3.csv`); perimetri ricalcolati contro `meta$P`. | `recalc_validity.log` |
| `recalc_independent.R` | **Ricalcolo indipendente dei confronti con la verità**, senza sorgere `tools/S1.3d_*.R`: V0 (sf) e F1 (geos `*_prec`) riscritti dalla definizione; O1 su N ∪ H (2 100 casi × 2) e S (F1); O2/O3 di F1 su O2 ∪ H (420) e S, V0 su S; bin enumerati dal riquadro della cella (non dal predicato), chiavi intere, τ come pre-registrato (τ_U, τ_ij, τ_i, τ_kb = g·perimetro esatto, spurio g·perimetro candidato, O3 = somma dei τ dei 16 sotto-pezzi); elenco coppie ricontrollato con `st_intersects`. | `recalc_independent_*.log` |
| `recalc_compare_indep.R` | Confronto riga per riga ricalcolo indipendente ↔ res. | `recalc_compare_indep.log` |
| `recalc_truth_spot.py` | Controllo a campione dell'arbitro con aritmetica esatta propria (Fraction) sui 24 nulli di K: area di ogni cella, identità del json, intersezione esatta (Sutherland–Hodgman) delle coppie violate da V0 con celle convesse. | `recalc_truth_spot.log` |
| `recalc_seeds.py` | Seed di H (rep 61..64) contro tutte le formule di seed di R3, S1.3, S1.3b (da grep del codice). | `recalc_seeds.log` |
| `recalc_b1_rows.R`, `recalc_misc*.R`, `recalc_vs_official.R` | B-1 riga per riga; dettagli F5/V0/F4; confronto finale con i CSV ufficiali (letti solo qui). | omonimi `.log` |

## Numeri ricalcolati
| quantità | dichiarato | ricalcolato | |
|---|---|---|---|
| Eleggibilità (20 righe op × candidato) | `S1.3d_eligibility.csv` | identica, 0 differenze in 12 colonne | ✓ |
| Decisione | F1, F1, F1 | F1, F1, F1 (classe garantiti eleggibili = {F1} in O1, O2, O3) | ✓ |
| CP-3 F1 su H | 0/300, 0/600, 0/300 | res: 0/300, 0/600, 0/300; indipendente: 0/300, 0/600, 0/300 | ✓ |
| V0 in N, errore silenzioso O1 | 70 (15 unioni, 65 con coppie) | 70 (15 unioni, 65 coppie; 10 entrambe, 5 solo unione); indipendente: 70/15/65 | ✓ |
| Copertura V0 invalida in N | 1108/1800 | 1108/1800 (61.6 %); di cui 20 dei 70 nulli con errore hanno copertura **valida** | ✓ |
| Copertura F4 valida in N | 1800/1800 | 1800/1800 (H: 300/300) | ✓ |
| F2 celle invalide in N | 2, riparate | 2 (A5_r1_CSR_12, A1_r4_RSArule_13), still_invalid 0, 0 errori silenziosi | ✓ |
| C-0 | 1800/1800 | 1800/1800 (n identico; c1, c2_symdiff identici bit per bit, c2_overlap max rel 1.4e-16) | ✓ |
| C-5 | 1986/1986, 0 parità (brief); 2286/2286 nel verdicts finale | 1986/1986 senza H (N 1800, CP1P 180, U5 5, S 1); 2286/2286 con H; 0 parità, 0 ambigui, 0 errori | ✓ |
| CP-1 V0 (O1) | 5/346, esclusi 14 | 5/346 (tutti CP1P), esclusi 14 (CP1T; 15 celle invalide in ingresso, ricalcolate) | ✓ |
| B-1 | 62/62 per candidato | 62/62 previsioni B e 98/98 righe B identiche per ogni candidato; F1/F2/F3 differiscono solo nella riga C-6 di S1.3 (186/1200 → 1/1200), fuori da B | ✓ |
| F4 celle invalide (D-1) | «6 celle in 6 nulli» (D-1, intermedio) | 13 in N (13 nulli A5) + 1 in H; il report in bozza dice 13 | doc |
| Tempi S (mediana di 3) | perf.csv | F1 O1 81.82 s, O2s2 347.97 s; F5 O2s2 697.11 s; F5 O1 **nessuna riga PERF** | ✓ |

## Rilievi (ognuno con evidenza)
**Confermati**
- R-01 Flag dei res coerenti con le quantità grezze: 15 630 righe O1 applicabili, τ_U ricalcolato da meta con max rel 0; 0 discordanze su union/symd/silent/fail;
  tabelle `viol` coerenti con i conteggi (0) e con τ_ij = g(P_i+P_j), τ_i = g·P_i (0 differenze); O2/O3: 0 discordanze fail/silent, nessun max_ratio > 1 con 0 violazioni.
- R-02 Ricalcolo indipendente (4 201 righe O1, 1 266 righe O2/O3): **0 differenze** in union/symd/coppie ij, ji/celle/eccezioni/silent/fail; union_err identico bit per bit
  4201/4201; O2/O3 n_truth, n_cand, n_viol, n_missing, n_spurious, max_ratio identici. Include S: F1 O1/O2/O3 e V0 O2/O3 (quest'ultimo valida il riconfronto D-3: 0 violazioni).
- R-03 Tolleranze: `tools/S1.3d_compare.R` usa τ_U = g(P_reg+ΣP), τ_ij = g(P_i+P_j), τ_kb = g·perimetro esatto, O3 = somma dei τ esatti dei sotto-pezzi
  (`tools/S1.3d_run.R` riga 115; `S1.3d_perf.R` riga 36) con g = `G_GRID` = 1e-6: come la pre-registrazione. Perimetri di meta: max rel 5e-15 dal ricalcolo. Confronto simmetrico (stessa τ e stessa verità per (i,j) e (j,i)).
- R-04 K: 24 nulli = unione esatta dei 15 nulli con unione fallita in S1.3 (c2 > 1e-9) e delle 11 coppie (10 di P1 + 71/175 di A5_r5_RSA_15); tutti i 26 casi sono rilevati da V0 (24/24 nulli silenziosi).
- R-05 Arbitro a campione (K): 200 944 aree di cella esatte = verità (differenza 0.0 in float); 19 coppie violate da V0 con celle convesse: intersezione esatta S-H = verità (diff 0.0); identità union = Σa − Σov esatta.
- R-06 Ingressi: 0 celle invalide in N, H, S (l'analisi ufficiale le assume valide senza verificarle; l'assunzione regge); CP1T 14 casi invalidi = esclusi dichiarati.
- R-07 D-3 ha toccato solo S: tutti i res di N, U, CP1 hanno mtime precedente alla modifica di `compare.R` (2026-10-15 02:04:04); dopo: H (calcolato dopo la decisione, per disegno) e solo `S_cperf_{V0,F2,F3,F4}_O2s{2,8}.rds` + `S_cperf_F2_validity.rds`; i log perf non sono stati sovrascritti (mtime tutte 2026-10-14). F1, F5, F6 su S non riconfrontati (indici interi), coerente con D-3.
- R-08 Seed di H: 300 seed distinti (BASE + 1000·i + {161..164, 261..263, 361..363}); 0 collisioni con i seed di generazione dei nulli di R3, S1.3, S1.3b.
- R-09 Robustezza della decisione: anche senza D-1 (F4 eleggibile in O1) o con F5 completato su S, F1 resta: servirebbe t ≤ 40.9 s (F1 81.82 s); F4 O1 2873 s, F5 O1 > 6 h 40 min. In O2/O3 F5 ha errori silenziosi (O2: 7 righe op×caso in 5 casi di N∪U∪S; CP-1: 50 righe FAIL in O2, 30 in O3).

**Discrepanze / metodo (nessuna cambia il verdetto)**
- R-10 (discrepanza) F5 in O1 su S: nessuno dei 3 log `perf_F5_O1_*` ha una riga PERF (rep 2 e 3 uccise con SIGTERM a 2:47:14 e 24:38; `perf_F5_O1_1.log` è **vuoto, 0 byte, mtime 02:04:53**). La rep 1 è stata **rilanciata** il 2026-10-15 alle 02:04:52 (processo `S1.3d_perf.R F5 O1 1`, figlio di un `cmd.sh` di job; ancora in corso alle 05:25, > 3 h), troncando il log precedente. D-2 dice che la rep 1 «prosegue» e che il tempo «è riportato su una sola ripetizione»: allo stato finale non c'è alcun tempo (t_S = NA, has_S = FALSE). Il rilancio non è dichiarato e, quando finirà, scriverà `res/S_cperf_F5_O1.rds` e il log nella cartella del run ufficiale **dopo** il commit `1431d75`: rieseguire l'analisi cambierebbe l'ingresso (F5-O1 potrebbe diventare has_S = TRUE). La decisione non cambia (R-09). Anche D-3 («F5 non è ≥ 2× più veloce») va letto come «F5 non eleggibile per mancanza di S».
- R-11 (metodo) `S1.3d_verdicts.csv`: le righe P-F2..P-F5 verificano «eleggibile O1», non il testo pre-registrato. F4: la previsione «copertura valida 100 %; 0 errori in O1» è confermata nelle sue parti (1800/1800; 0 silenziosi O1) ma la riga è «smentita» per le 13 celle invalide (clausola di eleggibilità, non della previsione). F5: la previsione «0 errori silenziosi» è smentita in O2/O3 (N: 3 nulli A4 a 2 µm; U5_cocirc; S: 62 pezzi a 2 µm, 7 a 8 µm) ma il file non ha una riga F5-O2/O3; P-F5 risulta «smentita» solo per has_S = FALSE in O1, dove nessun errore è stato misurato. La bozza del report descrive correttamente entrambe; l'esito «chiusa con riserve» è già determinato da P-V0cov (61.6 % < 80 %).
- R-12 (discrepanza, documentazione) D-1 riporta «6 celle F4 invalide in 6 nulli di A5» (letto a 486/1800); il valore finale è 13 celle in 13 nulli A5 di N più 1 in H (A5_r4_RSArule_63). Il report usa 13; le dev_notes andrebbero annotate.
- R-13 (metodo, minore) Seed di H: 57 coincidenze numeriche con seed di **altri usi** di S1.3b (10 con `st_sample` di C-3, fino a 51 con `s3b_perturb` di CP-2 — limite superiore, solo per i 60 job selezionati), dovute a SEL_SEED = BASE + 2. Flussi RNG usati per operazioni diverse su ROI diverse: nessun effetto plausibile su H, ma «base di seed mai usata» è vera solo per i seed di generazione.
- R-14 (metodo, minore) Commento di testa di `tools/S1.3d_compare.R` (riga 4: τ di O3 «approssimata») non descrive il codice, che calcola la somma esatta dei τ dei sotto-pezzi (conforme). `S1.3d_perf.csv`: le colonne t_prepare e t_op sono mediane separate e non sommano a t_total (F4 O1: 1660.29 + 1216.18 = 2876.47 contro 2873.20); la regola usa t_total (mediana dei totali), corretto.
- R-15 (metodo, potenziale) `cmp_o2` porta a 0 le aree NA/NaN del candidato: un NaN su un pezzo esatto ≤ τ passerebbe in silenzio. Nessun caso osservato (F1 e V0 ricalcolati identici, n_cand coerente).
- R-16 (metodo, minore) `checks/validity_F1F3.csv` non copre H (2 369 chiavi); ricalcolato: 0 celle F1 invalide in H (e in tutte le 2 669 chiavi; 1 cella vuota in U3_sliver, area esatta entro τ).

## Deviazioni di questa revisione
- Il ricalcolo indipendente dalle geometrie copre V0 e F1 (O1 su N ∪ H ∪ S per F1; O2/O3 di F1 su O2 ∪ H ∪ S; V0 O2/O3 su S). F2, F3, F4, F5, F6 e i casi U/CP1 sono
  ricontrollati solo dai res (flag ri-derivati dalle quantità grezze), non rieseguiti.
- L'arbitro (verità esatta) non è stato ricalcolato per intero: controllo esatto a campione su K (aree delle celle, identità, 19 coppie convesse).
