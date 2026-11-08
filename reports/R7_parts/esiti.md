### Cancelli

| id | esito | dettaglio |
|---|---|---|
| G-1 | superato | batteria 10/10, mutanti rilevati 7/7 (anche dopo C-3) |
| G-2 | superato | 79/79 controlli; F_box contro i valori di E-0 entro 0.0049 (intestino r3-r5, 0.3499 contro 0.345 scritto); ogni bin unico in un solo pezzo; nel cuore valori sui pezzi identici |
| G-3 | superato | 342/342: quantità per area di R4 e intervalli di `R4_refs.csv` riprodotti, scarto relativo massimo 2.5e-12 |
| G-4 | superato | 248 file letti, md5 identici prima e dopo (`results/R7/R7_md5_before.txt`, `R7_md5_after.txt`) |
| G-5 | superato dopo le correzioni C-1...C-6 | quattro revisori indipendenti (sezione «Revisione avversariale»); i ricalcoli tornano con le correzioni |

### Regole

- **R-1, la correzione**: allargamento calibrato (C-a), per pareggio (12 celle tessuto × scenario su 18 dentro Bradley per C-a e per C-b, anche con i κ di C-1; `R7_P1_cells_v3.csv`). Dopo P-2 la regola si applica per tessuto: **stroma, linfonodo e cervello** hanno almeno un campo smentito e prendono κ dalla struttura misurata (campo UMI per le quantità per bin, campo nuclei per quelle dei nuclei, il maggiore dei due per le altre).
- **R-2**: la condizione scatta (stroma r2-r5: quota di nuclei condivisi 0.168 contro quota di area 0.280); si usano le matrici specifiche per bin e per nuclei in tutti i tessuti (lettura dichiarata, deviazione DV-2).
- **R-3**: applicata alle righe a due stadi di R4; non alle finestre (C-5).

Fattori finali (bootstrap percentile; `R7_kappa_final_v4.csv`):

| tessuto | per bin | per nuclei | altre | ± 2 SE, per bin / nuclei | minimo-massimo, per bin / nuclei | fonte |
|---|---|---|---|---|---|---|
| intestino | 1.857 | 1.854 | 1.837 | 1.586 / 1.583 | 1.033 / 1.033 | geometria |
| tumore del colon | 1.655 | 1.666 | 1.662 | 1.413 / 1.422 | 1.004 / 1.005 | geometria |
| stroma del tumore | 3.193 | 2.817 | 3.193 | 2.740 / 2.417 | 1.253 / 1.188 | struttura misurata |
| linfonodo | 2.571 | 2.966 | 2.966 | 2.206 / 2.550 | 1.147 / 1.216 | struttura misurata |
| cervello di topo | 1.818 | 3.865 | 3.865 | 1.555 / 3.324 | 1.029 / 1.372 | struttura misurata |
| cuore | 1.624 | 1.627 | 1.625 | 1.387 / 1.387 | 1.000 / 1.000 | geometria |
| aree simulate (5) | 1.632 | | | 1.391 | 1 | poche unità |
| finestre a mano (4) | 1.979 | | | 1.594 | 1 | poche unità |

### Previsioni

| id | tessuto | atteso | intervallo di rilevanza (fonte) | osservato | esito |
|---|---|---|---|---|---|
| P-1 | sei tessuti, valori lognormali con CV 0.624 | copertura dentro | 0.925-0.975 (Bradley) | 0.907-0.912 (± 0.004) | **smentita** in 6 su 6 |
| P-1 | sei tessuti, valori t con 3 gradi | copertura dentro | 0.925-0.975 (Bradley) | 0.959-0.966 (± 0.003) | confermata in 6 su 6 |
| P-2 | intestino | dentro per i due campi | 0.925-0.975 (Bradley) | UMI 0.947 [0.599, 0.949]; nuclei 0.949 [0.947, 0.948] | UMI indeterminata; nuclei confermata |
| P-2 | tumore del colon | dentro | 0.925-0.975 | UMI 0.942 [0.935, 0.946]; nuclei 0.948 [0.940, 0.951] | confermata, confermata |
| P-2 | stroma del tumore | dentro | 0.925-0.975 | UMI 0.846 [0.796, 0.893]; nuclei 0.884 [0.814, 0.930] | **smentita**; indeterminata |
| P-2 | linfonodo | dentro | 0.925-0.975 | UMI 0.900 [0.853, 0.909]; nuclei 0.864 [0.579, 0.917] | **smentita**, **smentita** |
| P-2 | cervello di topo | dentro | 0.925-0.975 | UMI 0.943 [0.739, 0.945]; nuclei 0.742 [0.736, 0.922] | indeterminata; **smentita** |
| P-3 | cinque tessuti, 7 quantità | scarto entro la soglia | ± metà della semi-ampiezza corretta | rapporto massimo 0.50 (linfonodo, rapporto nucleo/corona); scarto relativo massimo 7.3 % (linfonodo, frazione di bin vuoti a 2 µm) | confermata in 35 su 35 |

P-1 e P-2 sono smentite come scritte (P-2 attendeva 5 tessuti su 5). Nel cervello, campo nuclei, il modello del semivariogramma ha la scala al limite superiore (L = 2 000 µm): la densità dei nuclei non si stabilizza entro 720 µm e la correlazione fra aree lontane è un'estrapolazione; il bordo superiore della banda è 0.922 contro 0.925. Nell'intestino e nel cervello, campo UMI, l'indeterminazione viene da repliche bootstrap con L al limite; le stime puntuali (0.947, 0.943) sono dentro. Tabelle: `results/R7/coverage/R7_P1_cells_v3.csv`, `results/R7/cp1_v2/R7_P2.csv`, `results/R7/pieces/R7_P3.csv`.

### Controprove

- **CP-1 = P-2**: fallisce nello stroma, nel linfonodo e nel cervello: lì la sola geometria non basta, e la regola R-1 calibra κ sulla struttura misurata. La struttura ha scala 270-480 µm nello stroma e nel linfonodo (figura 3); le aree che si toccano sono correlate 0.28-0.73 invece della quota in comune 0.08-0.39 (figura 2).
- **CP-2 = P-3**: superata. Si correggono gli intervalli, non i valori.
- **CP-3, cuore**: superata (G-2; κ_range 1 esatto dopo C-3).

### Misure

- **M-a, copertura del metodo di oggi** (20 000 repliche, `R7_kappa.csv`): bootstrap sulla media di 5 aree 0.798 intestino, 0.835 tumore, 0.760 stroma, 0.770 linfonodo, 0.816 cervello, 0.842 cuore; 5 aree indipendenti 0.837; 4 finestre 0.796. Minimo-massimo come intervallo della media: da 0.895 (stroma) a 0.940 (cuore), 0.872 con 4 finestre.
- **M-b, verdetti passati** (`results/R7/report/R7_transitions_by_step.csv`): 7 852 righe di 29 step; 6 624 con almeno una componente corretta; 709 cambiano fra i tre esiti (confermata → indeterminata 168, smentita → indeterminata 197, smentita → confermata 175, indeterminata → confermata 147, indeterminata → smentita 22, confermata → smentita 0); 20 cambi di regole o controprove (`R7_rule_changes.csv`). Tutti i verdetti confrontabili si riproducono dai numeri delle tabelle.
- **M-c, riferimenti corretti** (`results/R7/Mc/R7_Mc_refs.csv`): 386 componenti di 200 righe di riferimento (bootstrap 209, minimo-massimo 128, finestre 18, lascia-fuori 27, fisse 4); i 42 bersagli del confronto finale in `results/R7/report/R7_S25_targets.csv`.
- **M-d, nuclei condivisi nelle validazioni lasciando fuori un'area** (`R7_Md_leakage.csv`): fino a 0.596 (linfonodo f2), 0.570 (linfonodo p2), 0.478 (stroma r2); tumore al più 0.099; cuore 0.
- **M-e**: deviazione standard fra pezzi disgiunti / fra aree da 0.90 a 1.59 (`R7_P3.csv`, colonna `sd_ratio`).
- **M-f**: semi-ampiezza a due stadi / uno stadio, mediana degli UMI per cellula: 0.86-1.23 (`R7_Mf_two_stage.csv`).
