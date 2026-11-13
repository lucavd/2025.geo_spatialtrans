
**PH-1, mescolanza dei tipi fra vicini** (`tools/F3.1_PH1_mixing.py`, `results/F3.1/F3.1_PH1_mixing*.csv`; figura 4). Per ogni nucleo con tipo, quota dei 6 nuclei con tipo più vicini che sono di un altro tipo; stessa definizione sui due lati. Nel simulato è più alta in tutti i tessuti; dentro ogni lato i marcatori estranei crescono con la mescolanza (tranne il linfonodo, dove i marcatori non separano i linfociti B dai T). Ripesando le cellule simulate per livello di mescolanza con la distribuzione del vero, si spiega gran parte dello scarto di M-3 in intestino e tumore:

| tessuto | vicini di altro tipo, vero | simulato | M-3 vero | M-3 simulato | M-3 simulato ripesato | scarto spiegato |
|---|---|---|---|---|---|---|
| intestino | 0.22 | 0.47 | 0.240 | 0.363 | 0.247 | 95% |
| tumore del colon | 0.06 | 0.35 | 0.045 | 0.151 | 0.082 | 65% |
| stroma | 0.27 | 0.69 | 0.216 | 0.332 | 0.307 | 22% |
| linfonodo | 0.36 | 0.58 | 0.674 | 0.765 | 0.768 | -3% |
| cervello | 0.40 | 0.47 | 0.236 | 0.307 | 0.291 | 22% |
| cuore | 0.31 | 0.40 | 0.169 | 0.203 | 0.183 | 59% |

**PH-2, composizione dei primi tipi di RCTD** (`F3.1_PH1_mixing_rctd_composition.csv`). Nel simulato RCTD chiama il tipo dominante più spesso del vero: intestino epitelio 0.87 contro 0.67 e stroma (lamina propria e muscolo) 0.08 contro 0.25, benché nel simulato le cellule stromali siano il 21 % come fra i nuclei veri; stroma del tumore fibroblasti 0.23 contro 0.60; cervello neuroni 0.75 contro 0.55; linfonodo plasmacellule 0.30 contro 0.13 e linfociti T 0.21 contro 0.42.

**PH-3, semina a compartimenti** (livelli Q6 di S1.6, 30 aree × 1 seme, stessi semi di k = 0; `tools/F3.1_run_PH3.sh`, `tools/F3.1_PH3_analysis.py`): la semina a compartimenti riporta la mescolanza dei vicini al livello del vero (quota dei vicini di altro tipo, semina uniforme → compartimenti: intestino 0.47 → 0.23 (vero 0.22); tumore del colon 0.35 → 0.10 (vero 0.06); stroma 0.69 → 0.31 (vero 0.27); linfonodo 0.58 → 0.37 (vero 0.36); cervello 0.47 → 0.36 (vero 0.40); cuore 0.40 → 0.32 (vero 0.31)). Le misure dei metodi si avvicinano al vero in intestino, tumore, stroma e linfonodo: fra i giudizi smentiti con la semina uniforme, lo scarto si chiude dal 37 al 95 % per singoletti, concordanza e coerenza spaziale (nell'intestino la concordanza passa oltre il vero), e di circa la metà per i marcatori estranei (intestino 0.52, tumore 0.59; linfonodo 0.17; nello stroma, non smentito, 0.48); 5 dei 17 giudizi smentiti fuori dal cuore non lo sono più (posizione descrittiva, un seme). Restano: UMI assegnati in eccesso in intestino e tumore (densità: con i compartimenti peggiorano appena), marcatori estranei nell'intestino (0.298 contro 0.241) e nel linfonodo, il cervello (mescolanza già vicina al vero; nessuna misura si muove), e la concordanza dell'intestino che supera il vero (0.907 contro 0.843-0.882). UMI per quadratino e cellule per area invariati. **Anomalia**: nell'area r3 del cervello la simulazione a pezzi sul livello a compartimenti si ferma («.s27_cdt(): triangolazione fallita anche dopo la semplificazione (geometria 168)»): il cervello ha 4 aree su 5 in PH-3 (voce di BACKLOG).

| misura | tessuto | vero | semina uniforme | compartimenti | frazione dello scarto chiusa | esito con la semina uniforme | posizione dei compartimenti rispetto al vero (descrittiva) |
|---|---|---|---|---|---|---|---|
| M-1 | intestino | 0.676 | 0.743 | 0.756 | -0.19 | smentita | smentita |
| M-2 | intestino | 470 | 447 | 453 | 0.26 | confermata | confermata |
| M-3 | intestino | 0.241 | 0.361 | 0.298 | 0.52 | smentita | smentita |
| M-4 | intestino | 1.000 | 0.968 | 0.980 | 0.37 | smentita | smentita |
| M-5 | intestino | 0.864 | 0.836 | 0.907 | 2.50 | indeterminata | smentita |
| M-6 | intestino | 0.909 | 0.819 | 0.854 | 0.40 | smentita | smentita |
| M-1 | tumore del colon | 0.744 | 0.817 | 0.827 | -0.14 | smentita | smentita |
| M-2 | tumore del colon | 1136 | 1094 | 1101 | 0.18 | confermata | confermata |
| M-3 | tumore del colon | 0.045 | 0.152 | 0.088 | 0.59 | smentita | smentita |
| M-4 | tumore del colon | 0.954 | 0.815 | 0.947 | 0.95 | smentita | indeterminata |
| M-5 | tumore del colon | 0.990 | 0.922 | 0.977 | 0.82 | smentita | smentita |
| M-6 | tumore del colon | 0.936 | 0.762 | 0.917 | 0.89 | smentita | indeterminata |
| M-1 | stroma | 0.586 | 0.421 | 0.495 | 0.45 | confermata | confermata |
| M-2 | stroma | 669 | 657 | 763 | 8.47 | confermata | confermata |
| M-3 | stroma | 0.222 | 0.334 | 0.281 | 0.48 | indeterminata | confermata |
| M-4 | stroma | 0.985 | 0.848 | 0.934 | 0.62 | smentita | smentita |
| M-5 | stroma | 0.890 | 0.802 | 0.873 | 0.80 | smentita | confermata |
| M-6 | stroma | 0.886 | 0.647 | 0.791 | 0.60 | smentita | smentita |
| M-1 | linfonodo | 0.984 | 0.998 | 0.990 | 0.62 | confermata | confermata |
| M-2 | linfonodo | 102 | 103 | 102 | 1.15 | confermata | confermata |
| M-3 | linfonodo | 0.672 | 0.765 | 0.749 | 0.17 | smentita | smentita |
| M-4 | linfonodo | 0.999 | 0.993 | 0.994 | 0.10 | confermata | confermata |
| M-5 | linfonodo | 0.600 | 0.405 | 0.510 | 0.54 | smentita | indeterminata |
| M-6 | linfonodo | 0.472 | 0.339 | 0.413 | 0.56 | smentita | indeterminata |
| M-1 | cervello | 0.387 | 0.348 | 0.361 | 0.33 | confermata | confermata |
| M-2 | cervello | 4514 | 3298 | 3399 | 0.08 | confermata | confermata |
| M-3 | cervello | 0.232 | 0.306 | 0.304 | 0.04 | smentita | smentita |
| M-4 | cervello | 1.000 | 1.000 | 0.999 | 1.14 | confermata | confermata |
| M-5 | cervello | 0.904 | 0.944 | 0.935 | 0.22 | indeterminata | indeterminata |
| M-6 | cervello | 0.775 | 0.824 | 0.820 | 0.09 | smentita | smentita |
| M-1 | cuore | 0.177 | 0.157 | 0.163 | 0.31 | confermata | confermata |
| M-2 | cuore | 150 | 143 | 150 | 1.05 | confermata | confermata |
| M-3 | cuore | 0.164 | 0.199 | 0.197 | 0.07 | smentita | indeterminata |
| M-4 | cuore | 1.000 | 1.000 | 1.000 |  | confermata | confermata |
| M-5 | cuore | 0.809 | 0.844 | 0.855 | -0.28 | indeterminata | indeterminata |
| M-6 | cuore | 0.750 | 0.708 | 0.733 | 0.60 | confermata | confermata |

