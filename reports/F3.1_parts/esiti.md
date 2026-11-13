
Tabelle: `results/F3.1/F3.1_verdicts.csv` (una riga per giudizio), `F3.1_CP1.csv`, `F3.1_CP2.csv`, `F3.1_CP3.csv`, `F3.1_R1.csv`, `F3.1_R1_summary.csv`, `F3.1_truth.csv`; per area e per esecuzione `F3.1_real_runs.csv`, `F3.1_sim_runs.csv`, `F3.1_sim_areas.csv`, `F3.1_simnuc_runs.csv`. Figure in `reports/F3.1_parts/`.

### Cancelli

| id | esito | dettaglio |
|---|---|---|
| G-1 | passato | 209 asserzioni, 0 fallite (`results/F3.1/F3.1_G1.log`); la prima esecuzione si era interrotta per un errore nel codice delle prove (riga commentata per metà), corretto prima del rilancio |
| G-2 | passato | 581 controlli su 115 aree (30 vere, 60 simulate, 25 con i profili dal nucleo) e i 6 riferimenti, 0 falliti; un controllo della bozza era una tautologia ed è stato sostituito (vedi Deviazioni) |
| G-3 | passato dopo una correzione | prima esecuzione: RCTD non riproducibile (stessa area simulata dell'intestino: singoletti 0.9672 contro 0.9687; esportazione e bin2cell identici). Causa: RCTD estrae a caso i blocchi su cui stima il rumore e lo script non fissava il seme. Correzione: seme dichiarato 20261118 prima di RCTD, due esecuzioni della stessa area identiche, RCTD rifatto su tutte le aree. Dopo: 72 confronti su 72 identici (un'area per tessuto, simulata da zero e vera, 3 processi invece di 6, più i file del riferimento) |
| G-4 | passato | md5 identici prima e dopo per i 270 file letti (`results/F3.1/F3.1_md5_before.txt`, `_after.txt`) |
| G-5 | passato | due revisori indipendenti (sotto, «Revisione avversariale»): 78 misure su 78 entro 1e-9 (scarto massimo 1e-15); 1 119 quantità statistiche su 1 121; le 2 diverse sono esiti al limite spiegati dal caso del bootstrap |

### Previsioni

Conteggio per misura (36 giudizi: 14 confermati, 4 indeterminati, 18 smentiti):

| misura | confermata | indeterminata | smentita |
|---|---|---|---|
| M-1 | 4 | 0 | 2 |
| M-2 | 6 | 0 | 0 |
| M-3 | 0 | 1 | 5 |
| M-4 | 3 | 0 | 3 |
| M-5 | 0 | 3 | 3 |
| M-6 | 1 | 0 | 5 |

| misura | tessuto | vero | intervallo di rilevanza | simulato | intervallo del simulato | esito |
|---|---|---|---|---|---|---|
| M-1 | intestino | 0.676 | 0.667-0.684 | 0.743 | 0.738-0.747 | smentita per eccesso |
| M-2 | intestino | 470 | 414-526 | 446 | 440-453 | confermata |
| M-3 | intestino | 0.241 | 0.226-0.255 | 0.363 | 0.360-0.367 | smentita per eccesso |
| M-4 | intestino | 1.000 | 0.990-1.000 | 0.969 | 0.966-0.972 | smentita per difetto |
| M-5 | intestino | 0.864 | 0.842-0.882 | 0.835 | 0.825-0.843 | indeterminata |
| M-6 | intestino | 0.909 | 0.898-0.918 | 0.819 | 0.814-0.824 | smentita per difetto |
| M-1 | tumore del colon | 0.744 | 0.717-0.776 | 0.817 | 0.815-0.820 | smentita per eccesso |
| M-2 | tumore del colon | 1136 | 923-1352 | 1096 | 1071-1117 | confermata |
| M-3 | tumore del colon | 0.045 | 0.035-0.057 | 0.152 | 0.149-0.154 | smentita per eccesso |
| M-4 | tumore del colon | 0.954 | 0.940-0.968 | 0.815 | 0.809-0.820 | smentita per difetto |
| M-5 | tumore del colon | 0.990 | 0.987-0.993 | 0.923 | 0.918-0.927 | smentita per difetto |
| M-6 | tumore del colon | 0.936 | 0.920-0.951 | 0.763 | 0.757-0.769 | smentita per difetto |
| M-1 | stroma | 0.586 | 0.389-0.713 | 0.419 | 0.416-0.423 | confermata |
| M-2 | stroma | 669 | 403-837 | 653 | 642-668 | confermata |
| M-3 | stroma | 0.222 | 0.118-0.334 | 0.333 | 0.330-0.335 | indeterminata |
| M-4 | stroma | 0.985 | 0.971-0.999 | 0.849 | 0.843-0.854 | smentita per difetto |
| M-5 | stroma | 0.890 | 0.835-0.940 | 0.804 | 0.795-0.812 | smentita per difetto |
| M-6 | stroma | 0.886 | 0.845-0.929 | 0.649 | 0.644-0.654 | smentita per difetto |
| M-1 | linfonodo | 0.984 | 0.969-0.999 | 0.998 | 0.998-0.998 | confermata |
| M-2 | linfonodo | 102 | 56-149 | 103 | 102-103 | confermata |
| M-3 | linfonodo | 0.672 | 0.630-0.708 | 0.765 | 0.763-0.767 | smentita per eccesso |
| M-4 | linfonodo | 0.999 | 0.989-1.000 | 0.993 | 0.993-0.994 | confermata |
| M-5 | linfonodo | 0.600 | 0.456-0.700 | 0.406 | 0.403-0.409 | smentita per difetto |
| M-6 | linfonodo | 0.472 | 0.421-0.518 | 0.340 | 0.338-0.341 | smentita per difetto |
| M-1 | cervello | 0.387 | 0.326-0.453 | 0.349 | 0.343-0.354 | confermata |
| M-2 | cervello | 4514 | 1534-6564 | 3344 | 3239-3435 | confermata |
| M-3 | cervello | 0.232 | 0.179-0.285 | 0.307 | 0.305-0.309 | smentita per eccesso |
| M-4 | cervello | 1.000 | 0.990-1.000 | 1.000 | 1.000-1.000 | confermata |
| M-5 | cervello | 0.904 | 0.873-0.941 | 0.943 | 0.939-0.947 | indeterminata |
| M-6 | cervello | 0.775 | 0.754-0.796 | 0.824 | 0.820-0.827 | smentita per eccesso |
| M-1 | cuore | 0.177 | 0.139-0.211 | 0.158 | 0.152-0.163 | confermata |
| M-2 | cuore | 150 | 80-226 | 143 | 137-151 | confermata |
| M-3 | cuore | 0.164 | 0.136-0.192 | 0.203 | 0.199-0.207 | smentita per eccesso |
| M-4 | cuore | 1.000 | 0.990-1.000 | 1.000 | 1.000-1.000 | confermata |
| M-5 | cuore | 0.809 | 0.783-0.832 | 0.847 | 0.830-0.862 | indeterminata |
| M-6 | cuore | 0.750 | 0.685-0.816 | 0.708 | 0.698-0.721 | confermata |

Contro le previsioni scritte: P-2 e P-4 come previsto in tutti i tessuti; P-1 smentita anche nel tumore del colon (oltre l'intestino); P-3 smentita anche nel cervello e nel cuore; P-5 indeterminata nell'intestino (prevista smentita); P-6 smentita per eccesso nel cervello (prevista indeterminata) e confermata nel cuore (prevista indeterminata). Gli scarti delle smentite vanno da 1.4 a 25 semi-ampiezze dell'intervallo del vero (`gap_in_halfwidths`).

Osservazione non prevista (descrittiva): le 5 aree simulate di un tessuto variano fra loro molto meno delle 5 aree vere (figura 1: punti blu quasi sovrapposti); per questo l'intervallo del simulato è stretto, e una conferma vuol dire «il simulato cade nella variabilità del vero», non che ne riproduca la variabilità.

### Controprove

- **CP-1 (le misure vedono gli errori di assegnazione): passata** in tutti i tessuti e su entrambi i lati per M-3 e M-5 (intestino M-3 real 0.179 ± 0.004; intestino M-5 real 0.341 ± 0.005; tumore del colon M-3 real 0.150 ± 0.018; tumore del colon M-5 real 0.175 ± 0.018; stroma M-3 real 0.135 ± 0.009; stroma M-5 real 0.197 ± 0.014; linfonodo M-3 real 0.057 ± 0.006; linfonodo M-5 real 0.187 ± 0.007; cervello M-3 real 0.224 ± 0.013; cervello M-5 real 0.336 ± 0.010; cuore M-3 real 0.094 ± 0.030; cuore M-5 real 0.093 ± 0.036; intestino M-3 sim 0.143 ± 0.002; tumore del colon M-3 sim 0.236 ± 0.001; stroma M-3 sim 0.375 ± 0.003; linfonodo M-3 sim 0.072 ± 0.001; cervello M-3 sim 0.154 ± 0.003; cuore M-3 sim 0.161 ± 0.007; intestino M-5 sim 0.130 ± 0.002; tumore del colon M-5 sim 0.233 ± 0.002; stroma M-5 sim 0.509 ± 0.005; linfonodo M-5 sim 0.109 ± 0.002; cervello M-5 sim 0.276 ± 0.003; cuore M-5 sim 0.245 ± 0.010).
- **CP-2 (le misure senza verità seguono la verità): passata, debolmente.** Correlazione di Spearman sulle 10 esecuzioni simulate per tessuto, positiva per M-3 in 4 tessuti su 6 e per M-5 in 5 su 6 (criterio: almeno 4). Con 10 punti per tessuto le correlazioni sono rumorose; in cervello e cuore M-3 non segue la verità.

| tessuto | M-3 contro verità (Spearman) | M-5 contro verità (Spearman) |
|---|---|---|
| intestino | 0.85 | -0.12 |
| tumore del colon | 0.85 | 0.71 |
| stroma | 0.16 | 0.71 |
| linfonodo | 0.85 | 0.21 |
| cervello | -0.26 | 0.36 |
| cuore | -0.08 | 0.02 |

- **CP-3 (contaminazione contata due volte): fallita, quindi l'ipotesi di E-3 è smentita.** Con i profili dal solo nucleo la quota di marcatori estranei del simulato si avvicina al vero di almeno metà dello scarto solo nel cervello (1 tessuto su 4, criterio 3 su 4). Anche i profili veri generati dal solo nucleo hanno quasi la stessa quota di marcatori estranei di quelli con la corona (`F3.1_truth.csv`, colonna `foreign_profile`: intestino 0.221 contro 0.224, tumore 0.037 contro 0.044, linfonodo 0.634 contro 0.645). L'UMI per quadratino resta invariato (intestino 29.15 contro 29.14; tumore del colon 47.41 contro 47.41; linfonodo 11.56 contro 11.56; cervello 44.32 contro 44.32; cuore 3.60 contro 3.60).

| tessuto | misura | vero | simulato (k = 0) | profili dal solo nucleo | frazione dello scarto chiusa |
|---|---|---|---|---|---|
| intestino | M-3 | 0.241 | 0.361 | 0.369 | -0.06 |
| intestino | M-5 | 0.864 | 0.836 | 0.839 | 0.13 |
| intestino | M-6 | 0.909 | 0.819 | 0.819 | 0.00 |
| tumore del colon | M-3 | 0.045 | 0.152 | 0.140 | 0.11 |
| tumore del colon | M-5 | 0.990 | 0.922 | 0.921 | -0.01 |
| tumore del colon | M-6 | 0.936 | 0.762 | 0.764 | 0.01 |
| linfonodo | M-3 | 0.672 | 0.765 | 0.757 | 0.08 |
| linfonodo | M-5 | 0.600 | 0.405 | 0.415 | 0.05 |
| linfonodo | M-6 | 0.472 | 0.339 | 0.339 | 0.00 |
| cervello | M-3 | 0.232 | 0.306 | 0.259 | 0.65 |
| cervello | M-5 | 0.904 | 0.942 | 0.938 | 0.09 |
| cervello | M-6 | 0.775 | 0.823 | 0.849 | -0.53 |
| cuore | M-3 | 0.164 | 0.199 | 0.181 | 0.50 |
| cuore | M-5 | 0.809 | 0.844 | 0.825 | 0.54 |
| cuore | M-6 | 0.750 | 0.708 | 0.698 | -0.24 |

- **CP-4 (il cuore come controllo negativo): fallita.** La coerenza spaziale (M-6) del cuore simulato è confermata (0.708 contro 0.685-0.816): la misura non vede le cellule allungate del cuore vero. Era stato messo in conto dopo la fase A (attesa indeterminata). Vuol dire che M-6 non serve a giudicare l'anisotropia; non tocca le altre conclusioni.

### Regola R-1

Giudizi smentiti fuori dal cuore: 17. Attribuzione (`F3.1_R1.csv`): disposizione dei tipi 7, non attribuiti 7, densità dei preset 2, profili 1. **Limite da sistemare per primo per la regola: la disposizione dei tipi** (7 giudizi, scarto medio 7.7 semi-ampiezze). Registrato come decisione tecnica in [F3.1: il limite del simulatore da sistemare per primo è la disposizione dei tipi cellulari (DEC-087)](../docs/decisions/DEC-087.md).

| misura | tessuto | scarto | limite attribuito |
|---|---|---|---|
| M-1 | intestino | +0.067 | densità dei preset |
| M-3 | intestino | +0.122 | non attribuito |
| M-4 | intestino | -0.031 | disposizione dei tipi |
| M-6 | intestino | -0.090 | disposizione dei tipi |
| M-1 | tumore del colon | +0.074 | densità dei preset |
| M-3 | tumore del colon | +0.106 | non attribuito |
| M-4 | tumore del colon | -0.139 | disposizione dei tipi |
| M-5 | tumore del colon | -0.068 | non attribuito |
| M-6 | tumore del colon | -0.173 | disposizione dei tipi |
| M-4 | stroma | -0.137 | disposizione dei tipi |
| M-5 | stroma | -0.086 | non attribuito |
| M-6 | stroma | -0.237 | disposizione dei tipi |
| M-3 | linfonodo | +0.093 | non attribuito |
| M-5 | linfonodo | -0.195 | non attribuito |
| M-6 | linfonodo | -0.133 | disposizione dei tipi |
| M-3 | cervello | +0.074 | profili dai nuclei con la corona |
| M-6 | cervello | +0.049 | non attribuito |

### Verità del simulato (descrittiva)

| tessuto | purezza delle cellule di bin2cell | recupero | quota da cellule di altro tipo | accuratezza del primo tipo di RCTD |
|---|---|---|---|---|
| intestino | 0.54 | 0.38 | 0.19 | 0.89 |
| tumore del colon | 0.74 | 0.60 | 0.08 | 0.95 |
| stroma | 0.85 | 0.36 | 0.10 | 0.81 |
| linfonodo | 0.59 | 0.59 | 0.23 | 0.45 |
| cervello | 0.86 | 0.29 | 0.04 | 0.86 |
| cuore | 0.85 | 0.13 | 0.06 | 0.81 |

