
| cosa | perché | cosa serve per chiuderlo | step destinatario | voce di BACKLOG |
|---|---|---|---|---|
| triangolazione fallita nel cervello r3 con la semina a compartimenti | difetto del codice dello Step 2, fuori dal perimetro | riprodurre con `F31_LAYER=Q6 Rscript --vanilla tools/F3.1_sim.R <dir> 20261114 A5_r3 1 0` e correggere `.s27_cdt()` | step della semina a compartimenti predefinita (o primo step sul codice dello Step 2) | BL-259 |
| linfonodo: plasmacellule in eccesso e concordanza bassa con RCTD, non spiegate dalla mescolanza | serve guardare i profili dei linfociti B e delle plasmacellule e la composizione per area (follicoli e paracorteccia) | rifare le misure dopo la semina a compartimenti; se resta, confrontare i profili simulati dei B con il riferimento | step della semina a compartimenti predefinita | BL-260 |
| cervello: RCTD chiama neuroni il 75 % dei quadrati simulati contro il 55 % dei veri; coerenza spaziale in eccesso | non dipende dalla mescolanza né dai profili | confrontare la composizione del preset del cervello con i tipi chiamati sul vero, area per area | da assegnare (con la densità dei preset) | BL-261 |
| le aree simulate di un tessuto variano fra loro molto meno delle vere | i metodi vedrebbero un tessuto più omogeneo del vero | misurare la variabilità fra aree delle quantità principali e decidere se il simulatore deve produrla | da assegnare | BL-262 |
| la coerenza spaziale dei tipi non vede le cellule allungate del cuore (CP-4) | serve una misura orientata per giudicare il controllo negativo | definire una misura dell'allungamento delle zone dello stesso tipo | da assegnare (lavoro sul cuore) | BL-263 |
| definizioni delle misure per cellula (quadratini fuori dal tessuto; marcatori comuni a due tipi nel linfonodo) | segnalate dalla revisione; non spostano esiti | precisarle prima di riusare le misure | prossimo step che rifà le misure dei metodi | BL-264 |

