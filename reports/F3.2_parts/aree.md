Il confronto a vista di Luca è stato ritirato ([DEC-095](../docs/decisions/DEC-095.md)). Due regole automatiche hanno fallito la prova sulle aree di R2 (12 e 18 su 30); la terza (26 su 30, non cieca) funziona sulle sezioni di taratura ma non si trasferisce alle altre (kappa con i lettori 0,78 contro 0,26). Le aree tenute sono quelle su cui concordano regola e due lettori a vista (l'agente, registrato prima della terza regola, e un sotto-agente alla cieca) sulle sezioni di taratura e sui tessuti nuovi, e i due lettori sulle altre sezioni dello stesso tessuto. Accordo fra i lettori 0,907, kappa 0,882 su 107 candidate. Tabella completa: `results/F3.2/F3.2_area_agreement.csv`; coordinate: `results/F3.2/F3.2_selected_*.csv`.

| livello | sezione | classe | aree |
|---|---|---|---|
| L1 | intestino di topo | mucosa | 4 |
| L1 | tumore del colon (P2) | tumore | 3 |
| L1 | tumore del colon (P2) | stroma | 0 |
| L1 | linfonodo (500 µm) | follicolo | 2 |
| L1 | linfonodo (500 µm) | paracorticale | 0 |
| L1 | cervello di topo | corteccia | 1 |
| L1 | cuore | miocardio | 5 |
| L2 | tumore P1 | tumore, stroma | 4, 3 |
| L2 | tumore P5 | tumore, stroma | 3, 1 |
| L2 | tumore 6,5 mm (paziente incerto) | tumore, stroma | 5, 3 |
| L2 | cervello di topo, seconda sezione | corteccia | 4 |
| L3 | rene umano | tessuto | 5 |
| L3 | fegato di ratto | tessuto | 5 |

Prova del lato dell'area (solo dati di taratura, `results/F3.2/F3.2_size_check.csv`): 20 casi su 24 entro un quarto della variabilità fra aree (19 su 23 informativi), tutti e 4 nel linfonodo. Paziente del tumore da 6,5 mm (`results/F3.2/F3.2_sexcheck.csv`): XIST assente dal pannello, geni del cromosoma Y 1,1, 0, 0 e 3,7 UMI per 10 000 in P2, P1, P5 e 6,5 mm: non decidibile.
