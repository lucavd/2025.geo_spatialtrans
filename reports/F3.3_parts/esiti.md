Tabelle in `results/F3.3/`: `F3.3_summary.csv` (previsioni e controprove), `F3.3_cells.csv` (una riga per tessuto × misura × insieme), `F3.3_area_values.csv` (valori per area), `F3.3_varcomp.csv`, `F3.3_r4a.csv`, `F3.3_bridge_cg.csv`, `F3.3_cp2.csv`, `F3.3_sens_insd.csv`, `F3.3_bl264_a4.csv`, `F3.3_D1.csv`, `F3.3_D2.csv`. Figure in `reports/F3.3_parts/`. Calcolo: catena delle aree di verifica dal 2026-11-19 18:48 alle 20:50 (4 processi), riproducibilità circa 1 ora (2 processi).

### Cancelli

| id | esito | dettaglio |
|---|---|---|
| G-1 | passato dopo una correzione decisa da Luca | (a)-(e) 205 controlli su 205 (`F3.3_G1_a_b_c_d_e_f.csv`): ritagli, Cellpose, Space Ranger, StarDist, consenso, nuclei e flag keep, maschere valide, regioni, esportazioni identici a R2, R3b, S1.6 e F3.1 sulle 30 aree di taratura; misure uguali a F3.1 su 90 esecuzioni (scarto massimo 1,3e-15). (f) **fallita la prova di S\* con uno spostamento di 1 h** (0,557 invece di 1): correzione [DEC-098](/decisions/DEC-098.md); le prove degli stimatori nuovi passano (`F3.3_G1_f.csv`: quadrato medio senza spostamento 0,007 e −0,007 h², spostamento di 1 h stimato 0,98 e 1,01) |
| G-2 | passato | 1 669 file letti con md5 identici prima e dopo (`/mnt/micron/geo_spatialtrans/F3.3/md5/`); dei file di verifica letti, 2 con md5 pubblicato uguali al manifesto di F3.2 e 2 con il solo ETag (già controllato in F3.2); controllo di fuga 0 occorrenze con 15 lettori di F3.3 ammessi; riga del registro delle letture scritta prima (commit della pre-registrazione `eedd71a`) |
| G-3 | passato dopo una correzione dello strumento | prima esecuzione 8 su 9 identiche; il linfonodo (aree di 500 µm, StarDist senza tasselli) differiva per 22 pixel su 895 933 (M-5 0,5327 contro 0,5304). Con StarDist deterministico (`F33_SD_DET=1`) le 2 aree di verifica del linfonodo sono state rifatte da StarDist in poi: 9 su 9 identiche (`F3.3_G3.json`) |
| G-4 | passato | 54 esecuzioni su 54 (`F3.3_G4.csv`) |
| G-5 | passato | due revisori indipendenti: misure 160 controlli su 162, statistiche 2 817; le discrepanze sono spiegate sotto («Revisione avversariale») |

### Previsioni

| id | quantità | atteso | osservato (intervallo) | esito |
|---|---|---|---|---|
| P-1 | spostamento del vero, aree nuove della stessa sezione (30 celle) | circa 0 | −0,20 h (−0,81; 1,12) | non rilevato a questa profondità |
| P-2 | spostamento del vero, altre sezioni o pazienti (18 celle) | sopra 0,5 | 2,38 h (2,07; 2,73) | smentita (spostamento grande), come atteso |
| P-3 | degrado del simulato, L1 (18 celle non smentite dentro campione) | circa 0 | −0,19 h (−1,14; 1,21) | non rilevato a questa profondità |
| P-3 | degrado del simulato, L2 (8 celle) | sopra 0,5 | 2,94 h (2,33; 3,43) | smentita (degrado grande), come atteso |
| P-4 | persistenza dei giudizi smentiti, L1 | ≥ 0,8 | 10 su 12 (0,83) | confermata |
| P-4 | persistenza, L2 | ≥ 0,8 | 19 su 26 (0,73) | smentita |
| P-5 | condizione 4 di R-3, L1 | rispettata | 0 misure confermate diventano smentite | rispettata |
| P-5 | condizione 4 di R-3, L2 | violata | UMI per cellula nello stroma smentiti nel paziente 1 e nel tumore da 6,5 mm | violata, come atteso |
| P-6 | ρ fra pazienti e fra aree (P1, P2, P5) | indeterminato | 12 righe su 12 indeterminate; con il tumore da 6,5 mm, 4 su 12 sopra 0,5 (UMI per cellula e singoletti nel tumore, singoletti e vicini dello stesso tipo nello stroma) | come atteso con 3 pazienti |
| P-7 | aree vere di verifica dentro il 90 % centrale del simulato | ≤ 0,10 | 18 su 228 (0,079) | come atteso |

Statistiche scritte nella pre-registrazione, giudicate com'erano scritte ([DEC-098](/decisions/DEC-098.md), punto 3): S\* L1 0,04 (−0,09; 0,16) confermata; G\* L1 −0,01 (−0,23; 0,23) confermata; S\* L2 1,46 (1,29; 1,59) smentita; G\* L2 1,84 (1,39; 2,31) smentita. Gli esiti di L2 coincidono; in L1 le statistiche scritte «confermano» dove quelle corrette dicono «non rilevato»: è il difetto per cui sono state sostituite.

**Ponte dei geni comuni** ([F3.3: scelte tecniche della misura fuori campione (DEC-097)](/decisions/DEC-097.md), punto 4): la restrizione ai 18 024 geni comuni alle quattro sezioni del tumore e al simulato (non 18 107 come in F3.2: là si contavano le matrici grezze) sposta le misure di P2 al più di 0,074 h (UMI per cellula nel tumore): il simulato resta quello di F3.1.

**Scomposizione descrittiva della tabella per cella** (non è un'indagine nuova: legge `F3.3_cells.csv`). In L2 il degrado sta nelle celle che dentro campione andavano bene: UMI assegnati (3,2 h) e UMI per cellula (2,6 h) nel tumore e nello stroma, concordanza con il nucleo nel cervello FFPE (4,4 h). I trascritti per cellula del tumore valgono 1 120 in P2, 640 in P1, 2 000 in P5, 270 nel tumore da 6,5 mm (simulato 1 096); nello stroma 663, 207, 567, 215 (simulato 653); nel cervello 4 514 nella sezione di taratura, 1 417 nella sezione FFPE. Fra i 7 giudizi smentiti che non restano tali in L2, 6 sono nello stroma degli altri pazienti (marcatori estranei, concordanza, vicini dello stesso tipo): lì il vero è più mescolato che in P2 e più vicino al simulato (vicini dello stesso tipo: P2 0,886, P1 0,752, P5 0,716, simulato 0,649; tumore da 6,5 mm 0,885).

### Controprove

- **CP-1, cuore: passata.** Anisotropia del campo degli UMI sulle 5 aree nuove, vero meno simulato 0,112 (errore standard 0,018), oltre 0,032 e oltre 2 errori standard (dentro campione 0,097). Negli altri tessuti, aree nuove: intestino 0,037, tumore 0,051, cervello −0,019 (il linfonodo non ha il valore simulato per i quarti).
- **CP-2, falso fuori campione: passata.** Media di S\* sui 10 modi di dividere le aree di taratura in 3 + 2, sui cinque tessuti con 5 aree: 0,007 (per tessuto da −0,032 a 0,045). Interpretazione dichiarata prima del calcolo: criterio sulla media dei 10 modi, non su ognuno.

### Documentazione

- **D-1** (`F3.3_D1.csv`): per tessuto da 40 a 80 valori scalari e 5 tabelle del simulatore vengono dalle 5 aree di taratura, cioè da 5 unità indipendenti di una sola sezione.
- **D-2** (`F3.3_D2.csv`): dei 14 giudizi confermati di F3.1, 6 sono UMI per cellula, costruiti dal fattore di profondità tarato in S2.5 sull'UMI per quadratino delle stesse aree; 4 sono UMI assegnati, «in parte» (densità e taglia dei nuclei dalle stesse aree); 4 no. Fuori campione: gli UMI per cellula restano senza smentita in L1 (spostamenti da −0,9 a +0,3 h) e sono smentiti in 5 insiemi L2 su 7.
