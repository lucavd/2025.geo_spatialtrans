Commit dei cancelli: `c4832c2`. Ok di Luca sulle previsioni: 2026-11-08, «Ok, come proposte (Recommended)» (chat), soglie comprese (criterio di Bradley per le coperture; metà della semi-ampiezza corretta per P-3). Commit: questo commit, «R7: ok di Luca sulle previsioni (2026-11-08)».

### Metodi

- **M-1, geometria.** Tre matrici per tessuto: quota di area in comune fra i quadrati (F_box); quota dei bin di 2 µm in tessuto condivisi (F_bin, dai file per area di R4, posizione globale = posizione nell'area + angolo dell'area, in pixel interi); quota dei nuclei condivisi (F_nuc, nuclei Cellpose-HED di R2 tenuti, scala 1, centroide globale). Correlazione per coppia n_ij / √(n_i n_j). F_box è la matrice primaria; F_bin e F_nuc servono a controllarla (regola R-2).
- **M-2, pezzi disgiunti** (l'alternativa indicata nel [step sulle aree sovrapposte prima di S2.5 (DEC-048)](../docs/decisions/DEC-048.md)). Ogni punto dell'unione delle 5 aree va all'area che lo contiene con il centro più vicino (a parità, l'area con l'indice più basso). Ogni nucleo e ogni bin conta una volta. Nel cuore i pezzi coincidono con le aree.
- **M-3, simulazione di copertura.** Valori per area x = μ + σ·z, z normale con matrice di correlazione data, 5 aree (4 per le finestre, indipendenti perché disgiunte). Metodo di oggi: bootstrap percentile della media con enumerazione esatta delle 5⁵ estrazioni (B infinito; equivalenza con `R4_common.boot_one_stage` a B = 2000 nel cancello G-1). Candidati: **C-a**, allargamento calibrato: le due semi-ampiezze dell'intervallo registrato (stima − inferiore, superiore − stima) moltiplicate per κ, il fattore che porta la copertura a 0.95 con la correlazione F_box del tessuto; **C-b**, intervallo t con l'effetto della correlazione: x̄ ± t(n−1; 0.975) · s/√n · k, con k = √((1 + (n−1)ρ̄)/(1 − ρ̄)), ρ̄ media delle correlazioni fuori diagonale. Per gli intervalli minimo-massimo (descrittivi: «fra aree vere») si corregge solo la zona comune: κ_range = rapporto fra l'ampiezza attesa di 5 aree indipendenti e quella con F_box. Repliche 20 000 per ogni caso (errore Monte Carlo circa 0.0015 a 0.95); seme 20261108.
- **M-4, correlazione dentro le aree** (base di CP-1). Campi: UMI medi per bin di 2 µm in tessuto e nuclei per mm² di tessuto, su tessere di 48 µm di una griglia globale, solo nell'unione delle 5 aree (tessere con almeno metà dei bin in tessuto; nuclei presi dal pezzo disgiunto di ciascuna area, così nessuno conta due volte). Semivariogramma empirico fra 48 e 720 µm in 15 classi, coppie di tessere della stessa area; modello γ(h) = c0 + c1 (1 − exp(−h/L)) adattato per minimi quadrati pesati (pesi N_h/γ(h)²), L fra 10 e 2 000 µm. Correlazione implicita fra le medie di due aree: media delle covarianze del modello fra le tessere delle due aree, normalizzata. Incertezza: bootstrap sulle 5 aree (si ricampionano i semivariogrammi per area e si riadatta), B = 2000, seme 20261108, 4 000 repliche di copertura per bootstrap.
- **M-5, riesame dei verdetti passati.** Per ogni verdetto degli step da R2 a S2.4 nelle tabelle dei verdetti o delle previsioni, se la stima o l'intervallo di rilevanza viene da aree o finestre: stesso giudizio con l'intervallo corretto dalla regola R-1 (κ sul lato per area, κ_range sui margini minimo-massimo), con la funzione dei verdetti di R4 (`R4_common.verdict`). Tabella delle transizioni per step e tessuto. Ogni transizione che tocca una scelta oggi in uso (preset, modello, riferimento adottato) è un dubbio da chiudere.

### Regole di decisione dei dati (scritte prima del calcolo)

- **R-1, la correzione** (BL-186, decisa dai dati per scelta di Luca). Si adotta, fra C-a e C-b, quella con la copertura dentro l'intervallo di Bradley (0.925-0.975) nel maggior numero di celle tessuto × scenario, sugli scenari di P-1 (normale, lognormale, t con 3 gradi di libertà); a parità, C-a (si applica agli intervalli registrati senza i valori per area e ne conserva l'asimmetria). Se P-2 è smentita in un tessuto, in quel tessuto κ si calibra sulla correlazione implicita di M-4 (stima puntuale) invece che su F_box.
- **R-2, quale matrice di zone comuni.** Se F_bin e F_nuc differiscono da F_box di più di 0.05 in una coppia, per le quantità per bin e per nucleo si usa la matrice specifica; altrimenti F_box per tutte.
- **R-3, due stadi.** Per le righe con bootstrap a due stadi (aree, poi cellule) si usa il κ a uno stadio; il rapporto fra le semi-ampiezze a due e a uno stadio si misura su una riga di R4 per tessuto (misura M-f).

### Cancelli (al massimo 5)

| id | controllo | tolleranza | comando |
|---|---|---|---|
| G-1 | batteria `tools/R7_test.py`: quote di zona comune su quadrati sintetici noti; formula di k contro simulazione bruta di Var(x̄) ed E[s²]; intervallo t con aree normali indipendenti copre 0.950; enumerazione esatta uguale in distribuzione a `R4_common.boot_one_stage` (B = 200 000); pezzi disgiunti: ogni punto in un solo pezzo, unione dei pezzi = unione delle aree; funzione dei verdetti uguale a `R4_common.verdict`. Almeno 6 mutanti (F ignorata, gradi di libertà n invece di n−1, centro più lontano, κ su un solo lato, bootstrap senza reinserimento, κ_range invertito), tutti rilevati | numerica (copertura entro 3 errori Monte Carlo) | `.venv/bin/python tools/R7_test.py` |
| G-2 | geometria e pezzi sui dati veri: F_box uguale ai valori di BL-215 ed E-0 entro 0.005; cuore senza zone comuni; ogni nucleo e ogni bin in un solo pezzo; nel cuore valori sui pezzi identici bit per bit a quelli sulle aree | 0.005; identità | `tools/R7_run.sh G2` |
| G-3 | riproduzione: le quantità per area di R4 ricalcolate dai file per area senza pezzi uguali a `results/R4/R4_roi.csv`; gli intervalli registrati in `results/R4/R4_refs.csv` (bootstrap a uno stadio) ricalcolati dai valori per area con i semi di R4 | 1e-9 relativo | `tools/R7_run.sh G3` |
| G-4 | integrità: md5 prima e dopo di tutto ciò che si legge (aree, nuclei, file per area di R4, maschere, tabelle dei verdetti, `BIO_REFERENCES.md`) | identità | `tools/R7_md5.sh before/after` |
| G-5 | revisione avversariale: sotto-agenti indipendenti ricalcolano dal codice e dai dati κ e κ_range per tessuto, le coperture di P-1 e P-2, gli scarti di P-3 e le transizioni di M-5 | numerica | report, sezione revisione |

### Previsioni

Unità indipendente per l'incertezza: replica della simulazione (P-1); area di R2, con bootstrap sulle 5 aree (P-2); nessuna per P-3 (scarto calcolato esattamente sugli stessi dati). B = 2000; seme = 20261108. Intervallo di rilevanza per le coperture: criterio liberale di Bradley, 0.925-0.975, scelto da Luca il 2026-11-08 (regola 3 del protocollo: non è un valore biologico; è la convenzione già usata in MMM).

| id | quantità | archetipo | atteso | intervallo di rilevanza | fonte dell'intervallo | 2 errori standard tipici |
|---|---|---|---|---|---|---|
| P-1 | copertura dell'intervallo corretto con R-1 quando i valori per area non sono normali: lognormali con CV alto (90° percentile dei CV fra aree delle quantità di `R4_roi.csv`) e t con 3 gradi di libertà (marginali, correlazione F_box) | tutti, A1-A6 | dentro l'intervallo in 6 tessuti su 6 per entrambi gli scenari | 0.925-0.975 | Bradley | 0.003 (Monte Carlo) |
| P-2 (controprova CP-1) | copertura dell'intervallo corretto con κ calibrato su F_box quando le aree hanno la correlazione implicita nella struttura del tessuto misurata dentro le aree (M-4), per UMI per bin e per nuclei | A1-A5 (nel cuore nessuna zona comune: la correzione per le zone comuni vale 1) | dentro l'intervallo, per entrambi i campi, in 5 tessuti su 5 | 0.925-0.975 | Bradley | non noto prima del calcolo; se l'intervallo bootstrap è più largo di metà dell'intervallo di rilevanza, l'esito sarà indeterminato |
| P-3 (controprova CP-2) | scarto del valore del tessuto quando ogni nucleo e ogni bin conta una volta (pezzi disgiunti) invece che nelle aree: UMI per bin a 2 e 8 µm, frazione di bin vuoti a 2 µm, geni per bin a 8 µm, rapporto nucleo/corona non corretto, area nucleare mediana, densità dei nuclei | A1-A5 (nel cuore è 0 per costruzione, cancello G-2) | dentro l'intervallo per tutte le 7 quantità in 5 tessuti su 5 | ± metà della semi-ampiezza dell'intervallo corretto della stessa quantità (circa 1 errore standard) | variabilità fra aree misurata sui dati veri (regola 2 del protocollo), fattore 0.5 scelto da Luca il 2026-11-08 | nessuno: lo scarto è esatto per questi dati |

### Controprova

- **CP-1 = P-2.** Se la struttura del tessuto dentro le aree rendesse le aree sovrapposte molto più correlate della sola quota in comune, la correzione calibrata sulla geometria sarebbe troppo piccola e la copertura vera scenderebbe sotto 0.925. Allora R-1 calibra κ su quella correlazione.
- **CP-2 = P-3.** La correzione cambia solo gli intervalli. Se contare due volte le zone comuni spostasse il valore stesso dei riferimenti, andrebbero rifatti i valori, non solo allargati gli intervalli.
- **CP-3, cuore (controllo negativo per le zone comuni).** Senza zone comuni la parte di correzione dovuta alle zone comuni vale esattamente 1, e i pezzi coincidono con le aree (cancello G-2).

### Misure (senza intervallo; servono alle decisioni)

- M-a: κ, κ_range e k per tessuto e per numero di unità (5 aree, 4 finestre); copertura del metodo di oggi a 20 000 repliche.
- M-b: transizioni dei verdetti passati (M-5), per step, tessuto e tipo di transizione.
- M-c: riferimenti corretti per il confronto finale (righe attive con unità area o finestra): valore, intervallo registrato, intervallo corretto.
- M-d: quota dei nuclei dell'area tenuta fuori che stanno anche nelle aree di stima, per le validazioni lasciando fuori un'area (S1.4, S1.5c, S1.5d, S2.2).
- M-e: deviazione standard fra pezzi disgiunti contro quella fra aree, per le 7 quantità di P-3.
- M-f: rapporto fra le semi-ampiezze del bootstrap a due e a uno stadio, una riga di R4 per tessuto.
