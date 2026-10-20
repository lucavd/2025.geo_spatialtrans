---
type: Step
title: "Nel linfonodo e nel cuore, da quale insieme di nuclei veri prendiamo posizioni, aree e vicini?"
description: "Sceglie per linfonodo e cuore una sola maschera nucleare (i nuclei di Space Ranger con la forma di StarDist, oppure solo StarDist), confrontandole con il conteggio manuale dei nuclei e dei contatti di Luca, e controlla con i punti manuali se l'eccesso di contatti simulati nel linfonodo dipende davvero dai centri dei nuclei."
step_id: "R3b"
status: draft
branch: "step1-cell-layer"
resource: "reports/R3b.html"
depends_on: ["S1.4d (chiusa, protocollo v2)", "R2c (chiusa con riserve, protocollo v1)", "R2b (chiusa con riserve, protocollo v1)"]
pi_approval: pending
generated: { by: geo-agent/skill-v2, at: 2026-10-20T14:10:00+02:00 }
sources:
  - id: roadmap-s7
    resource: "git:step1-cell-layer@b5aaa76:docs/ROADMAP.md (sezione 7, prossima sessione; sezione 6, righe 2026-10-07 R2b, 2026-10-08 BL-060, 2026-10-20 R2c)"
    title: "Perimetro di R3b e decisioni passate sul segmentatore di linfonodo e cuore"
  - id: backlog
    resource: "git:step1-cell-layer@b5aaa76:docs/BACKLOG.md (BL-064, BL-100, BL-115, BL-121, BL-130)"
    title: "Voci di backlog assegnate a R3b"
  - id: r2b-match
    resource: "git:step1-cell-layer@b5aaa76:results/R2b/R2b_matching.csv (righe rater Luca, A4 e A6)"
    title: "Oggetti dei segmentatori contro i punti manuali di Luca nelle 4 finestre per tessuto (già visto, non cieco)"
  - id: b046
    resource: "BIO_REFERENCES B-046 (linfonodo) e B-048 (cuore): densità manuale, IC Poisson 95 % del totale"
    title: "Densità nucleare contata a mano, con il suo intervallo"
  - id: b130
    resource: "BIO_REFERENCES B-130 (linfonodo) e B-132 (cuore); results/R2c/R2c_posthoc_ci.csv, righe «due stadi»"
    title: "Frazione di nuclei a contatto contata a mano da Luca, con intervallo bootstrap a due stadi"
  - id: r2c-seg
    resource: "git:step1-cell-layer@b5aaa76:results/R2c/R2c_s1.csv (seg_spaceranger, seg_stardist_he) e R2c_seg_focal.csv"
    title: "Contatti nelle maschere dei segmentatori sui nuclei focali di R2c (già visto, non cieco)"
  - id: s14d
    resource: "git:step1-cell-layer@b5aaa76:results/S1.4d/S1.4d_windows.csv (K2_tissue, K2_stardist) e docs/steps/S1.4d.md"
    title: "Contatti del nucleo ritagliato e gonfiato sulle finestre, con centri Space Ranger e StarDist (già visto, non cieco)"
  - id: protocol
    resource: "/decisions/DEC-002.md"
    title: "Protocollo di sessione versione 2"
---
# Dubbi ricevuti

Le voci di BACKLOG che la riga «Prossima sessione» di ROADMAP assegna a R3b[^roadmap-s7][^backlog]. Si leggono per prime.

| voce | cosa | cosa serve per chiuderla in R3b |
|---|---|---|
| BL-121 | nel linfonodo e nel cuore i nuclei simulati partono dai centri di Space Ranger, ma aree e vicini reali vengono da StarDist: tre grandezze da due maschere diverse | scegliere una sola maschera per posizioni, aree e vicini (regola di decisione sotto) e ricalcolare con quella i riferimenti di questi due tessuti |
| BL-130 | il nucleo ritagliato e gonfiato dà troppi contatti nel linfonodo con i centri Space Ranger (0.74 contro il manuale 0.55); con i centri StarDist 0.41. S1.4d ha attribuito l'eccesso ai centri troppo densi | controprova CP-1 con i punti manuali di Luca come centri: dice se l'eccesso viene dai centri o dal modello del nucleo |
| BL-100 | i parametri del modello di taglia di linfonodo e cuore sono stimati sui soli nuclei Space Ranger appaiati a un nucleo StarDist (71-76 % nel linfonodo, 47-77 % nel cuore), e gli esclusi hanno territori più piccoli | con una maschera unica ogni nucleo ha la sua area: il modello si stima su tutte le celle interne; si riporta di quanto cambiano i parametri |
| BL-115 | la misura dei contatti reali nel linfonodo e nel cuore è dominata dai vicini non appaiati | con una maschera unica non esistono vicini non appaiati; si riporta la frazione di contatti nella maschera scelta |
| BL-064 | forma vera delle cellule del cuore (i territori di Voronoi del cuore sono rotondi quanto quelli casuali) | **non in R3b**: rinviata a S1.5 per scelta di Luca (2026-10-20, in chat: «Rinviarla a S1.5»); voce aggiornata alla chiusura |

**Un fatto già nei file cambia la lettura di BL-130**[^r2b-match]. Nelle 4 finestre di linfonodo Space Ranger conta 417 oggetti contro i 404 punti di Luca (rapporto medio per finestra 1.04), StarDist 314 (0.77). Nel cuore 264 e 202 contro 242 (1.10 e 0.84). I centri di Space Ranger quindi non sono troppi rispetto al conteggio manuale; StarDist ne vede troppo pochi. Il valore di contatti di StarDist nel margine del manuale (0.41) potrebbe venire dai nuclei mancanti, cioè essere giusto per la ragione sbagliata. La controprova CP-1 serve a distinguere i due casi.

# Attese pre-registrate

Regole del protocollo v2[^protocol]. Commit dei cancelli: vedi `git log` (messaggio «R3b: pre-registrazione»). Ok di Luca sulle previsioni: 2026-10-20, in chat, alla domanda «Approvi le previsioni, i loro intervalli e la regola di decisione così come sono scritti, compresi gli esiti indeterminati dichiarati in anticipo (conteggio di «Space Ranger con la forma di StarDist» al bordo superiore; contatti nella maschera con poco potere su 4 finestre)?»: «Approvo così». Commit: vedi `git log` (messaggio «R3b: ok di Luca sulle previsioni»).

## Le maschere a confronto

Per ognuno dei due tessuti (linfonodo, 5 ROI di 1 mm², 4 con finestre manuali; cuore, idem):

- **CO, «i nuclei di Space Ranger con la forma di StarDist»** (candidata). Space Ranger decide quali nuclei esistono; dove un nucleo Space Ranger ha un corrispondente StarDist, la forma (pixel) è quella di StarDist; dove non lo ha, resta il poligono Space Ranger. I nuclei visti solo da StarDist non entrano. È la versione coerente della scelta del 2026-10-08 (posizioni Space Ranger, forma StarDist), con aree e vicini presi dalla stessa maschera.
- **SD, «solo StarDist»** (candidata). La maschera StarDist `2D_versatile_he` già usata come sensibilità in R3, S1.4 e S1.4d.
- **SR, «solo Space Ranger»** (riferimento, non candidata). Le sue forme sono state escluse da Luca il 2026-10-08 (quantizzate sulla griglia di 2 µm, maschere contigue). Si riportano P-1 e P-2 per mostrare l'effetto della contiguità; il nucleo simulato non si rigenera.
- **Stato attuale** (riferimento): centri Space Ranger, aree StarDist sui soli appaiati, cioè S1.4 e S1.4d. P-3 = 0.744 nel linfonodo, 0.049 nel cuore[^s14d].

Cellpose non è candidata: nel cuore trova meno di metà dei nuclei di Luca (rapporto 0.48)[^r2b-match], e la decisione di R2b ha già scelto Space Ranger e StarDist per questi due tessuti.

## Cancelli (al massimo 5)

| id | controllo | tolleranza | comando |
|---|---|---|---|
| G-1 | **batteria di test** del costruttore della maschera CO su etichette sintetiche (appaiamento uno a uno; due StarDist dentro un Space Ranger; uno StarDist su due Space Ranger; Space Ranger coperto del tutto da StarDist e quindi vuoto; Space Ranger senza corrispondente; identificativi unici; nessun pixel con due etichette; pixel StarDist appaiati identici all'ingresso; filtro del tessuto valido) e della tabella dei nuclei che ne deriva (centroide, area, eccentricità, orientazione contro `skimage.measure.regionprops`) | tutte PASS | `bash R/testing/run_R3b.sh test` |
| G-2 | **mutanti** del costruttore, ognuno deve far fallire almeno un'asserzione di G-1: M1 includere i nuclei solo StarDist; M2 dipingere Space Ranger sopra StarDist; M3 appaiamento non reciproco; M4 non togliere le etichette rimaste vuote; M5 saltare il filtro del tessuto valido | 5/5 rilevati | `bash R/testing/run_R3b.sh mutants` |
| G-3 | **regressione della catena** (maschera, tabella dei nuclei, generatori, territori, modello di taglia, nucleo ritagliato, raster, finestre; e il codice di P-1 e P-2): con la maschera StarDist in ingresso riproduce i file R3 `*_stardist_he_gen.parquet`, i territori e i parametri di S1.4 (`S1.4/real`, `S1.4/loo`, StarDist) e le righe `K2_stardist` di `S1.4d_windows.csv` (con φ = 0.816 come in S1.4d); con le maschere SR e SD riproduce `n_obj` e `tp` di `R2b_matching.csv` e `seg_contact` di `R2c_seg_focal.csv` | identico (scarto 0 o ≤ 1e-12 sui numeri in virgola mobile) | `bash tools/run_R3b.sh regress` |
| G-4 | **invarianti della maschera CO e integrità** sui 10 ROI: nessun pixel con due etichette; ogni nucleo StarDist appaiato identico pixel per pixel; nessun pixel Space Ranger dove c'è un nucleo StarDist appaiato; numero di etichette = appaiati + solo Space Ranger non vuoti; ogni centroide nel tessuto valido; rerun di 2 ROI identico bit per bit (md5); md5 dei file letti (maschere R2, annotazioni R2b e R2c, uscite R3, S1.4, S1.4c, S1.4d) uguali prima e dopo | 0 violazioni; md5 identici | `bash tools/run_R3b.sh check` |
| G-5 | **ricalcolo indipendente** della revisione avversariale: due revisori ricalcolano dal codice e dai dati la maschera CO di almeno 2 ROI, P-1…P-4, CP-1 e i conteggi di G-4 | scarto ≤ 1e-9 o identico | revisione |

## Previsioni

Unità indipendente (dichiarata ora): **finestra** per P-1 (4 per tessuto; stima = media dei 4 rapporti, errore standard = deviazione standard fra finestre / √4); **finestra e nucleo focale** per P-2 (bootstrap a due stadi come R2c: finestre, poi nuclei focali dentro la finestra; B = 2000; seme 20261021; intervallo percentile 95 %); **seme** per P-3 e CP-1 (10 semi, come S1.4d: finestre e centri fissi, l'incertezza del riferimento manuale è già nel suo intervallo); **ROI** per P-4 (5 per tessuto). Semi del nucleo simulato: `seed_fold + 100000 · s`, s = 1…10, con `seed_fold = 20261015 + 1000 · i_tessuto + 10 · k` come S1.4 e S1.4d, quindi appaiati a S1.4d. Estremi tagliati a [0, 1] per le frazioni.

**Non tutto è cieco.** I valori di SR e SD per P-1 e P-2 sono già nei file di R2b e R2c[^r2b-match][^r2c-seg]; P-3 di SD è noto con il φ di Space Ranger (0.41 linfonodo, 0.012 cuore)[^s14d] e qui cambia solo φ. Sono ciechi: tutto ciò che riguarda CO, P-3 di SD con il suo φ, P-4, CP-1.

| id | quantità | tessuto | maschera | atteso | intervallo di rilevanza | fonte dell'intervallo | 2 errori standard tipici |
|---|---|---|---|---|---|---|---|
| P-1 | rapporto fra nuclei della maschera e punti di Luca nelle stesse finestre (regole di R2b: oggetti con centroide nella finestra e fuori dalle esclusioni) | linfonodo | CO | 1.00-1.04, **indeterminata** (estremo superiore oltre 1.102) | 0.905-1.102 | B-046: IC Poisson 95 % del totale manuale / totale[^b046] | 0.06 |
| P-1 | idem | linfonodo | SD | 0.77, smentita (nota) | 0.905-1.102 | B-046 | 0.06 |
| P-1 | idem | cuore | CO | 1.05-1.09, **indeterminata** | 0.878-1.134 | B-048 | 0.06 |
| P-1 | idem | cuore | SD | 0.84, **indeterminata** (nota) | 0.878-1.134 | B-048 | 0.11 |
| P-2 | frazione dei nuclei focali di Luca (R2c) appaiati a un nucleo della maschera che toccano un altro nucleo della stessa maschera (pixel entro 2 px, operatore di S1.4c) | linfonodo | CO | 0.45-0.65, confermata o indeterminata | 0.35-0.75 | B-130, due stadi[^b130] | 0.15-0.20 |
| P-2 | idem | linfonodo | SD | 0.37, **indeterminata** (nota) | 0.35-0.75 | B-130 | 0.15-0.20 |
| P-2 | idem | linfonodo | SR (riferimento) | 0.98, smentita (nota: maschere contigue) | 0.35-0.75 | B-130 | 0.05 |
| P-2 | idem | cuore | CO, SD | 0.02-0.04, **indeterminata** (con 1-2 contatti su circa 50 focali l'estremo superiore supera 0.07) | 0.00-0.07 | B-132 | 0.05-0.08 |
| P-3 | frazione di nuclei simulati a contatto sulle finestre di R2c (regola dei focali di R2c, peso uguale per finestra, come P-1 di S1.4d), con centri, modello di taglia (M2 stimato a fold sulla maschera) e φ (q99 del rapporto area del nucleo / area del territorio delle celle interne della maschera) della stessa maschera | linfonodo | CO | 0.65-0.75, **indeterminata** o confermata | 0.35-0.75 | B-130 | 0.015-0.03 |
| P-3 | idem | linfonodo | SD | 0.38-0.45, confermata | 0.35-0.75 | B-130 | 0.03 |
| P-3 | idem | cuore | CO | 0.03-0.06, confermata | 0.00-0.07 | B-132 | 0.015 |
| P-3 | idem | cuore | SD | 0.005-0.02, confermata | 0.00-0.07 | B-132 | 0.015 |
| P-4 | frazione di nuclei simulati troncati dal tetto φ (area finale < area bersaglio), tutte le celle, media dei 10 semi per ROI | linfonodo, cuore | CO, SD | 0.01-0.04, confermata | 0-0.05 | ROADMAP §4, criterio di Luca, come P-2 di S1.4d | 0.01 |

**Esiti indeterminati dichiarati ora.** P-1 di CO nei due tessuti: 2 errori standard (0.06) sono circa metà dell'intervallo, e l'atteso sta vicino all'estremo superiore, come per Space Ranger. P-2 in entrambi i tessuti: l'intervallo bootstrap di 4 finestre e circa 50 nuclei focali è largo quanto l'intervallo di rilevanza, quindi P-2 può escludere solo una maschera molto sbagliata (come SR, contigua). Più finestre manuali stringerebbero entrambe le misure: è un'opzione per Luca, non parte di R3b.

Cosa direbbe una smentita.

- **P-1** smentita: la maschera ha troppi o troppo pochi nuclei rispetto al conteggio manuale; i territori, e quindi il rapporto nucleo/territorio e il tetto φ, sarebbero sbagliati in proporzione. Per SD nel linfonodo è attesa.
- **P-2** smentita: nella maschera i nuclei si toccano molto più o molto meno di quanto Luca veda al microscopio; i vicini reali di quella maschera non sono affidabili.
- **P-3** smentita: con quella maschera il nucleo simulato non riproduce i contatti manuali. Nel linfonodo, per CO, è possibile; CP-1 dice se la colpa è della maschera o del modello del nucleo.
- **P-4** smentita: il tetto φ stimato su quella maschera tronca più del 5 % dei nuclei; tetto e modello di taglia sarebbero incoerenti in quel tessuto.

Il cuore, controllo negativo previsto per il Voronoi isotropo, qui non deve fallire: le quantità di questo step (conteggio, contatti, troncamento) non dipendono dall'allungamento dei territori, e i contatti nel cuore sono rari.

## Controprova

**CP-1: i contatti simulati con i punti manuali di Luca come centri.** Nelle 4 finestre di ogni tessuto i centri sono i punti di Luca (R2b); fuori dalle finestre e dentro le sue esclusioni, i centri di CO. Modello di taglia e φ di CO (fold in cui il ROI è di test). Stesso nucleo ritagliato e gonfiato, stessi semi, stessa regola dei focali di P-3 (i nuclei simulati misurati sono quelli dei punti di Luca).

- Atteso nel linfonodo: 0.60-0.75, **indeterminata** o confermata, perché i centri di Space Ranger sono già quasi tanti quanti i punti di Luca. Nel cuore 0.02-0.05, confermata. Intervalli B-130 e B-132; 2 errori standard circa 0.02.
- **Se CP-1 è smentita nel linfonodo** (tutta sopra 0.75): anche con i centri veri il nucleo simulato tocca troppo i vicini. L'eccesso non viene dalla maschera e la spiegazione di S1.4d (centri troppo densi, BL-130) è falsa. P-3 esce dalla regola di decisione e si apre un dubbio sul modello del nucleo (taglia o forma), da chiudere o rinviare.
- **Se CP-1 è confermata e P-3 di CO è smentita**: i centri contano, e quelli di CO hanno un difetto (per esempio i falsi positivi di Space Ranger, 12 % nel linfonodo[^r2b-match]); dubbio aperto su CO.

## Regola di decisione (per tessuto, applicata dall'agente; Luca può porre il veto)

1. Una maschera candidata è **eleggibile** se né P-1 né P-2 sono smentite.
2. Se una sola è eleggibile, si sceglie quella.
3. Se lo sono entrambe: se CP-1 non è smentita, si preferisce quella con P-3 confermata; se P-3 è confermata per entrambe o per nessuna, oppure se CP-1 è smentita, si sceglie quella con il rapporto di P-1 più vicino a 1 (minimo |log rapporto|).
4. Se nessuna è eleggibile, la regola non sceglie: si portano a Luca le opzioni con le loro conseguenze.

Esito atteso: **CO in entrambi i tessuti** (nel linfonodo SD ha troppo pochi nuclei; nel cuore CO ha il rapporto più vicino a 1: |log 1.09| = 0.09 contro |log 0.84| = 0.17).

## Cosa si fa comunque con la maschera scelta (descrittivo, nessuna previsione)

Con la maschera scelta si ricalcolano per linfonodo e cuore, con il codice di R3, S1.4 e S1.4c: geometria dei territori (B-073…B-078, B-085…B-090), modello di taglia (B-100…B-102, B-106…B-108), intervalli di C8 (B-112, B-114), q99 del rapporto nucleo/territorio e contatti nella maschera (B-124, B-126). Le righe nuove sostituiscono le vecchie, che restano marcate come superate. Nei banchi la maschera primaria di questi due tessuti diventa quella scelta. Si riporta di quanto cambia ogni valore.

## Definizioni operative (scritte prima del calcolo)

- **Ingressi**: maschere native di R2 (`/mnt/micron/geo_spatialtrans/R2/masks/<A>_<roi>_{spaceranger,stardist_he}_native.npz`), flag `keep` di R2 (`results/R2/nuclei/*_native.parquet`), maschera del tessuto valido (filtro come lo stadio cp2 di S1.4c: centroide nel tessuto valido).
- **Appaiamento per CO**: tabella delle sovrapposizioni in pixel fra etichette Space Ranger e StarDist; la coppia (i, j) è appaiata se j è lo StarDist con la massima sovrapposizione con i, i è lo Space Ranger con la massima sovrapposizione con j, e la sovrapposizione è > 0 (pareggi: identificativo più piccolo). Si riportano le frazioni di sovrapposizione delle coppie.
- **Raster CO**: si dipingono prima gli Space Ranger senza corrispondente, poi gli StarDist appaiati sopra (vincono i pixel StarDist). Uno Space Ranger rimasto senza pixel si toglie (conteggio riportato). Nuove etichette 1…n in ordine di etichetta Space Ranger.
- **Tabella dei nuclei** di CO come R2 (`tools/R2_metrics.py`): centroide in µm, area in µm², eccentricità e orientazione da `regionprops` sui pixel; poi R3 (generatori), S1.4 (territori con `tessellate_voronoi()`, tutte le celle interne con la propria area, quindi `matched = TRUE` per tutte), S1.4 loo (M2, 5 fold), φ, S1.4d (nucleo ritagliato e gonfiato, raster con `tools/S1.4c_raster.py sim`, finestre con la regola di `tools/S1.4d_windows.py`).
- **P-1, P-2 su una maschera nuova**: codice di `tools/R2b_match.py` (appaiamento punto-oggetto: dentro la maschera, oppure centroide più vicino entro 3 µm, assegnazione 1:1) e di `tools/R2c_refs.py` (contatto = etichetta in una coppia a contatto di `contact_pairs(lab, 2)`).
- **CP-1**: punti di Luca in µm del ROI (`results/R2b/R2b_points.csv`, rater Luca); centri CO tolti dentro il quadrato della finestra salvo dentro le esclusioni di Luca, dove restano.
- Calcolo su lesexp-server, al massimo 6 worker a un thread (server condiviso con MMM), file `WORKERS` per cambiarli in corsa, `OMP_NUM_THREADS=1`; uscite stabili in `/mnt/micron/geo_spatialtrans/R3b/`, un file per unità, rilanciabili. Stima: meno di un'ora di calcolo (tempi di S1.4d: 50-130 s per fold di linfonodo, 3 s per il cuore).

## Note di sviluppo e definizioni operative aggiunte durante il calcolo (non cambiano previsioni, intervalli, regola)

- **D-1 md5 «prima».** L'istantanea md5 dei 343 file letti è stata presa alle 13:59 del 2026-10-20, dopo l'avvio della costruzione della maschera CO, che legge soltanto. Controprova: 120 di quei file sono anche nell'elenco md5 di chiusura di S1.4d e coincidono tutti; le altre maschere native di R2 sono lette e mai scritte dal codice di R3b. Dopo il calcolo: 343/343 identici.
- **D-2 Stima di P-2** come in R2c: media sui nuclei focali appaiati (non media dei rapporti per finestra); l'intervallo è il bootstrap a due stadi pre-registrato.
- **D-3 φ** arrotondato a tre decimali e limitato a 1, come i φ per tessuto di S1.4d.
- **D-4 Flag `keep` della maschera CO** ricalcolati con la regola di R2 (`R2_common.nuclei_table`: centroide nel tessuto, fuori dalle bolle, area 4-400 µm²) sulle etichette CO, non ereditati dalle etichette di origine. Le etichette scartate restano nel raster come vicini, come nelle maschere di R2.
- **D-5 Stadio loo**: si stima solo M2 (l'unico usato dal nucleo simulato); la regressione G-3 confronta M2, ROI di test e seme con S1.4.
- **D-6 CP-1**: tutti i punti di Luca nelle 8 finestre cadono nel tessuto valido (0 esclusi, `results/R3b/R3b_cp1_inputs.csv`).
- **D-7 Contatti nella maschera a livello di ROI** (la prima parte di B-124 e B-126, descrittiva) non ricalcolati: per la maschera CO si riporta P-2, sugli stessi nuclei del giudizio di Luca.
- **D-8 Strumenti**: codice nuovo in `tools/R3b_*.py`, `tools/R3b_*.R`; nessun file degli step passati modificato; le funzioni copiate da R3, S1.4 e S1.4d sono verificate dalla regressione G-3.

## Risultati (dai file di `results/R3b/`, prima della revisione; i verdetti non cambiano dopo)

| previsione | tessuto | maschera | stima | ± 2 SE o IC | intervallo | esito |
|---|---|---|---|---|---|---|
| P-1 rapporto nuclei / punti di Luca | linfonodo | CO | 1.039 | 0.983-1.095 | 0.905-1.102 | confermata (attesa indeterminata) |
| P-1 | linfonodo | SD | 0.768 | 0.711-0.824 | 0.905-1.102 | smentita (nota) |
| P-1 | linfonodo | SR (rif.) | 1.043 | 0.980-1.105 | 0.905-1.102 | indeterminata |
| P-1 | cuore | CO | 1.088 | 1.018-1.158 | 0.878-1.134 | indeterminata (dichiarata) |
| P-1 | cuore | SD | 0.842 | 0.727-0.956 | 0.878-1.134 | indeterminata (nota) |
| P-2 contatti nella maschera, focali di R2c | linfonodo | CO | 0.746 | 0.621-0.862 | 0.35-0.75 | **indeterminata** (attesa 0.45-0.65) |
| P-2 | linfonodo | SD | 0.367 | 0.184-0.528 | 0.35-0.75 | indeterminata |
| P-2 | linfonodo | SR (rif.) | 0.983 | 0.931-1.000 | 0.35-0.75 | smentita (nota) |
| P-2 | cuore | CO / SD / SR | 0.019 / 0.022 / 0.038 | 0-0.080 / 0-0.095 / 0-0.111 | 0-0.070 | indeterminate (dichiarate) |
| P-3 contatti simulati | linfonodo | CO | 0.690 | 0.672-0.709 | 0.35-0.75 | confermata |
| P-3 | linfonodo | SD | 0.409 | 0.385-0.432 | 0.35-0.75 | confermata |
| P-3 | cuore | CO / SD | 0.047 / 0.012 | 0.037-0.056 / 0-0.027 | 0-0.070 | confermate |
| P-4 troncati | linfonodo | CO / SD | 0.023 / 0.031 | 0.017-0.029 / 0.024-0.037 | 0-0.05 | confermate |
| P-4 | cuore | CO / SD | 0.016 / 0.010 | 0.009-0.022 / 0.003-0.016 | 0-0.05 | confermate |
| CP-1 punti di Luca come centri | linfonodo | manuale + CO | 0.715 | 0.697-0.734 | 0.35-0.75 | confermata |
| CP-1 | cuore | manuale + CO | 0.015 | 0.009-0.021 | 0-0.070 | confermata |

**Regola di decisione**: linfonodo, CO unica eleggibile (SD smentita in P-1) → **CO**; cuore, entrambe eleggibili, P-3 confermata per entrambe, CP-1 non smentita → rapporto di P-1 più vicino a 1 (|log| 0.084 contro 0.172) → **CO**. Tetti φ: CO 0.868 (linfonodo), 0.147 (cuore); SD 0.649, 0.166.

**Cancelli prima della revisione**: G-1 26/26, G-2 5/5, G-3 identico in tutte le voci (nuclei R2 21/21 colonne × 10 ROI; generatori, celle, momenti R3, territori S1.4, M2, nuclei simulati S1.4d, finestre S1.4d 80/80; R2b 646 punti e 8 finestre; R2c 119 focali per maschera; φ di Space Ranger 0.816 e 0.180), G-4: invarianti della maschera CO 0 violazioni su 10 ROI, nucleo simulato 0 fuori, 0 sovrapposti, 0 invalidi su 4.76 milioni di nuclei, md5 343/343.

## Dopo il risultato e la revisione avversariale (aggiunto dopo i verdetti; non cambia previsioni, intervalli né verdetti)

- **D-9 Definizione di «troncato» in P-4** (revisore «ricalcolo», F-2): area finale < (1 − 1e-6) · area bersaglio, come P-2 di S1.4d; coincide con `nucleus_capped` in 40/40 righe rigenerate. Senza tolleranza metà dei nuclei risulterebbe «troncata» per la tolleranza del solutore (1e-7): definizione non scritta nella pre-registrazione, ora fissata.
- **D-10 G-4 fallito come scritto per una voce** (revisori F-1, F-05): «ogni centroide nel tessuto valido» non vale sulla maschera CO salvata (59 etichette su 145 098), ed è stato verificato solo dopo. Lo stesso accade nelle maschere native di R2 (Space Ranger 60, StarDist 107): le maschere grezze non sono filtrate, il filtro si applica a valle. Dopo il filtro usato da P-2: 0; fra i generatori: 0 (`results/R3b/R3b_posthoc_valid.csv`). **Errore di scrittura del cancello**: la proprietà vale sugli oggetti usati nei calcoli, non sul raster grezzo. Nessun numero cambia. Decisione di Luca sotto.
- **D-11 Batteria rinforzata** (F-01…F-04, F-09): test T12-T17 (pareggio dal lato Space Ranger, filtro del tessuto valido in geometria asimmetrica con troncamento, StarDist che esce dal poligono, StarDist non appaiato che non cancella uno Space Ranger, resto Space Ranger piccolo tenuto, appaiamento per sovrapposizione e non per IoU). G-1 34/34; mutanti M1-M5 5/5; mutanti del revisore 12/12 (prima 4/12).
- **PH-1 Contatti nella maschera CO nel linfonodo** (`tools/R3b_posthoc.py`): sui 59 focali appaiati la maschera segna 44 contatti, Luca 32; 15 contatti dove Luca vede nuclei separati, 3 mancati. Il 47 % dei focali (28/59) ha un vicino a contatto di forma Space Ranger; fra i 44 focali a contatto la quota è 64 % (F-07). Sensibilità del revisore (F-06, non pre-registrata): 8 contatti esistono solo sulla linea di taglio fra uno StarDist e il resto Space Ranger che copriva; senza di essi P-2 sarebbe 0.610 [0.467, 0.750]. StarDist sugli stessi focali: 18 contatti, 8 dove Luca vede separati, 12 mancati.
- **PH-2 Quanti poligoni Space Ranger rimasti sono nuclei** (`tools/R3b_posthoc_srprec.py`, scritta dopo la figura di controllo): nelle finestre del linfonodo 102 oggetti di forma Space Ranger, 65 appaiati a un punto di Luca (precisione 0.64, contro 0.97 dei nuclei di forma StarDist); nel cuore 67 e 47 (0.70 contro 0.87). Il conteggio di CO è giusto anche per compensazione fra circa 37 falsi positivi e circa 36 nuclei non visti (linfonodo).
- **PH-3 Aree dei nuclei di forma Space Ranger e modello di taglia** (`tools/R3b_posthoc_area.R`): 27 % delle celle interne in entrambi i tessuti; area mediana 12.8 contro 16.0 µm² (linfonodo), 15.2 contro 20.5 (cuore); sullo stesso nucleo il poligono Space Ranger è 1.34 volte lo StarDist (mediana, IQR 1.18-1.51). M2 su tutti i nuclei CO: σ 0.288, β 0.535 (linfonodo); sui soli nuclei di forma StarDist: σ 0.210, β 0.455, come oggi (0.212, 0.442). Parte della differenza d'area viene dai territori più piccoli di questi nuclei (28.9 contro 33.0 µm², F-08): lo script confronta mediane grezze.
- **PH-4 Sensibilità: taglia dai soli nuclei di forma StarDist** (`tools/R3b_posthoc_sdarea.R`): contatti simulati nel linfonodo 0.736 ± 0.016 (contro 0.690), troncati 0.028; cuore 0.049. φ 0.791 e 0.149.
- **Riferimenti ricalcolati** (`tools/R3b_refs.R`, `results/R3b/R3b_refs.csv`): lo stesso codice sulle uscite originali ridà le 20 righe di BIO_REFERENCES di linfonodo e cuore; con CO le misure dei territori cambiano al più dell'1.6 %; N/C mediano del linfonodo 0.606 → 0.480 (il vecchio valore usava le aree dei poligoni Space Ranger); φ del cuore 0.180 → 0.147.
- Note del revisore senza effetto: F-10 resti Space Ranger in più pezzi o sotto il 50 % (al più 0.2 % dei nuclei contati); F-11 lato della finestra in CP-1 (60.0 contro 59.96 µm, 0 centri nella striscia); F-12 φ stimato sui 5 ROI compreso quello di test (come la modalità «tissue» di S1.4d).

## Revisione avversariale (cancello G-5)

Due revisori indipendenti al commit `fdac76a`, con codice proprio sul server (`/mnt/micron/geo_spatialtrans/R3b/review_recalc/`, `review_code/`).

- **Ricalcolo**: 65 controlli su 65 uguali all'autore (scarto massimo 2.2e-16): maschera CO identica bit per bit su 10/10 ROI con conteggi identici; flag `keep`; P-1 (8 finestre × 3 maschere); P-2 (357 righe focale × maschera, intervalli bootstrap identici); P-3 e CP-1 (240 righe finestra × seme dai raster); P-4 rigenerando 4 unità × 10 semi (40/40); 24 esiti e regola di decisione. Limiti dichiarati: P-4 rigenerato su 4 unità su 20, aggregato dalle righe dell'autore. Rilievi: F-1 (→ D-10), F-2 (→ D-9), F-3 e F-4 informativi.
- **Codice**: nessun errore che cambi un numero o un verdetto; maschera ricostruita identica su 10/10 ROI; funzioni copiate identiche agli originali salvo i percorsi; CP-1 misura solo nuclei nati dai punti di Luca (0 di origine CO in 80 raster); md5 343/343. Rilievo principale: batteria debole (4/12 mutanti del revisore) → D-11 (12/12). Altri: F-05 (→ D-10), F-06 e F-07 (→ PH-1), F-08 (→ PH-3), F-10…F-12 note. Limiti dichiarati: catena R non rieseguita, ingressi CO controllati su 3 ROI, verdetti dai CSV.
- G-5: ricalcolo indipendente concordante → superato.

# In breve

## In breve: nel linfonodo e nel cuore, da quale insieme di nuclei veri prendiamo posizioni, aree e vicini?

**A cosa serviva.** Nel linfonodo e nel cuore il simulatore prendeva i centri dei nuclei da Space Ranger, ma aree e vicini da StarDist: tre misure da due maschere diverse. Abbiamo confrontato con i tuoi conteggi manuali due maschere uniche: «i nuclei di Space Ranger con la forma di StarDist» e «solo StarDist». Volevamo anche capire perché nel linfonodo i nuclei simulati si toccano più di quanto vedi tu.

**Cosa abbiamo trovato.**

1. **La regola sceglie in tutti e due i tessuti i nuclei di Space Ranger con la forma di StarDist: ne conta quanti ne conti tu, mentre StarDist da solo ne perde troppi.**
   - Nelle finestre di linfonodo questa maschera trova 416 nuclei contro i tuoi 404; StarDist da solo 314, quasi uno su quattro in meno. Nel cuore 262 e 202 contro i tuoi 242.
   - Il conteggio giusto viene in parte da una compensazione. Dei 102 nuclei che restano con la forma di Space Ranger (quelli che StarDist non vede) 65 sono nuclei veri, spesso pallidi, e 37 non corrispondono a nessun tuo punto. Nel cuore 47 veri su 67.
2. **L'eccesso di contatti simulati nel linfonodo non viene dai centri dei nuclei: con i tuoi punti come centri resta lo stesso.**
   - Con i tuoi punti i nuclei simulati a contatto sono 0.72, con la maschera scelta 0.69, contro il tuo 0.55; tutti dentro il tuo margine, che arriva a 0.75.
   - Lo 0.41 di StarDist dipendeva dai nuclei che StarDist non vede. La spiegazione che avevamo dato quando abbiamo portato il nucleo ritagliato nel simulatore (centri di Space Ranger troppo fitti) era sbagliata.
   - Nel cuore tutte le versioni stanno nel margine (0.05 contro il tuo 0.02).
3. **I nuclei rimasti con la forma di Space Ranger hanno aree e contatti poco affidabili.**
   - Sono più di un nucleo su quattro in entrambi i tessuti, con la forma a scalini di 2 µm. Dove lo stesso nucleo è visto anche da StarDist, il poligono di Space Ranger è più grande di un terzo, con grandi oscillazioni.
   - Se li usiamo per la taglia, la variabilità della taglia dei nuclei del linfonodo cresce di più di un terzo; se la stimiamo sui soli nuclei con la forma di StarDist resta quella di oggi. I contatti simulati nel linfonodo passano da 0.69 a 0.74, al bordo del tuo margine.
   - Sui nuclei che hai giudicato nel linfonodo, la maschera segna 44 contatti su 59 contro i tuoi 32: quasi due contatti su tre passano da un poligono di Space Ranger. Il riferimento del simulatore per i contatti resta comunque il tuo conteggio.
4. **Con la nuova maschera i territori non cambiano; cambiano il rapporto fra nucleo e cellula nel linfonodo e il tetto del nucleo nel cuore.**
   - Area e forma dei territori cambiano meno del 2 %.
   - Nel linfonodo il nucleo occupa in mediana il 48 % della cellula invece del 61 %: il vecchio valore usava le aree dei poligoni di Space Ranger. Nel cuore il tetto al nucleo scende dal 18 % al 15 % della cellula.
5. **Un controllo di integrità è fallito per come l'avevo scritto, non per un difetto della maschera.**
   - Avevo chiesto che ogni nucleo della maschera avesse il centro nel tessuto valido. Sulla maschera salvata 59 nuclei su 145 000 hanno il centro fuori, come nelle maschere originali di Space Ranger (60) e di StarDist (107).
   - Il filtro del tessuto si applica dopo: nessuno di quei nuclei entra nei territori o nei conteggi. I due revisori indipendenti hanno rifatto la maschera e tutti i numeri e trovano gli stessi valori.

**Cosa devi decidere.**

1. **Adottiamo per linfonodo e cuore i nuclei di Space Ranger con la forma di StarDist, ma stimiamo la taglia dei nuclei solo su quelli con la forma di StarDist? Consiglio di sì.**
   - Posizioni, territori e vicini vengono da un'unica maschera, che conta i nuclei come te. La taglia si misura solo dove la forma è affidabile, e il modello della taglia resta quello di oggi.
   - Alternativa: la taglia da tutti i nuclei della maschera. Più coerente sulla carta, ma un nucleo su quattro porta un'area a scalini e circa un terzo di questi non è un nucleo; la variabilità della taglia nel linfonodo crescerebbe di più di un terzo.
   - Alternativa: solo StarDist. Anche i suoi contatti simulati stanno nel margine (0.41), ma nel linfonodo ha un nucleo su quattro in meno, quindi territori più grandi di circa un terzo.
2. **Accettiamo che il controllo di integrità fallito sia un errore di come l'avevo scritto, e non un difetto della maschera? Consiglio di sì.**
   - Le maschere originali da cui partiamo hanno lo stesso comportamento, e il filtro del tessuto esclude quei nuclei da ogni calcolo.
   - Alternativa: filtrare già la maschera salvata e ripetere i controlli. Non cambierebbe nessun numero.

**Comunque**, appena decidi aggiorneremo i riferimenti di linfonodo e cuore (territori, taglia, tetto, rapporto fra nucleo e cellula) e la maschera di questi due tessuti nei banchi. I nuclei con la forma di Space Ranger che non sono nuclei veri li riguarderemo più avanti, per esempio tenendo solo quelli con un segnale di ematossilina. Nel linfonodo i contatti simulati stanno nella parte alta del tuo margine anche con i tuoi punti come centri: per stringere il confronto servirebbero più finestre contate a mano, e lo teniamo da parte. La forma vera delle cellule del cuore passa a quando scriveremo la composizione dei tessuti.

Appena decidi questi 2 punti, lo step è chiuso.

# Dubbi rinviati

| cosa | perché | cosa serve per chiuderlo | step destinatario | voce di BACKLOG |
|---|---|---|---|---|
| Forma vera delle cellule del cuore | richiede una segmentazione dei profili cellulari, non dei nuclei | misura indipendente dell'allungamento dei cardiomiociti | S1.5 (scelta di Luca, 2026-10-20) | BL-064 |
| Falsi positivi fra i poligoni Space Ranger rimasti: precisione 0.64 nel linfonodo, 0.70 nel cuore (contro 0.97 e 0.87 dei nuclei di forma StarDist) | il conteggio della maschera CO è giusto anche per compensazione; i falsi positivi dividono i territori vicini | filtro dei poligoni Space Ranger solo con evidenza di ematossilina (come gli «esclusivi» di R2), valutato sulle finestre manuali | da assegnare | BL-141 |
| Linfonodo: contatti simulati 0.69-0.74 (anche 0.72 con i punti di Luca) contro il manuale 0.55, margine 0.35-0.75 da 4 finestre | dentro il margine, ma sempre nella parte alta; il margine è largo | più finestre manuali di contatto nel linfonodo, oppure un nucleo meno «riempitivo» nei tessuti densi | da assegnare | BL-142 |

# Verifiche di Luca

- Perimetro, alla domanda «La voce sulla forma vera delle cellule del cuore (BL-064) è assegnata a R3b oppure S1.5. La includiamo in R3b?»: «Rinviarla a S1.5» (trascritto dalla chat, 2026-10-20).
- Ok sulle previsioni, sui loro intervalli e sulla regola di decisione: «Approvo così» (trascritto dalla chat, 2026-10-20).

# Riproduzione

```bash
cd ~/2025.geo_spatialtrans && git checkout <commit «R3b: chiusura»>
bash R/testing/run_R3b.sh test && bash R/testing/run_R3b.sh mutants
bash tools/run_R3b.sh
quarto render reports/R3b.qmd
```

[^roadmap-s7]: Perimetro di R3b e decisioni passate sul segmentatore di linfonodo e cuore
[^backlog]: Voci di backlog assegnate a R3b
[^r2b-match]: Oggetti dei segmentatori contro i punti manuali di Luca nelle 4 finestre per tessuto (già visto, non cieco)
[^b046]: Densità nucleare contata a mano, con il suo intervallo
[^b130]: Frazione di nuclei a contatto contata a mano da Luca, con intervallo bootstrap a due stadi
[^r2c-seg]: Contatti nelle maschere dei segmentatori sui nuclei focali di R2c (già visto, non cieco)
[^s14d]: Contatti del nucleo ritagliato e gonfiato sulle finestre, con centri Space Ranger e StarDist (già visto, non cieco)
[^protocol]: Protocollo di sessione versione 2
