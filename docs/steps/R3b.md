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

# In breve

(da scrivere a numeri verificati)

# Dubbi rinviati

| cosa | perché | cosa serve per chiuderlo | step destinatario | voce di BACKLOG |
|---|---|---|---|---|
| Forma vera delle cellule del cuore | richiede una segmentazione dei profili cellulari, non dei nuclei | misura indipendente dell'allungamento dei cardiomiociti | S1.5 (scelta di Luca, 2026-10-20) | BL-064 |

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
