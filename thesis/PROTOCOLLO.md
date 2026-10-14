# Protocollo di tesi — Pierangela

**La fonte dei profili di espressione cambia le conclusioni di un benchmark?**
*Studio di simulazione pre-registrato (ADEMP) su dati di trascrittomica spaziale sintetici*

Scuola di specializzazione in Statistica medica · supervisore: Luca Vedovelli
Branch: `thesis/expression-sources` del repository `lucavd/2025.geo_spatialtrans`
Periodo: **19 ottobre – 31 dicembre 2026**, a tempo pieno. Il 31/12/2026 la tesi è **scritta e consegnata** a Luca.
Versione del protocollo: 1.0 (2026-10-14)

---

## 0. Come usare questo documento

Questo documento è la tua guida per tutta la tesi. Ogni settimana ha un obiettivo, un prodotto da consegnare e un criterio che dice quando è "fatta" (sezione 6).

Quando un dubbio non trova risposta qui, la decisione si prende con Luca e si scrive in `thesis/DECISIONS.md` con la data. Una decisione che non è scritta lì non vale. Vale anche per le decisioni di Luca.

Due principi del progetto valgono anche per te:
- **Non siamo infallibili, né tu né io.** Ogni affermazione ha bisogno di dati, di un grafico, di una conferma e di una controprova.
- **Nessun numero "a memoria".** Un valore entra nella tesi solo con una fonte: un DOI, oppure lo script e il dataset che lo hanno prodotto.

## 1. La domanda in breve

Un simulatore di dati single-cell o spaziali deve decidere da dove prendere l'espressione "vera" di ogni tipo cellulare. Le scelte possibili sono quattro:
- un riferimento scRNA-seq reale;
- un modello parametrico costruito da zero;
- un simulatore model-based (scDesign3);
- un modello generativo addestrato su archivi pubblici (GEM-1 di Synthesize Bio).

Crowell et al. (2023) hanno mostrato che la scelta del simulatore incide sui confronti fra metodi di analisi single-cell. Nessuno l'ha ancora misurato per la **deconvoluzione spaziale**, cioè la stima della composizione in tipi cellulari di ogni spot.

**Obiettivo 1 — Fedeltà.** Quanto i profili per tipo cellulare di ciascuna fonte si avvicinano a profili reali **held-out** dello stesso tessuto e dello stesso tipo cellulare?

**Obiettivo 2 — Sensibilità del benchmark.** La classifica e l'errore dei metodi di deconvoluzione cambiano con la fonte? Quale fonte dà la classifica più vicina a quella ottenuta su pseudo-spot costruiti con **cellule reali held-out** (l'àncora reale)?

Il risultato serve al progetto principale: deciderà quale fonte di profili userà il simulatore nella fase successiva (Step 2).

## 2. Decisioni prese

| Data | Decisione |
|---|---|
| 2026-10-14 | Tema, titolo e i due obiettivi della sezione 1 |
| 2026-10-14 | Periodo dal 19/10 al 31/12/2026, a tempo pieno; tesi scritta e consegnata entro il 31/12/2026 |
| 2026-10-14 | Bracci dello studio: E, P, M, G, C−, E′ (sezione 5.2) |
| 2026-10-14 | GEM-1 (braccio G) è **opzionale**: entra solo se il gate G0 passa entro venerdì 30/10 (sezione 9) |
| 2026-10-14 | Tre tessuti: colon/tumore del colon, tessuto linfoide, cervello |
| 2026-10-14 | Metodi di deconvoluzione: RCTD, CARD, cell2location (sulla GPU del server), NNLS come baseline. SPOTlight è facoltativo |
| 2026-10-14 | Un unico modello di rumore binomiale negativo per tutte le fonti; livello spaziale fatto di pseudo-spot con composizione Dirichlet; il cell layer del simulatore principale **non** si usa |
| 2026-10-14 | Split dei dati sempre per donatore o campione, mai per cellula |
| 2026-10-14 | Account Synthesize Bio personale di Pierangela; la richiesta di permesso di pubblicazione la invia lei, con Luca in copia |
| 2026-10-14 | Branch `thesis/expression-sources` con push riservato; nessun merge in altri branch durante la tesi |
| 2026-10-14 | Regola di taglio in caso di ritardo: si toglie prima SPOTlight, poi il braccio G, poi il tessuto cervello. Non si tagliano mai i bracci E, P, M, C− ed E′ né i tessuti colon e linfoide |

## 3. Regole di lavoro (non negoziabili)

1. **Pre-registrazione.** Ipotesi, soglie, numero di repliche e misure di prestazione si scrivono **prima** di guardare i risultati, nel documento `thesis/protocol/ADEMP.md`, approvato da Luca con la data. Le attese non si cambiano dopo il risultato. Se il risultato le smentisce, la discrepanza si scrive ed è un risultato.
2. **Una sessione, un obiettivo.** Ogni sessione di lavoro (T1…T6, sezione 6) si chiude con un report Quarto a 4 sezioni (sezione 7). Idee e problemi laterali vanno in `thesis/BACKLOG.md`, non nel codice.
3. **Controprova sempre.** Per ogni risultato chiediti quale test lo smentirebbe se fosse falso, e fallo.
4. **Revisione avversariale.** Prima di chiudere un report, una persona (o un agente) diversa da chi ha scritto il codice ricalcola i numeri chiave dagli script e dai dati. Le discrepanze trovate si riportano nel report.
5. **Riproducibilità.** Ogni figura e ogni tabella si rigenera con un solo comando. Seed fissi e dichiarati. Nessun file modificato a mano fra uno script e l'altro.
6. **Unità statistica.** Il donatore o campione è l'unità per lo split. Lo spot simulato è l'unità per le prestazioni dentro una replica. La replica Monte Carlo è l'unità per l'incertezza.
7. **I dati reali sono in sola lettura.** Non si modificano i file scaricati: si derivano nuovi file in `data/derived/`.

## 4. Ambiente di lavoro

### 4.1 Repository e branch
- Repository: `https://github.com/lucavd/2025.geo_spatialtrans`.
- Il tuo branch è `thesis/expression-sources`, nato da `dev`. **Puoi fare push solo su questo branch.** Gli altri branch sono protetti: non puoi aggiornarli né fare merge.
- **Tutto il tuo lavoro sta nella cartella `thesis/`.** Non modificare file fuori da `thesis/`: il resto del repository è il simulatore principale e ti serve al massimo come lettura.
- A fine tesi Luca valuterà, in una revisione dedicata, cosa portare nel progetto principale.

### 4.2 Sul server (lesexp-server)
- Worktree del branch: `/mnt/micron/geo_spatialtrans/tesi_pp1209/repo`.
- Dati (fuori da git): `/mnt/micron/geo_spatialtrans/tesi_pp1209/data/` con le sottocartelle `raw/` (scaricati, sola lettura), `derived/` e `sim/` (output delle simulazioni).
- Il server è condiviso. Usa al massimo **32 core** e lancia i job con `nice`. La GPU (RTX 4090) si usa per cell2location: prenota con Luca gli slot nelle settimane 5–8.
- R 4.6 con pacchetti da Posit Package Manager (binari, nessuna compilazione). Gli script si lanciano con `Rscript --vanilla` e `.libPaths()` esplicito in testa, con la libreria della tesi in `thesis/renv/`.
- Python per cell2location: ambiente `uv` dedicato in `thesis/.venv-c2l/`, separato da quello del progetto principale.
- I job lunghi scrivono su percorsi stabili e devono poter ripartire da dove si sono fermati: il server a volte perde la rete.

### 4.3 Struttura della cartella `thesis/`
```
thesis/
  PROTOCOLLO.md          questo documento
  DECISIONS.md           decisioni datate
  BACKLOG.md             idee e problemi rimandati
  protocol/ADEMP.md      protocollo pre-registrato (settimana 2)
  R/                     funzioni (una responsabilità per file)
  scripts/               pipeline numerate: 01_download.R, 02_split.R, ...
  reports/T1.qmd ...     report di sessione (sezione 7)
  writing/               capitoli della tesi
  tests/                 test delle funzioni
  synthesize_log/        query e risposte GEM-1 (sezione 9)
```

## 5. Disegno dello studio (ADEMP; Morris, White & Crowther 2019)

### 5.1 Aims
Gli obiettivi 1 e 2 della sezione 1.

### 5.2 Data-generating mechanisms
Si varia **un solo fattore principale**: la fonte della matrice delle medie μ (geni × tipi cellulari).

| Braccio | Fonte di μ | Ruolo |
|---|---|---|
| E | Medie empiriche dallo split di training dello scRNA reale | Riferimento "onesto" |
| P | Parametrico de novo: medie da una distribuzione gamma, marcatori con fold-change controllati | Massimo controllo, minimo realismo |
| M | scDesign3 addestrato sullo split di training | Via di mezzo model-based |
| G | GEM-1-sc in modalità "mean estimation", da metadati CL/UBERON | Modello generativo (opzionale, gate G0) |
| C− | Profili E con le etichette dei geni permutate dentro ogni tipo cellulare | Controllo negativo: deve fallire |
| E′ | Medie da un secondo donatore reale indipendente | Controllo positivo: il tetto ottenibile |

Fattori secondari, in griglia fissa:
- tessuto: colon/tumore del colon, tessuto linfoide, cervello;
- cellule per spot: 2 livelli;
- profondità di sequenziamento: 2 livelli.

I livelli esatti si fissano nell'ADEMP della settimana 2, partendo dai dati reali.

Generazione dei dati di una replica:
1. per ogni spot si estrae una composizione da una Dirichlet e si campionano le cellule;
2. per ogni cellula di tipo k, il conteggio del gene g segue una binomiale negativa con media μ_gk (scalata per la library size) e dispersione φ_g stimata dai dati reali di training. **Il modello di rumore è lo stesso per tutti i bracci**;
3. si sommano le cellule dello spot.

**Àncora reale:** pseudo-spot costruiti sommando **cellule reali** del donatore held-out, con la stessa composizione Dirichlet. La verità è nota, perché si conosce il tipo di ogni cellula.

### 5.3 Estimands
- Primario: la proporzione di ciascun tipo cellulare in ogni spot.
- Secondario: il recupero dei marcatori per tipo cellulare (obiettivo 1).

### 5.4 Methods
- RCTD (pacchetto spacexr);
- CARD;
- cell2location (Python, GPU);
- NNLS (baseline);
- SPOTlight (facoltativo).

Ogni metodo riceve come riferimento lo scRNA dello split di training, sempre lo stesso per tutti i bracci.

### 5.5 Performance measures
- **Obiettivo 1:**
  - correlazione delle log-medie fra fonte e held-out, per tipo cellulare;
  - sovrapposizione dei top-k marcatori;
  - separabilità dei tipi cellulari.
- **Obiettivo 2:**
  - RMSE e divergenza di Jensen–Shannon fra proporzioni stimate e vere;
  - classifica dei metodi per ciascun braccio;
  - concordanza delle classifiche con l'àncora reale (Kendall τ) e fra bracci (W di Kendall).
- **Incertezza:** errore standard Monte Carlo per ogni misura. Il numero di repliche si sceglie nell'ADEMP sulla precisione desiderata della misura principale, usando il pilota della settimana 2.

### 5.6 Controprove previste
1. **C−** deve dare fedeltà bassa e classifiche diverse dall'àncora. Se le classifiche non cambiano nemmeno con profili privi di senso, il benchmark è insensibile al livello di espressione: è un risultato e va scritto.
2. **E′** è il tetto. Una fonte sintetica che lo supera indica un probabile leakage (lo stesso donatore in training e in test, o dati di GEM-1 sovrapposti). Va investigata prima di andare avanti.
3. **Split per donatore** verificato con un test automatico: nessun identificativo di donatore compare in due split.

## 6. Timeline settimana per settimana

Ogni riga dice cosa consegni e quando è "fatta". La tesi si scrive **mentre si lavora**: alla fine di ogni settimana c'è un pezzo di testo in `thesis/writing/`.

| Sett. | Date | Sessione | Lavoro | Prodotto ed è "fatta" quando… |
|---|---|---|---|---|
| 1 | 19–23 ott | **T1 Avvio** | Mezza giornata con Luca su server, git e worktree. Letture (sezione 10). Account Synthesize Bio e mail di richiesta (sezione 9). Inventario dei dataset scRNA pubblici per i tre tessuti, con almeno 2 donatori indipendenti ciascuno; verifica delle licenze | `thesis/data_inventory.csv` (dataset, DOI/URL, donatori, tipi cellulari annotati, tecnologia, data di deposito). Uno script di prova gira dal worktree. Bozza di 2–3 pagine dell'introduzione |
| 2 | 26–30 ott | **T2 Protocollo** | Download e split per donatore. Pilota dichiarato (i suoi numeri non entrano nei risultati) per misurare tempi, costo GPU di cell2location e variabilità. **Esito del gate G0 entro venerdì 30/10.** Scrittura dell'ADEMP | `thesis/protocol/ADEMP.md` **approvato da Luca con la data**: ipotesi, soglie, livelli dei fattori, numero di repliche giustificato con l'errore standard Monte Carlo. Esito di G0 in DECISIONS. Il capitolo Metodi è in bozza (è l'ADEMP riscritto) |
| 3–4 | 2–13 nov | **T3 Obiettivo 1** | Costruzione delle matrici μ per E, P, M, (G), C−, E′; misure di fedeltà contro l'held-out nei tre tessuti | Report `T3.html` a 4 sezioni con revisione avversariale. Bozza della sezione Risultati 1 |
| 5–6 | 16–27 nov | **T4 Simulazione** | Generatore degli pseudo-spot e àncora reale; test (conservazione dei conteggi, composizioni vere salvate, split rispettato). Lancio delle repliche | Report `T4.html`: check di correttezza superati, controprova C− sul generatore. Repliche complete su disco |
| 7–8 | 30 nov–11 dic | **T5 Obiettivo 2** | Deconvoluzione con i 4 metodi su tutte le repliche (cell2location sulla GPU); misure di prestazione, classifiche, concordanza con l'àncora | Report `T5.html` con revisione avversariale. Bozza della sezione Risultati 2 |
| 9 | 14–18 dic | **T6 Stesura** | Discussione, limiti, conclusioni; figure finali rigenerate da un solo comando | **Bozza completa della tesi a Luca entro venerdì 18/12** |
| 10 | 21–23 dic | — | Revisione sui commenti di Luca | Seconda versione |
| — | 24–27 dic | — | Pausa | — |
| 11 | 28–31 dic | — | Ultima rilettura, controllo che ogni numero del testo corrisponda a un output; pacchetto di riproducibilità | **Tesi consegnata il 31/12/2026.** Branch con tag `tesi-v1.0` |

**Punti di controllo con Luca:** un incontro breve ogni lunedì (cosa è fatto, cosa blocca) e un'approvazione formale al termine di T2 (protocollo) e prima di T6 (risultati chiusi).

**Se sei in ritardo:** il ritardo si dichiara al lunedì, non a fine mese. Si applica la regola di taglio della sezione 2 e si scrive in DECISIONS.

## 7. Report di sessione (Quarto, `thesis/reports/Tn.qmd` → HTML)

Ogni report ha quattro sezioni fisse, più una chiusura:
1. **Pre-registrazione.** Cosa ti aspettavi e con quali soglie, scritto prima del run, con la fonte.
2. **Check C — correttezza del software.** I test che il codice fa ciò che deve (dimensioni, conservazione dei conteggi, seed riproducibili, split per donatore). Tabella PASS/FAIL.
3. **Check B — plausibilità biologica.** I profili e i dati simulati hanno senso? Marcatori noti nei tipi giusti, distribuzioni di library size e di zeri confrontate con il reale, per tessuto. Qui la tua esperienza da biologa è decisiva.
4. **Controprova.** Il test che smentirebbe il risultato se fosse falso (C−, E′, un'alternativa), con il suo esito.

**Chiusura:** verdetto (chiusa / chiusa con riserve / non chiusa) ed esito della revisione avversariale.

## 8. Struttura della tesi scritta

1. **Introduzione** (8–10 pagine): trascrittomica spaziale e deconvoluzione; il ruolo della simulazione nei benchmark; Crowell 2023; i modelli generativi; la domanda.
2. **Metodi** (10–12 pagine): l'ADEMP completo; dati e split; generatore; metodi di deconvoluzione; misure di prestazione ed errore standard Monte Carlo; software e riproducibilità.
3. **Risultati** (12–15 pagine): Obiettivo 1, Obiettivo 2, controprove.
4. **Discussione** (6–8 pagine): cosa implicano i risultati per chi costruisce benchmark; limiti (pseudo-spot senza geometria, rumore uguale per tutti i bracci, tre tessuti); cosa resta da fare.
5. **Appendici:** protocollo pre-registrato con la data di approvazione; deviazioni dal protocollo; log delle query a GEM-1; istruzioni per riprodurre tutto.

Ogni numero nel testo deve avere accanto lo script o la tabella da cui viene. In settimana 11 si fa un controllo sistematico numero per numero.

## 9. Synthesize Bio / GEM-1 — istruzioni

**Settimana 1:**
1. Crea un account personale su app.synthesize.bio e genera una API key. Non condividerla e non metterla in git: va in `~/.Renviron` sul server, che è fuori dal repository.
2. Scrivi a support@synthesize.bio (Luca in copia). Chiedi se è consentito usare GEM-1 in una tesi e in una pubblicazione scientifica che confronta i suoi profili con dati reali e con altri simulatori. Chiedi anche la data di cutoff dei dati di addestramento dei modelli single-cell e i limiti di uso dell'API.

**Gate G0 (entro venerdì 30/10).** Il braccio G entra se **tutte** queste condizioni sono vere:
- c'è un permesso scritto a pubblicare;
- l'accesso all'API funziona entro la quota disponibile;
- il vocabolario del modello copre i tessuti e i tipi cellulari dei tre dataset (l'elenco si scarica dalla pagina dei metadati di Synthesize Bio);
- è nota la data di cutoff. Per il braccio G servono donatori held-out depositati **dopo** il cutoff, altrimenti si misura la memoria del modello.

Se una condizione manca, il braccio G esce dallo studio, si scrive in DECISIONS e il resto prosegue senza cambiamenti.

**Cose da sapere (pacchetto rsynthbio 4.2.0, verificato sul sorgente CRAN):**
- il modello single-cell `gem-1-sc` funziona solo in modalità `"mean estimation"`: restituisce profili medi, senza rumore di misura. Per noi va bene, perché il rumore lo aggiungiamo noi, uguale per tutti i bracci;
- i Termini d'uso prevedono che input e output possano essere usati dall'azienda per addestrare i suoi modelli: invia **solo metadati pubblici**;
- ogni chiamata va registrata in `thesis/synthesize_log/` (in git) con data, `model_id`, la query JSON completa e l'md5 della risposta. La risposta va salvata in `data/derived/synthesize/` (parquet o rds, fuori da git). Il modello può cambiare nel tempo, e il log è l'unico modo per riprodurre il braccio G.

## 10. Letture

Obbligatorie, settimana 1:
- Morris TP, White IR, Crowther MJ. *Using simulation studies to evaluate statistical methods.* Statistics in Medicine 2019;38:2074–2102. doi:10.1002/sim.8086. È la base dell'ADEMP: leggila tutta, con particolare attenzione alle sezioni sugli estimandi, sulle misure di prestazione e sull'errore standard Monte Carlo.
- Crowell HL, Morillo Leonardo SX, Soneson C, Robinson MD. *The shaky foundations of simulating single-cell RNA sequencing data.* Genome Biology 2023;24:62. doi:10.1186/s13059-023-02904-1 (vedi anche l'Author Correction, doi:10.1186/s13059-024-03329-0).
- Il preprint di GEM-1: *Generative genomics accurately predicts future experimental results.* bioRxiv 2025. doi:10.1101/2025.09.08.674753. Leggilo in modo critico: gli autori sono fondatori dell'azienda.

Da trovare e leggere in settimana 1, scegliendo con Luca:
- l'articolo di scDesign3;
- gli articoli originali di RCTD, CARD e cell2location;
- almeno un benchmark indipendente di metodi di deconvoluzione spaziale.

Per ciascuno scrivi in `thesis/writing/letture.md` cinque righe: domanda, dati, risultato principale, limite, cosa ci serve.
