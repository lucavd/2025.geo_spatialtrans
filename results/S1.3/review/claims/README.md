# S1.3 — Revisione delle affermazioni (RA-claims)

Revisore: sotto-agente «revisore delle affermazioni», 2026-10-09. Base: commit `37e68ad`, branch `step1-cell-layer`.
Ho scritto solo in questa cartella. Nessuna modifica al pacchetto, ai risultati, al report o a git.

## Metodo
1. **Ricalcolo dei verdetti** (`01_recompute.R`, codice mio, senza `source` di `tools/S1.3_analysis.R`). Legge `S1.3_test_*.csv`, i 40 rds reali, i 1 800 rds nulli, i 30 rds cp3, i csv `cp1/cp2/cp2b/d1/perf/c10a_syn/c2_coverage`, `R3_null_summary.csv` e `R3_roi_summary.csv` (primario: spaceranger in A4/A6, cellpose_rgb altrove). Ricostruisce le 188 righe di `S1.3_verdicts.csv` e le confronta campo per campo (`rc_compare_verdicts.csv`).
   Le previsioni di check B le ho trascritte dal **testo** della pre-registrazione, non dal codice.
   Controprove a valle:
   - C-5 ricalcolato dalle tabelle per cellula di R3 (`/mnt/micron/geo_spatialtrans/R3/real/*_cells.parquet`), con il cv_loc dalla lambda_loc di R3;
   - riassunti dei nulli ricalcolati dalle celle salvate (repliche 1–2).
2. **Diagnostica mirata** (`02_diagnostics.R`): scomposizione di K-2/K-3, attribuzione delle cause di C-5 e C-6, intervalli citati nel report, caso A5 r5.
3. **Previsioni codificate e testo**: confronto riga per riga (`predictions_coding_check.csv`).
4. **Affermazioni del report**: 54 affermazioni controllate contro i dati (`report_claims_check.csv`).
5. **Tempi**: `git log` e `git diff 72126a5 37e68ad`; mtime delle uscite grezze.

## Comandi per rifare
```bash
cd ~/2025.geo_spatialtrans
R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.3/review/claims/01_recompute.R > results/S1.3/review/claims/rc_log.txt 2>&1   # ~3 min
R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.3/review/claims/02_diagnostics.R > results/S1.3/review/claims/rc_diag_log.txt 2>&1
```
File prodotti:
- `rc_my_verdicts.csv`, `rc_compare_verdicts.csv`: ricalcolo dei verdetti e confronto;
- `rc_roi_real.csv`, `rc_c5_from_cells.csv`, `rc_c6.csv`, `rc_K2_decomposition.csv`, `rc_engine_choice.csv`: check C, C-5, C-6 e scelta del motore;
- `rc_arbitration.csv`, `rc_checkB_roi.csv`, `rc_G1_predictions_uncoded.csv`, `rc_ks_nsides.csv`, `rc_cp3.csv`: arbitro, check B, previsioni G1, KS, CP-3;
- `rc_perf_ratio.csv`, `rc_mutants.csv`, `rc_null_summary_from_cells.csv`, `rc_aux.txt`: prestazioni, mutanti, riassunti dei nulli, numeri ausiliari.

## Esiti
### 1. Ricalcolo dei verdetti
**188/188 righe identiche** in valore, esito, previsione e conferma. L'unica differenza testuale è una mia annotazione voluta: CP-3 riga 2, previsione vacua.
Le controprove concordano:
- C-5 dalle celle concorda con i vettori memorizzati in 40/40 ROI;
- i riassunti dei nulli ricalcolati coincidono con quelli memorizzati in 180/180 file (max rel 0).

### 2. Previsioni codificate e testo
Le soglie, la regola ≥ 4/5 ROI e tutte le previsioni G2 e G3 (comprese le 8 d'ordine), CP-1…CP-3 e D-1 coincidono col testo. Differenze:
- **G1 non codificato**, benché la pre-registrazione scriva «valgono le colonne CSR»: 23/24 confermate, A5 B-1 smentita (RA-claims-07);
- **C-5 con rel_summary in più** rispetto al testo: geos 9/40 invece di 11/40, esito invariato (RA-claims-06);
- CP-3 riga 2 confermata in modo vacuo (RA-claims-10);
- previsioni dell'addendum non codificate (RA-claims-17);
- interpretazioni minori su ordine, D-1 e M1 (RA-claims-21).

### 3. Testo del report
I numeri sono quasi tutti corretti: 40 delle 54 affermazioni sono confermate (alcune con nota); le altre sono discrepanze, sopra-interpretazioni o numeri non verificabili (`report_claims_check.csv`). Problemi:
- un numero sbagliato: «31 volte più veloce»; il dato è 131× (RA-claims-05);
- descrizione imprecisa delle discordanze dell'arbitro (RA-claims-04);
- sopra-interpretazioni nella sintesi del check B (RA-claims-08) e in CP-3, «uscite dalla regione» non misurate (RA-claims-09);
- un sottoinsieme chiamato «regioni piccole» che non lo è (RA-claims-13);
- attribuzioni di C-6 presentate come accertate (RA-claims-11);
- tre numeri senza file di supporto (RA-claims-15);
- commit di riproducibilità incompleto (RA-claims-20).

Le analisi post hoc (CP-2b, copertura per punti, sensibilità Bsens) sono dichiarate come tali. La deviazione 6 è dichiarata, ma senza la sua entità (RA-claims-16).

### 4. Percorso della scelta del motore
L'addendum è stato scritto prima del confronto e la regola è stata applicata correttamente dal codice: nessun motore eleggibile. Il report, però, **non riporta la conseguenza che la regola stessa prescrive** («S1.3 non chiusa; si torna al design»). Scrive invece che il confronto «non ha prodotto una scelta automatica» (RA-claims-02).

L'argomento «K-2/K-3 cadono per cause indipendenti dal motore» regge solo in parte (RA-claims-03):
- **K-2**: la parte dovuta ai frammenti è davvero comune ai motori (142/142 flag diversi stanno su celle con frammenti). Però:
  - anche un C-5 ristretto alle celle senza frammenti fallirebbe, per geos in 2/40 ROI e per deldir in 1/40;
  - deldir cade anche su C-1..3 in 5 ROI.
- **K-3**: sulle stesse 90 tassellazioni geos passa 90/90 e deldir 75/90. Il fattore comune è lo strumento di misura (l'unione), con un tasso che dipende dal motore.
- **Omissioni**:
  - deldir cade anche su K-1;
  - la clausola di arbitrato codificata (K_arb) dà FALSE anche per GEOS, per un caso di conteggio dei lati (RA-claims-04).

Nel merito, i dati favoriscono GEOS in modo netto:
- arbitro 25/26, di cui 16 casi con area sbagliata da deldir;
- C-perf 131×;
- K-3 a parità di ingressi 90/90 contro 75/90;
- eventi di robustezza sui ROI reali 91 contro 312.

La scelta è quindi ben fondata. È il percorso che non è riportato fedelmente.

### 5. «Chiusa con riserve» o «non chiusa»?
**Secondo la pre-registrazione la sessione è «non chiusa».** L'addendum lo stabilisce senza ambiguità per il caso «nessun eleggibile», e il caso si è verificato (riga `scelta` = FAIL). A questo si aggiungono:
- due check C pre-registrati falliti (C-5, C-6), attribuiti a un errore di definizione: per C-5 l'attribuzione è provata, per C-6 è un'inferenza;
- C-2 misurato con l'unione, fallito su 15/1 800 nulli;
- una controprova senza potenza (CP-2);
- una controprova smentita (CP-3).

Nel merito, invece, la funzione appare corretta:
- 286/286 asserzioni e 6/6 mutanti;
- CP-1 confermata;
- C-1 (Σ aree) entro 1e-9 ovunque;
- copertura per punti senza lacune né doppie coperture alla risoluzione disponibile;
- previsioni di check B confermate 36/38 (G1 23/24).

Le cadute residue riguardano gli strumenti di misura e le definizioni, non il codice.

Le strade compatibili coi principi del progetto sono due. Spetta a Luca sceglierne una:
- **(a) «Non chiusa».** Una sessione S1.3b breve, senza codice nuovo, pre-registra prima di eseguire:
  - C-5 sulle sole celle senza frammenti, oppure contro un R3 ricalcolato con D-S1.3.2;
  - C-6 con una prova diretta dell'attribuzione;
  - C-2 misurato con predicati o copertura;
  - la verifica dei 3 ROI con residui non-frammento (A1 r5, A4 p3).
  I numeri di S1.3 sono già noti, quindi questi criteri vanno dichiarati post hoc. Il valore della ripetizione sta nel misurare anche le cause oggi inferite.
- **(b) «Chiusa con riserve, in deroga alla regola 4 dell'addendum, per decisione del responsabile».** La deroga va scritta nel Verdetto, fra le deviazioni e in ROADMAP §6, con data e motivazione. Le riserve vanno estese a: la regola 4, K_arb = FALSE per GEOS, le previsioni dell'addendum smentite e C-6 attribuito per inferenza.

Senza l'una o l'altra correzione, il verdetto attuale non è coerente con la pre-registrazione. Io propendo per (a), perché la pre-registrazione aveva previsto proprio questo caso. Accetto (b) se la deroga è esplicita.

## Limiti di questa revisione
- Non ho ricalcolato geometrie: non ho rieseguito tassellazioni, l'arbitro o la copertura per punti. Ho lavorato sulle uscite salvate e sulle tabelle per cellula di R3.
- Le celle dei nulli sono salvate solo per le repliche 1–2. Il ricalcolo dei riassunti dalle celle copre quindi 180/1 800 file; gli altri 1 620 sono verificati solo come ingresso dei verdetti.
- Tre numeri del report non stanno in nessun file e non li ho potuti verificare.
