Commit provato finale **T = `fea12bd`** (due commit provati prima: `038f150`, superato dalle correzioni della revisione; `9102831`, superato dal secondo giro). Log in `/mnt/micron/geo_spatialtrans/S2.8b/official3/` (`fea12bd`), `official2/` (`9102831`), `official/` (`038f150`). Tabelle della revisione in `reports/S2.8b_parts/`.

### Cancelli

| id | esito | valore letto dai log |
|---|---|---|
| G-1 | PASS | `fea12bd`: 392/392 riferimenti e voci di copertura (357 riferimenti, 35 di copertura, 4 documenti); mutanti 4/4 (il quarto, percorso troncato, aggiunto dopo la revisione T01; i tre pre-registrati 3/3) |
| G-2 (a) | **FAIL come scritto** | `cell_df` dell'avvio rapido con md5 `57b560944ae0ffe69987d2a73b5f1bb6` in due esecuzioni, non `9f361ea9…` di S1.6b; 0 ERROR; l'oggetto è `identical()` a quello di S1.6b e i due file differiscono in 1 byte, il decimo (versione di R che scrive il file: 4.6.1 in S1.6b, 4.6.2 oggi; R aggiornato da apt il 2026-11-05). Errore di pre-registrazione 1 |
| G-2 (b) | PASS | avvio rapido dello Step 2: md5 di `counts_2um` e `counts_8um` `9224353f024651e39ab20c3d46b77069` nelle due esecuzioni, controlli interni TRUE |
| G-2 (c) | PASS | due rese della guida nel clone: 5 impronte e 6 controlli, uguali fra loro e alla guida resa nel commit, tutti i controlli TRUE (al primo commit provato FAIL per l'errore dello strumento C-1) |
| G-3 | **FAIL come scritto** | batterie dello Step 2 21, 36, 39, 37, 40, 33, 78, 46 PASS e 0 FAIL; batteria dello Step 1 (`run_S1.6.sh test`): `GATE G-1 FAIL` per un solo controllo fuori dai FAIL ammessi, T7b di `test_S1.5d.R` (impronte del pannello di S1.5c salvate con R 4.6.1, vedi «Anomalie»); gli altri FAIL sono quelli noti di S1.6b (T13 di S1.5d, T12 di S1.5c, E10 ed E10b di S1.5). Errore di pre-registrazione 4 (stessa causa di G-2 a) |
| G-4 | PASS | `dev` = `step2-visium-grid` su GitHub (`78c88e9`), avanzamento dal vecchio `dev` `9250749` (0 commit propri, nessun commit di unione); `main` = `ec29fd6` come all'apertura; `R/`, `inst/`, `renv.lock`, `renv/`, `.Rprofile` e strumenti d'installazione identici a `48a5f18` salvo `R/testing/run_S2.8b.sh`; `DESCRIPTION` diverso solo in `Version` (0.2.0), `Date`, `Description`; dopo `fea12bd` solo `docs/`, `reports/`, `CHANGELOG.md` (`official3/G4_merge.log`; rifatto dopo il commit che registra questo esito) |
| G-5 | PASS dopo le correzioni | due revisori, tre giri; 149 controlli registrati, 42 discrepanze (una era un errore del revisore), tutte chiuse: corrette, dichiarate, o rinviate in BACKLOG (sotto); controllo OKF 0 documenti in errore |

### Controprove

- **CP-1** superata: con lo strumento di `fea12bd`, il README e il CLAUDE.md di `48a5f18` coprono 1 voce su 35 (solo `cell_layer`).
- **CP-2** superata: la copia della guida con `piece_um = 500` non si rende, e si ferma sull'errore del lato dei pezzi («multiplo di 56»); al primo commit provato la controprova accettava qualunque errore (revisione T07, corretto).

### Anomalie viste durante il lavoro

- **R aggiornato sul server** (4.6.1 → 4.6.2, 2026-11-05, apt): spiega G-2 (a) e la batteria dello Step 1 (controllo T7b di `test_S1.5d.R`, impronte del pannello di S1.5c scritte con 4.6.1): rimettendo a 4.6.1 il byte della versione, i tre livelli di T7b danno le impronte salvate (`/mnt/micron/geo_spatialtrans/S2.8b/probe/T7b.log`). Entra anche nelle impronte della ripresa (cartelle `out_dir` scritte prima non si riprendono). BL-257.
- **Lo strumento dei documenti di S1.6b non leggeva i blocchi ```` ```{r} ````**: rifatto con quello di S2.8b, la guida per Daniele dello Step 1 ha un link rotto verso [la decisione sulla composizione del linfonodo (DEC-015)](../docs/decisions/DEC-015.md), scritto con un percorso relativo sbagliato in `results/S1.6b/S1.6b_limits.md`; BL-258.
- **Lo script dei cancelli cambiato durante il calcolo**: al primo commit provato gli script giravano dalla copia di lavoro, che ho modificato mentre la batteria era in corso; bash ha letto lo script nuovo ed è uscito con un errore di sintassi alla fine. Dal terzo commit provato i modi del clone girano dalla copia dello script presa nel clone.
- **Carico del server**: durante il terzo giro altri processi (lavoro della tesi) hanno portato il carico a circa 35; la batteria è più lenta, i risultati no.

### Durata

Tre giri di prove nel clone: `038f150` 10:24-11:05, `9102831` 10:52-11:32, `fea12bd` 11:36-12:39 (2026-11-12); resa della guida circa 3 minuti (2 min 53 s, `probe/render.log`); batterie dello Step 2 da circa 2 a circa 16 minuti ciascuna nel terzo giro (orari dei log; il README dice «fino a circa 25 minuti»: da correggere, BL-258).
