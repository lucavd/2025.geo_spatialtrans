---
type: Rule
title: "S5, A3 sullo Xenium: il simulatore di oggi contro lo Xenium in quadratini, sui geni del pannello"
description: "Regola scritta e committata prima di leggere l'espressione di X1 e prima di restringere i bracci al pannello: braccio Xenium in quadratini, metro sul pannello, cancelli, controprova, attese, esito."
status: stable
generated: { by: geo-agent/progetto-master, at: 2026-11-24T11:30:00+01:00 }
sources:
  - id: paper
    path: docs/PAPER.md
  - id: s1
    path: docs/sonda/REGOLA_S1.md
  - id: s3
    path: docs/sonda/REGOLA_S3.md
  - id: s3b
    path: docs/sonda/REGOLA_S3b.md
  - id: s3c
    path: docs/sonda/REGOLA_S3c.md
  - id: xen
    path: docs/riesame/DECISIONI_B_XENIUM.md
  - id: orient
    path: results/sonda/S5/orient.json
---
# Perché

Il contratto (`docs/PAPER.md`[^paper]) chiede A3: il simulatore batte lo Xenium in quadratini sulla stessa sezione, sui geni del pannello, con uno scarto D_xen − D_sim ≥ δ = 1 h. Lo Xenium in quadratini (trascritti Xenium veri messi sulla griglia da 2 µm, come in ENACT) è il modo più semplice oggi per avere dati simili a Visium HD con la verità per cellula: se il simulatore non lo batte, il lettore non ha motivo di usarlo. Questa sessione dà la prima misura di A3. Oggi il simulatore è il banale B (`REGOLA_S3.md`[^s3]), quindi A3 oggi = D_p(XB) − D_p(B).

Decisione di Luca, 2026-11-24 in chat: «Oggi solo A3; proroga della sonda da 7 a 9». A5 e A4 (con A2) avranno una sessione ciascuna. La proroga si registra dopo il riesame che lo strumento chiede per una seconda proroga.

Solo X1. X2 non ha passato la registrazione (S2: non giudicabile), quindi A3 si misura sulla sola X1 e lo si dichiara come limite.

# Già noto prima del commit (esplorazione senza espressione)

`tools/sonda/S5_explore_orient.py` e `S5_explore_local.py` (sola H&E hires e maschere in_tissue; `results/sonda/S5/orient.json`[^orient], `local_ecc.json`):

- La griglia di array di X1 ha lo stesso orientamento di quella di H2 (rotazione −0,03°, scala 1,002; H1 era ruotata di 180°). Dice delle maschere 0,95.
- Le 10 aree di V riportate in X1 (quadrati di 750 µm allineati agli assi, come le aree di T in S3b[^s3b]) cadono tutte nella metà di sviluppo di X1 (x < m − 50 µm, m = 3300 µm, `REGOLA_S1.md`[^s1]); frazione del quadrato vero coperta ≥ 0,9967; tessuto di X1 ≥ 0,93.
- Residuo locale X1 -> H2 (ECC di sola traslazione su tessere da 1 mm): mediana 8,2 µm, 95° percentile 19,2 µm; controllo positivo di 40 µm ritrovato a 39,5 µm. La correlazione di fase sulle stesse tessere dava 240 µm di mediana: è un errore di quello strumento, non della registrazione (l'ECC e la sovrapposizione delle immagini concordano).
- X2: ruotata di circa 5,7° rispetto a H2; non usata.
- Geni del pannello: 289 identificativi Ensembl nel pannello Xenium di X1, 284 fra le feature della matrice filtrata di H2 (e di H1). G_p = questi 284.

# Bracci (tutti ristretti a G_p, stesse 10 aree, stessa catena)

- Vp, Tp, Bp, Pp: le esecuzioni esistenti di V (H2, `S3_sonda/runs/V`), T (H1 omologo, `S3b_sonda/runs/T`), B (2 semi, `S3_sonda/runs/B`) e P (permutato, `S3_sonda/runs/P`) con X.parquet e genes.tsv ristretti a G_p (stesso ordine dei geni di genes.tsv); quadratini, nuclei e cells.parquet copiati. Nessuna nuova simulazione: per B la restrizione di una multinomiale è la multinomiale sul profilo ristretto, alla profondità ristretta.
- XB, Xenium in quadratini (braccio di A3): sezione X1, aree omologhe di questa regola. Quadratini = quelli da 2 µm di X1 con centro nel quadrato dell'area (in_tissue e in_reg come `S3_export.py`). Trascritti Xenium con qv ≥ 20 sui geni di G_p, portati sulla griglia con la catena di S1 (composita della variante scelta), più la correzione affine locale dell'area (sotto); ogni trascritto va al quadratino del centro più vicino. label = nucleo Xenium (`nucleus_boundaries.parquet`, stessa trasformazione) che contiene il centro del quadratino, tutti i nuclei (XD-9[^xen]); cells.parquet con cell_type vuoto (tipi trasferiti dalla catena, come V). Diradamento come XD-11: f per tessuto = media sulle 5 aree di Vp degli UMI per quadratino in_reg / la stessa media di XB grezzo; ogni conteggio (quadratino, gene) diventa Binomiale(n, f), `default_rng([20261130, i])`, i = 0...9 (tumore 1...5, stroma 1...5). Se f > 1 in un tessuto, XB resta alla sua profondità e lo si dichiara.
- XV, Visium HD dopo Xenium su X1 (descrittivo): stesse aree di XB, esportazione di `S3_export.py` su X1 (StarDist sull'H&E di X1, conteggi della matrice grezza), ristretta a G_p. Serve a separare l'effetto della tecnologia (XB contro XV, stessa sezione) da quello della sezione (XV contro Vp). Nessuna decisione ne dipende.

# Correzione affine locale di XB (per area, prima di esportare)

Metro di S1 invariato (`tools/sonda/S1_registration.py`: mappe per gene standardizzate, riquadri L = 200 µm, finestra ± 12 µm attorno allo spostamento grossolano, quadrica, validità ≥ 500 conteggi), su area + 240 µm di margine, Visium HD di X1 sui geni di G_p (solo metà di sviluppo, controllo di metà attivo). Affine a 6 parametri dei centri dei riquadri verso centri + s, stimata su tutti i riquadri validi e non di bordo; il segno si verifica con il cancello G-XR (riquadri rifatti dopo la correzione). Residuo dichiarato: scacchiera di S1 (stima sui pari, valutazione sui dispari e viceversa). Limite dichiarato: la correzione usa l'espressione del Visium HD di X1, non quella di H2; nessuna misura M-1...M-6 di X1 entra nella stima.

# Metro sul pannello

Invariato salvo i geni: catena di `REGOLA_S3.md` (bin2cell, RCTD doppietto seme 20261118, riferimento HLCA `L` invariato; RCTD e misure usano i geni comuni fra dati e riferimento, quindi G_p), tabella di `S3_table.py` con genes = "panel", h_p con la costruzione di S3b (κ = 1, B = 2000, `default_rng([20261115, i, j])`) sulle 5 aree di Vp, D_p con lo stimatore appaiato di `REGOLA_S3c.md`[^s3c] (funzioni di `S3c_D.py` senza modifiche), bootstrap per coppie stratificato per tessuto (B = 2000, seme 20261123), intervallo (ii) con f = 1,584. M-4 si toglie per un tessuto se la media dei singoletti di Vp è ≥ 0,99, deciso sulla tabella di Vp prima di D. Confronti con la regola di S3c: uguale se |Δ| + w ≤ 1 h, diverso se |Δ| − w > 1 h.

- Tetto del pannello D_p(T), pavimento del pannello D_p(B); A1 resta sui geni completi (non cambia qui).
- A3 = D_p(XB) − D_p(B), con w e se dalla distribuzione bootstrap della differenza (stessi indici per tutti i bracci). Regge se A3 − 2 se ≥ 1 (regola dello strumento).

# Cancelli (se uno fallisce i numeri non si leggono)

| id | controllo | esito richiesto |
|---|---|---|
| G-XR | registrazione Xenium -> X1 per area, dopo la correzione: residuo sulla scacchiera, 95° percentile; riquadri rifatti dopo la correzione, mediana di \|s\|; controllo positivo di S1 (spostamenti (1; 0), (0; 2,5), (1,5; −1,5) µm) | ≤ 2 µm in 10/10 aree; ≤ 0,5 µm in 10/10; ritrovati entro 0,3 µm per asse in 10/10 |
| G-SUB | restrizione: per ogni esecuzione di Vp, Tp, Bp, Pp, conteggi = originale ristretto a G_p (somma e posizioni identiche) | 50/50 |
| G-HALF | ogni quadratino di X1 usato (XB, XV, correzione) ha x < m − 50 µm | 0 eccezioni |
| G-MET | C-1p, C-4 (P contro T «diverso», P peggio) e C-5 (B0 contro B1 non «diverso») di S3c sulla tabella del pannello | passano |
| G-MD5 | md5 dei file letti (dati di X1, H1, H2, esecuzioni di S3, S3b, tabelle) prima e dopo | uguali |
| G-REV | un sotto-agente ricalcola dalla tabella h_p, D_p dei bracci e A3 con codice proprio | entro 1e-9 |

# Controprova

CP-XB: lo Xenium in quadratini è tessuto vero; deve stare più vicino a V del controllo permutato: D_p(XB) < D_p(P) con XB contro P «diverso». Se fallisce, il braccio XB è costruito male (registrazione, nuclei o diradamento) e A3 non si legge.

# Attese (scritte prima)

| id | attesa | giudizio |
|---|---|---|
| P1 | f fra 0,1 e 0,6 in entrambi i tessuti (Xenium più profondo di H2 sui geni del pannello) | confermata se dentro |
| P2 | tetto del pannello D_p(T) fra 0 e 1,0 h | sull'intervallo, come S4 |
| P3 | pavimento del pannello D_p(B) fra 2 e 4 h | sull'intervallo |
| P4 | D_p(XB) fra 1,0 e 3,0 h: oltre il tetto (manca lo spostamento dei trascritti e la cattura di Visium) ma sotto B (ha la mappa dei tipi vera, il 69 % del divario in S4) | sull'intervallo |
| P5 | A3 oggi < 0: il simulatore di oggi (il banale) non batte lo Xenium; A3 non regge | confermata se l'estremo superiore di A3 è < 0 |

Precisione attesa: w circa 0,3-0,5 h, come in S3c e S4.

# Esito e conseguenze

- Tutti i cancelli e CP-XB passano: A3 = D_p(XB) − D_p(B) e il suo se entrano in `results/misura.json` (sezione affermazioni) attraverso il check smoke, che rilegge `results/sonda/S5/D.json` e verifica l'identità. S5 si chiude come misura su A3.
- Livello obiettivo, scritto ora: se D_p(XB) − D_p(T) ≤ 1 h, nessun simulatore può battere lo Xenium di 1 h (non può scendere sotto il tetto): A3 è irraggiungibile come è scritta e il contratto va riposizionato a G1. Si scrive come allarme sul criterio; decide Luca.
- Un cancello fallisce: numeri provvisori, niente in `results/misura.json`; la sessione lo dichiara.

# Calcolo

Sul server al più 12 processi a un thread (file `WORKERS` di `/mnt/micron/geo_spatialtrans/S5_sonda`), sotto il tetto di 50 core e il limite termico. Pilota su un'unità per braccio con i tempi; se la corsa completa (70 esecuzioni: Vp 10, Tp 10, Bp 20, Pp 10, XB 10, XV 10) supera 4 ore ci si ferma e si chiede. Uscite in `/mnt/micron/geo_spatialtrans/S5_sonda/`, un file per unità, rilanciabile.

[^paper]: docs/PAPER.md
[^s1]: docs/sonda/REGOLA_S1.md
[^s3]: docs/sonda/REGOLA_S3.md
[^s3b]: docs/sonda/REGOLA_S3b.md
[^s3c]: docs/sonda/REGOLA_S3c.md
[^xen]: docs/riesame/DECISIONI_B_XENIUM.md, XD-9 e XD-11
[^orient]: results/sonda/S5/orient.json
