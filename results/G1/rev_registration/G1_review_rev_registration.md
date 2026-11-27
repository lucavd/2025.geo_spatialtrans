# Revisione indipendente del cancello di registrazione Xenium -> Visium HD (X1, polmone 10x)

Revisore: rev_registration, 2026-11-27. Repo letto in sola lettura (`~/2025.geo_spatialtrans`, ramo `riesame`, HEAD 17010b3); codice e uscite in `/mnt/micron/geo_spatialtrans/G1_review/rev_registration/`.

## Conclusione

- **X1 passa il cancello, con lo stesso verdetto di S1.** Il residuo al 95° percentile dopo l'affine locale (riquadri non usati nella stima) è 0,97 / 0,87 / 0,80 µm in R1 / R2 / R3, contro i dichiarati 0,84 / 0,68 / 0,39 µm; la soglia è 2 µm. Nessun riquadro supera 1,35 µm.
- **Gli scarti rispetto a S1 (da 0,13 a 0,41 µm) sono rumore del metro.** In ogni regione, i due valori (il mio e quello di S1) cadono dentro la banda al 95% attesa dal solo rumore di ciascun metro. Il mio metro è un po' meno preciso (errore standard radiale 0,41-0,46 µm, contro 0,29-0,41 µm di S1).
- **Senza affine il residuo è di circa 10 µm, come in S1.** È una traslazione quasi costante su tutta la sezione, circa (2; 9,6) µm. Sta nella seconda tratta (H&E -> quadratini di Space Ranger), non nella prima: i nuclei Xenium cadono sull'ematossilina dell'H&E entro 0,5 µm. La correzione affine per area viene già applicata a valle (regola S5, cancello G-XR).
- **X2: le affermazioni reggono sulle tabelle.** In S1 X2 ha errore standard da 1,13 a 3,44 µm. In S2 nessun metro arriva a ≤ 0,5 µm in tutte e 3 le regioni. Il filtro per gene tiene in X2 dallo 0 al 34% degli UMI. Per X1 l'intervallo dichiarato «67-90%» è impreciso: i valori sono 42-90% su tutti i lati dei riquadri e 70-90% a 200-250 µm.
- **I dati sono intatti.** Gli md5 sono uguali prima e dopo (12 file), coincidono con `MD5_download.txt` per gli 8 file scaricati, e i 4 file estratti coincidono con i membri dell'archivio `binned_outputs.tar.gz`.

![Revisione del cancello X1]({{artifact:art_19b40e7a-be95-4e16-96f6-a87306468adf}})

Cosa guardare: (a) gli spostamenti per riquadro dei due metri stanno sulla diagonale; (b) i residui dopo l'affine sono lontani da 2 µm e dentro le bande del solo rumore; (c) una deformazione non affine imposta viene ritrovata e fa fallire il cancello.

## Cosa ho ricalcolato e come

Regola: `docs/sonda/REGOLA_S1.md`. Dati grezzi in `/mnt/micron/geo_spatialtrans/data_lung/X1/`: trascritti, pannello, matrice di allineamento, cellule e nuclei dello Xenium; `square_002um/filtered_feature_bc_matrix.h5` e `tissue_positions.parquet` (2 e 8 µm) del Visium HD; le due immagini H&E. Nessuno script del progetto è stato importato o eseguito. Gli script li ho letti solo per le convenzioni (nomi delle colonne, segno di s, coordinate delle regioni).

1. **Catena.** Xenium µm / 0,2125 -> inversa di `he_imagealignment.csv` -> pixel H&E -> quadratini da 2 µm. La mappa proiettiva pixel -> quadratini l'ho stimata con il mio codice (minimi quadrati lineari a 8 parametri): residuo massimo 5e-12 quadratini. Un'affine darebbe 0,22 quadratini. Variante: ricalcolata, «inversa, px 0,2125» porta il 62,2% delle cellule nella griglia e il 99,9% di queste su tessuto, come in S1.
2. **Regioni.** Ricostruite da zero con il criterio della regola (quadrati da 1 mm nella metà di sviluppo, copertura ≥ 0,9, R1/R2/R3 per distanza). Esito identico a S1: R1 (1500; 3800), R2 (2200; 1300), R3 (0; 5700) µm, coperture 0,909 / 0,917 / 0,945, 260 candidati, m = 3300 µm.
3. **Stimatore diverso da S1** (cambia solo la funzione di somiglianza e il raffinamento; tutto il resto segue la regola):
   - i trascritti Xenium si spostano in modo continuo, su un raster da 0,25 µm con immagine integrale, e si ricontano sulla griglia da 2 µm;
   - le mappe per gene si lisciano con una gaussiana di 2 µm;
   - la somiglianza è la media delle correlazioni di Pearson per gene, pesata con √(nX·nV);
   - il picco è il massimo su griglia: 1 µm entro ± 12 µm dallo spostamento grossolano, poi 0,25 µm entro ± 1 µm, senza quadrica.

   Lo spostamento grossolano è mio: mappe totali lisciate, ricerca entro ± 100 µm, passo 2 µm. Restano come nella regola: L = 200 µm, filtro ≥ 10 conteggi per gene, validità ≥ 500, regola del bordo, affine a 6 parametri a scacchiera (stima sui pari, valutazione sui dispari, e viceversa), p95 con interpolazione lineare.
4. **Precisione.** Metà dei trascritti Xenium a caso e Visium assottigliato con una binomiale 1/2 (seme mio, 4242 + indice della regione). Errore standard radiale = rms(|s_A − s_B|/2).
5. **Test sintetico** prima dei dati veri: spostamenti noti (0; 0), (3,3; −1,7), (7,6; 9,1) µm ritrovati entro 0,25 µm; su tessuto sbagliato il picco finisce al bordo.

## Tabella dichiarato / ricalcolato

| quantità | ambito | dichiarato | ricalcolato | scarto | ok |
|---|---|---|---|---|---|
| p95 dopo affine locale (µm) | R1 | 0,84 | 0,97 | 0,13 | sì (≤ 2; dentro il rumore) |
| p95 dopo affine locale (µm) | R2 | 0,68 | 0,87 | 0,19 | sì |
| p95 dopo affine locale (µm) | R3 | 0,39 | 0,80 | 0,41 | sì (S1 al bordo inferiore della sua banda) |
| p95 senza affine (µm) | R1 / R2 / R3 | 10,44 / 9,04 / 11,22 | 10,34 / 9,32 / 11,52 | 0,10 / 0,28 / 0,31 | sì (fallisce senza affine in entrambi) |
| riquadri validi | R1 / R2 / R3 | 25 / 24 / 25 | 25 / 24 / 25 | 0 | sì (stessi riquadri) |
| regioni (angolo, µm) | R1-R3 | (1500;3800) (2200;1300) (0;5700) | identiche | 0 | sì |
| errore standard radiale a L = 200 (µm) | R1 / R2 / R3 | 0,41 / 0,34 / 0,29 | 0,46 / 0,44 / 0,41 | metro diverso | sì (≤ 0,5: stessa scelta L = 200) |
| errore standard radiale a L = 100 (µm) | R1 / R2 / R3 | 1,19 / 0,66 / 0,97 | 0,96 / 0,51 / 0,84 | metro diverso | sì (> 0,5: L = 100 scartato anche qui) |
| controllo positivo | 9 casi | 9/9 entro 0,3 µm | 9/9 + 9/9 (spostamenti fuori raster), errore ≤ 0,1 µm | | sì |
| controllo negativo, entro 2 µm | R1 / R2 / R3 | 0 / 0 / 0,04 | 0 / 0 / 0 (scambio di regione); 0 / 0 / 0 (permutazione spaziale) | | sì |
| geni condivisi | X1 | 286 (testo della regola) | 284 | 2 | no, solo nel testo (vedi discrepanze) |
| errore standard X2, metro S1 (µm) | 200-250 µm | 1,1-3,4 | 1,13-3,44 | < 0,05 | sì |
| X2, metri di S2 ≤ 0,5 µm in tutte e 3 le regioni | A2, A4, B2, B4 | nessuno | nessuno (minimo dei massimi 0,85 µm, A4 a 250 µm) | | sì |
| quota di UMI tenuta dal filtro, X2 | mediane per regione e L | 0-34% | 0,0-34,2% | < 0,001 | sì |
| quota di UMI tenuta dal filtro, X1 | mediane per regione e L | 67-90% | 42-90% (tutti gli L); 70-90% (L 200-250); 41-85% dal mio codice sulle metà | | no, intervallo impreciso |

Tabelle: [X1_gate_summary_review.csv](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/772e0ff0-ba7e-4810-b226-90d4e01db151/v1b4341c5_X1_gate_summary_review.csv), [X1_tiles_vs_S1.csv](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/9597f0cd-3d67-4f87-a3da-99e945a73b2f/vfbfe724a_X1_tiles_vs_S1.csv), [X1_se_by_L_review.csv](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/0ba04389-184f-449f-81dd-f891da0a3c6b/vf19ede52_X1_se_by_L_review.csv), [X1_positive_controls_review.csv](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/1acd8362-a682-4307-bd8e-8d928ab93305/vb78b2686_X1_positive_controls_review.csv), [X1_negative_controls_review.csv](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/02fcabfa-c480-4b66-acc8-3c213c7a5f19/v458bd660_X1_negative_controls_review.csv), [X1_warp_control_review.csv](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/1c1a1bfd-17dc-4130-8316-861889056c31/vfc54787b_X1_warp_control_review.csv), [X1_he_leg1_review.csv](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/6321256f-7ccf-40b9-aa13-b972ce51db84/v01717efd_X1_he_leg1_review.csv), [X1_gene_filter_share_review.csv](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/84149b0e-414f-48e2-ae34-5559016a3829/v4857c488_X1_gene_filter_share_review.csv), [X2_se_check_review.csv](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/ef774518-0a08-462d-95aa-235dfb882935/v1d181024_X2_se_check_review.csv), [md5_before.txt](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/6a287f87-5c36-4819-93b6-fc660d18ddbf/v51dd3aba_md5_before.txt), [md5_after.txt](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/6a287f87-5c36-4819-93b6-fc660d18ddbf/v51dd3aba_md5_before.txt).

## Controlli del mio metro

- **Positivo, traslazione.** Spostamenti imposti (3; 0), (0; 4), (−3,5; 3,5) µm: ritrovati esattamente in 9/9 casi. Questi valori sono multipli del raster da 0,25 µm, quindi verificano solo il segno e la contabilità. Con valori fuori raster, (3,1; 0), (0; 4,3), (−3,6; 3,45) µm: 9/9 entro 0,1 µm per asse. Il p95 dopo l'affine resta fra 0,64 e 0,96 µm: l'affine assorbe la traslazione, come atteso.
- **Positivo, deformazione non affine** (controllo aggiunto: verifica che il cancello possa fallire). Spostamento in y pari a 4·sin(2πx/400 µm). La media attesa per riquadro alterna ± 2,55 µm da una colonna all'altra; quella ritrovata correla 0,99 con l'attesa, con errore mediano 0,2-0,45 µm. Il p95 dopo l'affine sale a 3,7 / 3,7 / 4,3 µm: il cancello fallisce in tutte e 3 le regioni.
- **Negativo.** Xenium di una regione contro il Visium della regione successiva, oppure Visium permutato nello spazio dentro ogni riquadro: 0% dei riquadri entro 2 µm in tutte le regioni. Fra il 32 e il 68% dei riquadri finisce al bordo; i restanti hanno un p95 di 14-24 µm. Il picco di somiglianza mediano scende da 0,24-0,27 a 0,01-0,02.

## Discrepanze e spiegazioni

1. **p95 dopo l'affine più alto del dichiarato di 0,13-0,41 µm.** Spiegazione: precisione dei metri. Ho simulato 4000 volte il solo rumore, con gli stessi riquadri e lo stesso adattamento a scacchiera, all'errore standard misurato. Bande al 95%: mio metro 0,61-1,11 / 0,59-1,05 / 0,54-0,99 µm; metro di S1 0,55-0,99 / 0,45-0,83 / 0,38-0,70 µm. Ogni valore osservato cade nella banda del suo metro; il 0,39 di S1 in R3 è al bordo inferiore. Per riquadro, le differenze fra i due metri hanno mediana ≤ 0,26 µm e deviazione standard 0,22-0,46 µm. Il residuo dopo l'affine misura quindi il rumore del metro: a scala di 200 µm non si vede alcuna deformazione non affine. Il verdetto non cambia.
2. **«286 geni» nel testo di REGOLA_S1, 284 nei dati.** I 289 identificativi Ensembl del pannello sono tutti nella matrice grezza; 5 di essi (BCAS1, CD8B, DUOX1, FCN1, S100A7) mancano nella matrice filtrata, che è quella prescritta. Le regole successive (S5, S8) scrivono già 284. Effetto sul cancello nullo: quei geni non hanno UMI nella matrice usata.
3. **«67-90%» per X1.** Le mediane per regione e L, sulle metà, sono 0,42-0,67 a L = 100 e 0,70-0,90 a L = 200-250 (tabella `diag_counts_half_summary.csv`, che riproduco esattamente dalla tabella per riquadro). L'intervallo 67-90% non corrisponde a nessun sottoinsieme naturale. Il messaggio (X2 tiene molto meno di X1) resta vero.
4. **X2, sfumatura.** Una singola combinazione ha errore standard ≤ 0,5 µm: metro B2, R3, L = 250, 0,43 µm. Nessun metro però lo raggiunge in tutte e 3 le regioni, e in R2 solo 8 riquadri su 16 sono usabili. «Non giudicabile» resta corretto secondo la regola di S2.

## Diagnosi in più (nessuna scelta ne dipende)

- **Le due H&E** (quella dello Xenium e quella del Visium HD) sono identiche pixel per pixel nei ritagli delle 3 regioni: 0 pixel diversi su 48,5 milioni per regione.
- **Prima tratta** (Xenium -> H&E, matrice di 10x). Ho confrontato i centri dei nuclei Xenium con il canale dell'ematossilina, su una griglia da 1 µm. Spostamento per regione (−0,46; −0,13), (−0,27; −0,21), (−0,44; 0,37) µm, con deviazione standard fra riquadri ≤ 0,19 µm. Lo scarto di circa 10 µm sta quindi nella seconda tratta (H&E -> quadratini).
- **Un'affine unica per l'intera sezione** basta: rotazione −0,46 mrad, scala entro 150 ppm, traslazione circa (1,9; 9,6) µm. Lasciando fuori una regione alla volta, il p95 sulla regione esclusa è 0,89 / 1,02 / 0,96 µm. Anche la sola traslazione per regione dà p95 ≤ 0,98 µm.

## Cosa è girato

Server lesexp-server, 4 processi a un thread al massimo, Tctl finale 75,6 °C. Tempi:

- setup: 31 s;
- 33 esecuzioni per riquadro: circa 20 minuti in tutto, da 14 a 271 s ciascuna;
- 9 esecuzioni per i controlli fuori raster: circa 4 minuti;
- diagnosi H&E: 1 minuto.

Semi: 4242 (metà), 777 (permutazione), 12345 (sottocampione per la mappa pixel -> quadratini), 20261127 (bootstrap), 99 (simulazione del rumore).

Codice: [g1_common.py](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/4fef6d81-5195-40cd-ac71-234b2c455c4f/v7fa6f1a0_g1_common.py), [g1_setup.py](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/abf43a0c-7b8b-4922-8053-4df6200ceee7/vc2ecde13_g1_setup.py), [g1_tiles.py](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/6a86bf89-60e0-4e89-8b42-df2be86baf9c/v3616ff55_g1_tiles.py), [g1_he.py](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/6efab186-6075-45c3-be11-52223ad50748/v70bf1be9_g1_he.py), [test_synth.py](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/eaae51f9-1be9-4fab-8701-81ac4bac272d/vb315d4fc_test_synth.py), [run_tiles.sh](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/0a85fc9b-4b78-4834-bb67-3956a7459d2b/vc8e3ca3a_run_tiles.sh), [run2.sh](/Users/utente/.claude-science/orgs/8ef7b99f-0db6-4f39-9db4-aa2bf2f215a5/artifacts/proj_fcb66815598c/9eba7989-503b-4f64-a588-915b165470ef/v7f3c2312_run2.sh).

Limiti:

- Il mio metro usa un raster da 0,25 µm: quantizzazione ≤ 0,125 µm per asse.
- L'intervallo bootstrap del p95 sui riquadri (R1 0,66-0,98, R2 0,60-1,35, R3 0,63-0,87 µm) ignora la correlazione spaziale fra riquadri; il suo estremo superiore coincide con il residuo massimo.
- Nella diagnosi H&E la convenzione del centro del pixel è ambigua di mezzo pixel (0,14 µm).
