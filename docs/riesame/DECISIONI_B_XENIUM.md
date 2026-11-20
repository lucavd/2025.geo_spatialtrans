---
type: Decisions
title: "Riesame, Parte B: scelte di attuazione dei bracci Xenium (X, X1, XC, XC2)"
description: "Scelte tecniche imposte dal codice per preparare le esportazioni dei bracci Xenium, scritte prima di calcolare i bracci; nessuna guardando D."
status: draft
generated: { by: geo-agent/skill-v2, at: 2026-11-20T12:27:00+01:00 }
sources:
  - id: regola
    resource: "/riesame/REGOLA_B.md"
    title: "Riesame, Parte B: regola di esito scritta prima del calcolo"
---
# Scelte di attuazione dei bracci Xenium

Scritte il 2026-11-20 alle 12:27, prima di calcolare i bracci X, X1, XC e XC2 e prima di aprire i dati Xenium scaricati (lo zip era ancora in download). Non committate: la consegna di questo sotto-agente chiede di non committare perché altri lavorano sul repository; l'impronta md5 di questo file, presa prima del calcolo, è nel registro delle uscite (`/mnt/micron/geo_spatialtrans/riesame/B/decisioni_md5_before.txt`). Nessuna scelta guarda D né misure dei metodi. Ogni riga: data, scelta, motivo.

## Dati e geni

| id | data | scelta | motivo |
|---|---|---|---|
| XD-1 | 2026-11-20 | Trascritti dello zip GEO GSE280314 di P2 (`transcripts.parquet`), qualità `qv >= 20`, solo i 408 geni del pannello dentro i 18 024 comuni. Simbolo -> Ensembl con `gene_panel.json` (identificativo del pannello); LRMP entra come ENSG00000118308. Controlli negativi, codici non assegnati e sonde dei recettori T esclusi. | REGOLA_B, braccio X; i 408 identificativi sono `riesame/xenium/panel_ids_in_common.txt`. |
| XD-2 | 2026-11-20 | `genes.tsv` di ogni esportazione = i 408 identificativi Ensembl nell'ordine di `panel_ids_in_common.txt`; `X.parquet` con indici 0-based su quell'ordine. | Stessi geni per tutti i bracci (REGOLA_B). |

## Registrazione Xenium -> Visium HD

| id | data | scelta | motivo |
|---|---|---|---|
| XD-3 | 2026-11-20 | Trasformazione affine 2D dalle coordinate Xenium (µm, `x_centroid`, `y_centroid`) alle coordinate del H&E a piena risoluzione di V (`x_centroid_visium_scale` = colonna in pixel, `y_centroid_visium_scale` = riga in pixel), stimata sul file pubblicato da 10x (`Xenium_P2_cell_info.csv.gz`, GitHub 10XGenomics/HumanColonCancer_VisiumHD) con minimi quadrati robusti (pesi di Huber, soglia 1,345 volte la deviazione assoluta mediana dei residui, 20 iterazioni). | Registrazione pubblicata, non nostra; la trasformazione affine è la più semplice che la riproduce. |
| XD-4 | 2026-11-20 | Regola scritta prima: se il 95° percentile dei residui della trasformazione globale supera 2 µm, si usa una trasformazione affine locale per area, stimata con le sole cellule dentro il riquadro dell'area allargato di 250 µm. Qualità riportata in µm (pixel per 0,2737 µm/pixel) con la quota di cellule entro 2, 5, 10 µm. | La registrazione di 10x potrebbe non essere affine; una regola fissata prima evita di scegliere guardando i risultati. |
| XD-5 | 2026-11-20 | Prima di usare la trasformazione si verifica che il file GitHub e lo zip GEO siano nello stesso sistema di coordinate Xenium: abbinamento per `cell_id` e confronto dei centroidi; se gli identificativi non coincidono, distanza dal centroide GEO più vicino. Se non coincidono entro 1 µm per la maggior parte delle cellule, ci si ferma e lo si riporta. | Il file GitHub ha 1 105 615 righe contro 340 837 cellule dello zip GEO (letto il 2026-11-20 dalla sola intestazione e dal conteggio delle righe): potrebbe venire da un'altra segmentazione. |
| XD-6 | 2026-11-20 | Controlli indipendenti, solo riportati: densità dei nuclei Xenium per mm² nelle 10 aree contro i nuclei di V (`n_nuclei_keep` delle esportazioni di F3.1); quota di quadratini dell'area con il centro dentro una cellula Xenium contro la quota `in_tissue` e `in_reg` di V; correlazione incrociata fra la mappa degli UMI di V sui 408 geni e quella dei trascritti Xenium, a blocchi di 8 µm, con lo spostamento del picco. Nessuna scelta dipende da questi numeri, salvo un errore grossolano (picco oltre 16 µm), che si riporta e ferma il lavoro. | Controllo della registrazione con dati diversi da quelli che l'hanno stimata. |

## Griglia, nuclei e formato (braccio X)

| id | data | scelta | motivo |
|---|---|---|---|
| XD-7 | 2026-11-20 | Quadratini = quelli dell'esportazione di V di F3.1 (`F3.1/main/real/<area>/bins.parquet`): stessi `bin_id`, `row`, `col`, `x_um`, `y_um`, `in_tissue`, `in_filtered`, `in_reg`, copiati. Cambiano solo `label` e `label_dimg`. | Stessi indici dei quadratini di V, come chiede la consegna; le misure si calcolano sugli stessi quadratini dei due lati (scelta di S2.5). |
| XD-8 | 2026-11-20 | Trascritto -> quadratino: per ogni area, affine (`array_row`, `array_col`) -> (`pxl_row_in_fullres`, `pxl_col_in_fullres`) stimata ai minimi quadrati sui quadratini del riquadro (`tissue_positions.parquet` di V, sola lettura); il trascritto portato in pixel si riporta in (riga, colonna) frazionarie con l'inversa e si arrotonda all'intero più vicino; si conta solo se il quadratino è nel riquadro dell'area. I conteggi si scrivono per tutti i quadratini del riquadro che ne ricevono, non solo per quelli `in_filtered`. | Il quadratino contiene il punto se il centro più vicino nella griglia è il suo; con una mappa affine l'arrotondamento è esatto. |
| XD-9 | 2026-11-20 | `label` = nucleo Xenium (`nucleus_boundaries.parquet`, poligono portato in pixel di V con la trasformazione di XD-3/XD-4) che contiene il centro del quadratino; 0 altrimenti. Tutti i nuclei Xenium, senza il filtro `keep` di R2 (che è della segmentazione di V). Numero del nucleo = posizione (da 1) della cellula in `cells.parquet` dello zip; `cells.parquet` dell'esportazione elenca i nuclei con almeno un quadratino nel riquadro, con `cell_type`, `classifiable`, `cell_type_r6` vuoti come nelle aree di verifica di F3.3 (tipi trasferiti dopo dalla catena comune) e una colonna in più `xenium_cell_id`. | Formato delle aree di verifica di F3.3 (`tools/F3.3_export.py`, ramo non IN). |
| XD-10 | 2026-11-20 | `label_dimg` = `label`. | Lo scarto `d_image` di V descrive lo scarto fra i trascritti di V e il suo H&E; nello Xenium trascritti e nuclei hanno la stessa registrazione. |

## X1: diradamento

| id | data | scelta | motivo |
|---|---|---|---|
| XD-11 | 2026-11-20 | Un fattore f per tessuto (tumore, stroma), uguale per le sue 5 aree: f = media sulle 5 aree di V degli UMI per quadratino `in_reg` sui 408 geni, divisa per la stessa media dello Xenium (braccio X). Bersaglio di V dalle esportazioni di F3.1 (`F3.1/main/real/<area>/X.parquet`) ristrette ai 408 geni. Ogni conteggio (quadratino, gene) di X diventa Binomiale(n, f). Se f > 1 il diradamento non è possibile e lo si riporta senza esportare X1. | «Un diradamento binomiale unico» (REGOLA_B): un numero per tessuto; l'unità della misura è l'area, quindi la media si fa sulle medie delle aree. `in_reg` è l'insieme dei quadratini «in tessuto» delle misure di F3.1 e F3.3. |
| XD-12 | 2026-11-20 | Seme di X1: `numpy.random.Generator(PCG64(SeedSequence([20261124, i])))`, i = posizione dell'area nell'elenco A2_r1...A2_r5, A3_r1...A3_r5 (0...9). | Seme dichiarato, uno per area. |

## XC e XC2: spostamento, cattura, profondità

| id | data | scelta | motivo |
|---|---|---|---|
| XD-13 | 2026-11-20 | Parametri dal pacchetto senza modificarlo: `transcript_params(A)` (r, sigma_in, sigma_e, D, w, ell), `capture_params(A)` (strisce, b, w2, ell2), `stripe_ceiling(capture_params(A), 1e-5)` (tetto M*), letti con `Rscript` da `R/05b3`, `R/05b4`, `R/05b6`. | «Parametri del preset del tessuto» (REGOLA_B, braccio XC). |
| XD-14 | 2026-11-20 | Nascita = posizione Xenium del trascritto (qv >= 20, 408 geni), portata in µm del riquadro dell'area come i quadratini (`x_um`, `y_um`). Il campo di variabilità fine (S2.3) è spento: descrive dove nascono i trascritti, e qui le nascite sono quelle vere. Il rapporto nucleo/citoplasma r non si usa per lo stesso motivo. | Nello Step 2 il campo e r entrano nella nascita, non nello spostamento (`place_and_move_transcripts()`, passo 1). |
| XD-15 | 2026-11-20 | Trascritti nati nel riquadro dell'area allargato di 100 µm per lato; quelli che dopo lo spostamento cadono fuori dai quadratini del riquadro sono persi. | I vicini veri fuori dall'area portano trascritti dentro; con la prima coda (Gamma(2, ell = 24,7 µm)) la quota che si sposta oltre 100 µm è circa lo 0,35 % dei trascritti. |
| XD-16 | 2026-11-20 | Spostamento come `place_and_move_transcripts()` (codice riscritto in Python riga per riga, non nel pacchetto): distanza d dal bordo del tessuto sulla sottogriglia di 0,5 µm, tessuto = centro del sottopixel dentro le regioni di S2.3 dell'area (`S2.3/layers/<area>_regions.tsv`, le stesse di `in_reg` di V) con il riquadro aperto (fuori dal riquadro vale tessuto, come `frame_open = TRUE` di F3.1); con probabilità w coda (raggio Gamma(2, ell), angolo uniforme), altrimenti gaussiana con sigma_e con probabilità exp(-d/D), con sigma_in altrimenti; poi seconda coda con probabilità w2 (raggio Gamma(2, ell2) dal punto di nascita), ripiegata per riflessione sui lati del riquadro dell'area [0, lato]. | Stessa regola dello Step 2 con le aree di F3.1. |
| XD-17 | 2026-11-20 | XC2: sigma_in, sigma_e, ell, ell2 raddoppiati; w, w2, D invariati. | «Scala dello spostamento raddoppiata» (REGOLA_B, G-6): si raddoppiano le lunghezze, non le probabilità né la profondità della zona di bordo. |
| XD-18 | 2026-11-20 | Strisce con `simulate_stripes()` del pacchetto sulla griglia vera dell'area: riga i = `row` meno la riga pari più piccola del riquadro, colonna j idem (la parità segue `array_row`, `array_col` di V); normalizzazione a media 1 sui quadratini con il centro dentro le regioni dell'area; fattore tagliato al tetto M*: c_t = min(c, M*). | Come `simulate_visium_hd_pieces()` (DEC-075). |
| XD-19 | 2026-11-20 | Cattura: ogni conteggio (quadratino, gene) arrivato diventa Binomiale(n, min(1, k c_t)), con k unico per tessuto (uguale per le 5 aree e i 2 semi). Fondo: Poisson(b_p c_t) per quadratino, con b_p = b di `capture_params(A)` per la quota dei 408 geni negli UMI delle simulazioni S di F3.1 dello stesso tessuto (`F3.1/main/sim/<area>_k0`, `_k1`, media delle 10); geni del fondo estratti con la composizione dei trascritti Xenium nati nella finestra dell'area. | Nel simulatore il fondo è b c per quadratino su tutti i geni, con la composizione dei trascritti generati nell'area: sui 408 geni la sua parte è b per la quota del pannello. |
| XD-20 | 2026-11-20 | k di XC (e separatamente di XC2) risolve: media sulle 5 aree e 2 semi degli UMI attesi per quadratino `in_reg` (trascritti catturati attesi più fondo atteso) = bersaglio di V di XD-11. Si calcola sugli arrivi dopo lo spostamento e sulle strisce realizzate, poi si estrae la cattura. Se k M* > 1, il taglio min(1, ·) entra nel calcolo di k e si riporta la quota di quadratini e di trascritti attesi tagliati. | «Un fattore di profondità unico come per X1»; nel simulatore la profondità moltiplica i trascritti generati e il fondo resta b c. |
| XD-21 | 2026-11-20 | Semi di XC: per area i (0...9) e seme k (0, 1): `SeedSequence([20261125, i, k])` per numpy (spostamento, cattura, fondo) e, per `simulate_stripes()`, il seme intero = primo valore di `generate_state(1)` di `SeedSequence([20261126, i, k])` modulo 2^31-1. XC e XC2 usano gli stessi semi (numeri casuali comuni: il confronto fra i due vede solo la scala). | Seme dichiarato; G-6 confronta XC2 con X1, e i numeri comuni con XC tolgono rumore al confronto XC contro XC2. |
| XD-22 | 2026-11-20 | Uscite: `/mnt/micron/geo_spatialtrans/riesame/B/X/<area>/`, `X1/<area>/`, `XC/<area>_k<seme>/`, `XC2/<area>_k<seme>/`, con X.parquet, bins.parquet, cells.parquet, genes.tsv, export.json nel formato di `tools/F3.3_export.py`. Nessun metodo (bin2cell, RCTD) e nessuna misura eseguiti qui. | Consegna. |

## Esito della regola XD-5 (2026-11-20, Parte B, prima di ogni braccio Xenium)

| id | data | scelta | motivo |
|---|---|---|---|
| XD-23 | 2026-11-20 | La regola XD-5 ferma i bracci X, X1, XC, XC2. Il file di 10x (`Xenium_P2_cell_info.csv.gz`, 1 105 615 cellule, riquadro 10 869 x 20 941 µm) e lo zip GEO di P2 (340 837 cellule, riquadro 6 655 x 7 270 µm) non sono lo stesso esperimento: 83 `cell_id` in comune, con centroidi distanti in mediana 8 020 µm; cellule GitHub entro 1 µm da una cellula GEO 0,63 %, cellule GEO entro 1 µm da una GitHub 2,0 %. La trasformazione di 10x e' un'affine esatta (residuo al 95° percentile 0,0008 µm, scala 1,0025) ma vale per l'altra sezione. Nessuna registrazione stimata da noi: sarebbe una scelta nuova presa dopo aver visto i dati. Si applica la regola: «se lo Xenium non e' usabile, la Parte B si fa con V, T, S, B sui geni completi e decide solo sul generatore». Integrita' dello zip verificata (15 874 246 655 byte, CRC32 dei 26 membri, md5 in `data_xenium/P2/inventory.tsv`). Uscite: `riesame/B/registration.json`, `results/riesame/xen_registration.csv`. | Regola scritta prima (XD-5). |

## Criterio «stessa acquisizione» (Luca, 2026-11-20, scritto prima del controllo)

Richiesta di Luca in chat: «controlla se il file GEO è un ritaglio del file grande [...] è la stessa acquisizione se almeno il 99% degli identificativi coincide e le coordinate differiscono di un solo vettore costante. Se è così, lo scarto è un cambio di origine, non una registrazione stimata: la domanda sulla cattura torna decidibile. [...] Non calcolare nessun braccio Xenium prima del mio sì.»

| id | data | scelta | motivo |
|---|---|---|---|
| XD-24 | 2026-11-20 | Stessa acquisizione se valgono tutte e tre: (1) almeno il 99 % dei `cell_id` di `cells.parquet` dello zip GEO compare in `Xenium_P2_cell_info.csv.gz`, dopo aver reso uniforme il formato (testo, spazi, eventuale suffisso); (2) per le cellule comuni lo scarto (x, y) GitHub meno GEO e' un vettore costante: deviazione standard di ciascuna componente sotto 1 µm; (3) per le cellule comuni i conteggi per cellula (`transcript_counts`, `total_counts`) coincidono in almeno il 99 %. Gli identificativi dei trascritti non si confrontano: il file di 10x ha solo colonne per cellula. Se vale: la trasformazione di 10x si applica allo zip GEO dopo il vettore (cambio di origine) e la domanda sulla cattura torna decidibile, ma nessun braccio Xenium si calcola prima del si' di Luca. Se non vale: resta XD-23. Uscita: `results/riesame/xen_same_acquisition.json`. | Criterio di Luca; la condizione (3) e' aggiunta dall'agente prima del controllo perche' gli identificativi dei trascritti non esistono nel file di 10x. Nota di trasparenza: il conteggio grezzo di 83 `cell_id` comuni era gia' nel registro di XD-23 quando il criterio e' stato scritto. |

## Criterio «quale registrazione per quale pacchetto» (Luca, 2026-11-20, scritto prima di scaricare e di guardare)

Richiesta di Luca in chat: «Scarica da GitHub Xenium_P1_cell_info.csv.gz e Xenium_P5_cell_info.csv.gz [...] applica il criterio XD-24 contro il nostro pacchetto GEO di P2 [...] Poi leggi solo l'indice e cells.parquet degli zip GEO di P1 e P5 [...] e confrontali con i tre file di registrazione. Scrivi il criterio e fai il commit prima di guardare. [...] Nessun braccio Xenium prima del mio sì.»

| id | data | scelta | motivo |
|---|---|---|---|
| XD-25 | 2026-11-20 | Per ogni coppia (pacchetto GEO P1, P2, P5) x (file di registrazione 10x P1, P2, P5), 9 coppie, si applica XD-24 senza modifiche (identificativi >= 99 %, scarto con deviazione standard < 1 µm per componente, conteggi uguali >= 99 %). Corrispondenza = coppia che soddisfa tutte e tre. Un pacchetto con due corrispondenze e' ambiguo e si riporta come tale, senza scegliere. I file GitHub si verificano con lo sha1 git dell'API GitHub; degli zip GEO di P1 e P5 si leggono con richieste HTTP parziali solo l'indice e `cells.parquet` (CRC32 del membro verificato), riportando i byte scaricati. Script `tools/riesame_xen_match.py`; uscite `results/riesame/xen_match.csv`, `results/riesame/xen_match.json`. Nessun braccio Xenium prima del si' di Luca. | Criterio di Luca, gia' scritto in XD-24. |
