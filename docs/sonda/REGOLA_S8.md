---
type: Rule
title: "S8, A5: la classifica di cinque regole di assegnazione dei quadratini sul simulato banale e contro lo Xenium di X1"
description: "Regola scritta e committata prima di ogni calcolo di S8: metodi (ENACT 0.2.5 e bin2cell), verità, abbinamento, punteggio, concordanza, tetto riscritto senza X2, cancelli, controprova, esito."
status: stable
generated: { by: geo-agent/progetto-master, at: 2026-11-25T19:30:00+01:00 }
sources:
  - id: paper
    path: docs/PAPER.md
  - id: s5
    path: docs/sonda/REGOLA_S5.md
  - id: s3
    path: docs/sonda/REGOLA_S3.md
  - id: enact
    url: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11925495/
---
# Perché

A5 del contratto[^paper]: la classifica dei metodi sul simulato concorda con quella contro lo Xenium. Oggi il simulatore è il banale B[^s3], quindi questa è la prima misura, sul pavimento. X2 non è registrata: il tetto scritto (X1 contro X2) non esiste e si riscrive qui.

Decisioni di Luca, 2026-11-25 in chat: metodi = «Regole di ENACT più bin2cell»; margine = «Uno scambio di coppia: 0,2».

# Metodi (5), stessi nuclei, stesso raggio

Ingresso comune per area: griglia da 2 µm, conteggi grezzi su tutti i geni, `label` dei nuclei (B: nuclei veri del simulato; X1: StarDist sull'H&E di X1, esportazione XVfull di S5[^s5]).

- E-naive, E-area, E-gene, E-cluster: funzioni di ENACT 0.2.5[^enact] (`enact-SO`, `.venv-enact`) senza modifiche: nucleo = unione dei quadratini del nucleo, allargato di 2 quadratini (4 µm) e ritagliato nel Voronoi dei centri (`expand_nuclei_with_voronoi`); quadratini come quadrati da 2 µm; join spaziale `intersects`; regole `naive_assignment`, `weight_by_area_assignment`, `weight_by_gene_assignment`, `weight_by_cluster_assignment` (n_clusters 4, n_pcs 250, valori predefiniti); somma per cellula come `convert_adata_to_cell_by_gene`. Solo il codice di contorno (lettura dei nostri file) è nostro.
- b2c: bin2cell 0.3.5 `expand_labels(max_bin_distance = 2)` come nella catena (`tools/F3.1_b2c.py`), conteggi grezzi sommati per cellula.

# Bracci e verità

- B: 10 aree di H2 × 2 semi (k0, k1), esecuzioni di S3. Verità = conteggi veri per cellula: il territorio di ogni cellula si ricostruisce rieseguendo `S3_sim_B.py` con il suo seme fino ai territori; ogni UMI del simulato sta nel territorio della sua cellula, quindi conteggio vero = somma sui quadratini del territorio. Abbinamento = stessa etichetta.
- X1: le 10 aree omologhe di S5 (tutte nella metà di sviluppo). Verità = trascritti Xenium grezzi con qv ≥ 20 della cellula Xenium (`cell_id`), sui 284 geni del pannello G_p, senza diradamento. Abbinamento: per ogni nucleo StarDist, il nucleo Xenium (griglia `label` di XB, stessa griglia da 2 µm) con più quadratini in comune; tenuto se la scelta è reciproca e i quadratini in comune sono ≥ metà del nucleo più piccolo.

# Punteggio, classifica, concordanza

- Per cellula e metodo: coseno fra i conteggi assegnati e i conteggi veri sui 284 geni di G_p (in B come in X1, stesse condizioni). Cellule tenute: abbinate, con almeno un conteggio di G_p nella verità e in ognuno dei 5 metodi (insieme comune ai 5 metodi).
- Per area: media dei coseni. B: media dei due semi per area.
- Classifica: punteggio medio sulle aree. Concordanza = tau-b di Kendall fra i 5 punteggi di due bracci.
- Metà: tutte le divisioni delle 10 aree in due metà di 5 con tumore 3 + 2 e stroma 2 + 3 (100 divisioni, più le speculari: 200 coppie ordinate).
- Tetto riscritto: T5 = media sulle divisioni di tau(X1 su una metà, X1 sull'altra). È il sostituto di «X1 contro X2»: dati veri indipendenti contro dati veri.
- Simulato: S5c = media sulle stesse divisioni di tau(B su una metà, X1 sull'altra).
- A5 = S5c − T5 (verso alto), serve ≥ −0,2 (uno scambio di coppia su 10 con 5 metodi).
- Incertezza: B = 2000, `default_rng(20261201)`: a ogni replica una divisione a caso, poi ricampionamento con ripetizione delle aree dentro ciascuna metà e per tessuto (le metà restano disgiunte); se = deviazione standard delle repliche di A5. Stesso schema per T5 e S5c.
- Descrittivi: tau sulle 10 aree intere, punteggi per metodo, tessuto e braccio, quota di nuclei abbinati, UMI di G_p per cellula.

# Cancelli (se uno fallisce i numeri non si leggono)

| id | controllo | esito richiesto |
|---|---|---|
| G-ENACT | md5 dei file del pacchetto prima e dopo; per area: E-naive = somma dei quadratini non condivisi; E-area: pesi = area di sovrapposizione / 4 µm², ricalcolati con shapely a parte; E-cluster: totale = somma dei quadratini intersecati; E-gene: totale = somma dei quadratini intersecati meno i conteggi dei geni 0/0 contati a parte | 30/30 aree, tolleranza relativa 1e-9 |
| G-TRUTH | B: territori ricostruiti riproducono `label` esattamente; somma della verità = `umi_total` di `B_check.json` | 20/20 |
| G-MATCH | X1: nessuna cellula Xenium abbinata due volte; abbinamento identico se rieseguito | 10/10 |
| G-ID | una area di B e una di X1 rieseguite: punteggi identici bit per bit | 2/2 |
| G-REV | un sotto-agente senza storia ricalcola T5, S5c, A5 dalla tabella per area e i punteggi di 2 aree dai file | differenze < 1e-9 |

# Controprova

CP-shuf: in X1 i profili veri si permutano fra le cellule tenute della stessa area (`default_rng([20261201, i])`), poi si ricalcola T5. Attesa: T5 permutato con intervallo che contiene 0. Se resta lontano da 0, la classifica di X1 dipende da qualcosa che non è la verità (per esempio la quantità di UMI assegnata) e A5 non si legge come fedeltà.

# Esito (scritto prima)

1. Classifica vera non riproducibile: intervallo T5 ± 2 se contiene 0. A5 «non misurabile a questa profondità», registrata come tale. Livello: affermazione. A G1 si sceglie fra un punteggio più profondo (tipi cellulari) e il messaggio ristretto.
2. Altrimenti la regola dello strumento: regge se A5 − 2 se ≥ −0,2; non regge se A5 + 2 se < −0,2; altrimenti incerta.

Atteso: A5 ≈ −0,4, cioè non regge o incerta. Motivo: in B gli UMI sono uniformi nel territorio, nel vero si concentrano sui nuclei (voce 2 del bilancio). Le regole differiscono proprio sui quadratini di confine, il cui contenuto è diverso nei due casi. Nessun numero di letteratura: Crowell 2023 mostra solo che la fedeltà della classifica non segue il realismo.

# Uscite

`/mnt/micron/geo_spatialtrans/S8_sonda/{B,X1}/<unità>/` (un file per unità, rilanciabile); `results/sonda/S8/` (tabella per area, `A5.json`, md5). `check smoke` legge `A5.json` e lo scrive in `results/misura.json` (nessun ricalcolo).

[^paper]: docs/PAPER.md
[^s5]: docs/sonda/REGOLA_S5.md
[^s3]: docs/sonda/REGOLA_S3.md
[^enact]: ENACT: End-to-End Analysis of Visium High Definition (HD) Data, PMC11925495; codice github.com/Sanofi-Public/enact-pipeline

# Aggiunta prima del calcolo (2026-11-25, dopo la lettura del codice di ENACT 0.2.5)

Nessun numero di S8 era stato calcolato. Dettagli tecnici, nessuno cambia metodo, verità, punteggio o esito:

- Geni di ENACT: come `load_visiumhd_dataset` (use_hvg, n_hvg = 1000, scanpy `highly_variable_genes` su log1p dei conteggi normalizzati), con i 284 geni di G_p come marcatori: è il meccanismo di ENACT per i geni d'interesse; senza, i geni del pannello fuori dai 1000 sparirebbero dalle uscite. HVG calcolati per area (ENACT li calcola sul vetrino intero): limite dichiarato.
- Una sola tessera per area (`patch_size` 1e12): nessun confine di tessera dentro l'area.
- Identificativi dei quadratini non numerici («b<bin_id>»): ENACT rilegge i CSV intermedi e un identificativo numerico diventerebbe una posizione.
- Seme del generatore globale di numpy prima di ogni regola, 20261201 + 100 i + s (+ 50 per X1): la PCA di `weight_by_cluster` usa il solver randomizzato senza seme; KMeans ha già random_state = 0.
- bin2cell (b2c) ricalcolato in `S8_sonda` con `tools/sonda/S8_b2c.py` (stesse chiamate di `F3.1_b2c.py`), su tutti i quadratini dell'area come nella catena; ENACT sui quadratini in tissue, come sulla matrice filtrata.
- Cellule tenute anche con il criterio della catena: nucleo con almeno metà dei quadratini in `in_reg`.
- Codice: `tools/sonda/S8_A5.py`, `tools/sonda/S8_b2c.py`.
