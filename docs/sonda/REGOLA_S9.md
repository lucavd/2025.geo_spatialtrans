---
type: Rule
title: "S9, A4 con due catene mai usate in taratura (LuCA con RCTD; CellTypist), A2 sul banale, A5 col guadagno"
description: "Regola scritta prima del calcolo: catene U e C, tetti e soglie per catena, A4, A2 con incertezza per voce, A5 come guadagno sulla verità permutata, cancelli, controprove, esiti."
status: draft
generated: { by: geo-agent/progetto-master, at: 2026-11-27T12:45:00+01:00 }
sources:
  - id: paper
    path: docs/PAPER.md
  - id: s3
    path: docs/sonda/REGOLA_S3.md
  - id: s3c
    path: docs/sonda/REGOLA_S3c.md
  - id: s8
    path: docs/sonda/REGOLA_S8.md
  - id: rc
    path: docs/riesame/RIESAME_C.md
  - id: luca
    doi: 10.1016/j.ccell.2022.10.008
  - id: hlca
    doi: 10.1038/s41591-023-02327-2
---
# Perché

Ultima sessione della sonda (9 di 9). Il contratto[^paper] chiede una prima misura di ogni affermazione prima di G1. A2 e A4 non sono misurate. A4: «regge con una catena mai usata in taratura: D2 ≤ tetto della catena + δ». Il riesame[^rc] (§3) e STATO (minacce 1 e 2): i tipi del banale B vengono dalla stessa catena che misura (RCTD e classificatore con HLCA), e HLCA è tessuto sano, senza un tipo per le cellule tumorali.

Decisioni di Luca: 2026-11-26 in chat, «ok» al punteggio di A5 = guadagno sulla verità permutata, da pre-registrare qui e calcolare sui file di S8. 2026-11-27 in chat, per A4: «Tutte e due le catene» (LuCA con RCTD, e CellTypist).

# Già noto prima del commit

- A2 sulla catena di oggi: in `results/sonda/S3c/per_measure.csv` la voce più grande di T è 0,806 h (stroma, M-1); B la supera in 9 voci su 12 (due di poco: tumore M-2 0,83, tumore M-5 0,89). Nuovo: l'incertezza per voce.
- A5 col guadagno, esplorativo in S8: T5 1,000, S5c 0,400, A5 −0,600. Nuovo: errore standard, controprova.
- Catene U e C: nessun numero calcolato. LuCA scaricato prima del commit (md5 058d7069a9322dd4ac54d304ede2d2fc, 12 897 440 967 byte) e letti solo i metadati del riferimento (colonne, categorie); nessuna area toccata. Modello CellTypist letto (etichette, geni, md5).

# Bracci (invariati, sola lettura)

V = `S3_sonda/runs/V` (10 aree di H2), T = `S3b_sonda/runs/T` (10 aree omologhe di H1), B = `S3_sonda/runs/B` (10 aree × semi k0, k1), P = `S3_sonda/runs/P` (controllo positivo di S3), sotto `/mnt/micron/geo_spatialtrans/`. Si leggono le esportazioni (`bins.parquet`, `X.parquet`, `genes.tsv`, `cells*.parquet`, `blocks.csv`, `blocks_counts.mtx`, `export.json`). md5 di tutti i file letti prima e dopo. Uscite in `/mnt/micron/geo_spatialtrans/S9_sonda/{U,C}/` (un file per unità, rilanciabile); risultati in `results/sonda/S9/`.

# Catena U: stessa catena, riferimento LuCA

`tools/sonda/S3_chain.sh`, `tools/F3.1_b2c.py`, `tools/F3.1_rctd.R` (seme 20261118), `tools/sonda/S3_measures.py` invariati, con `S3_ROOT = S9_sonda/U` e `S3_REF = S9_sonda/U/ref`. Le esportazioni si collegano (symlink) nelle cartelle di U; le uscite della catena sono nuove.

Riferimento (`tools/sonda/S9_reference.py` = `S3_reference.py` con la sorgente come parametro):

- LuCA core atlas[^luca], CELLxGENE dataset 232f6a5a-a04c-4758-a6e8-88ab2e3a6e69, versione 46e0287b-9a33-4e83-99f3-8c044131bfdc; conteggi grezzi `raw/X` (interi).
- Etichetta (13 tipi): `ann_coarse` (12 tipi grossolani degli autori), con le cellule `cell_type_major == "Tumor cells"` tolte da «Epithelial cell» e messe nel tipo «Tumor cells». Motivo: grana vicina agli 11 tipi di HLCA `ann_level_2` (costo di RCTD simile) e un tipo per le cellule tumorali da un'annotazione curata (non `cell_type_predicted`, che è una previsione).
- Cellule: `tissue == "lung"`; saggi a UMI (fuori Smart-seq2); UMI ≥ 20 sui geni comuni; tipi con almeno 25 cellule; al più 3000 per tipo con `default_rng(20261123)`; marcatori di M-3 come S3 (50 per tipo).
- Geni: identificativi Ensembl comuni fra LuCA e le feature di H2.

# Catena C: CellTypist al posto del classificatore e di RCTD

Assegnazione e nuclei: bin2cell, uscite di U (identiche a quelle della catena di oggi: cancello G-B2C). Modello CellTypist 1.6.5 `Human_Lung_Atlas.pkl` (HLCA[^hlca], 61 etichette, 5017 geni, md5 bb198782a612415e56f6565fa9401d21). Ingresso: conteggi normalizzati a 10 000 e log1p, geni per nome di Space Ranger in comune con il modello. Etichette riportate a HLCA `ann_level_2` con la maggioranza delle cellule HLCA di quella `ann_finest_level` (tabella salvata). Nessun `majority_voting`. Codice: `tools/sonda/S9_celltypist.py`.

- Cellule: CellTypist sui conteggi del nucleo allargato (stessa matrice di `transfer_labels`), almeno 20 UMI; tipo = livello 2 dell'etichetta con probabilità massima.
- Blocchi da 8 µm (`blocks.csv`, `blocks_counts.mtx`, l'ingresso di RCTD): considerati i blocchi con almeno 20 UMI (UMI_min di RCTD). Tipi del blocco = tipi di livello 2 con probabilità ≥ 0,5 (`prob match`); «assegnato» se ne ha almeno uno (analogo di «non reject»); «primo» = livello 2 dell'etichetta con probabilità massima.
- M-1, M-2: come oggi (bin2cell). M-3: marcatori estranei con il tipo CellTypist della cellula, marcatori HLCA di S3. M-4c: fra i blocchi assegnati, quota con un solo tipo di livello 2 (analogo dei singoletti). M-5c: blocchi assegnati con un nucleo di almeno 4 quadratini: primo del blocco = tipo CellTypist della cellula del nucleo. M-6c: coppie di blocchi assegnati vicini (destra, sotto): stesso primo. Stesso codice di `F3_3_common.measures` per M-1, M-2, M-3.

# Metro per catena (invariato)

Per X ∈ {U, C}: tabella `area_values_X.csv` (`S3_table.py`), margini h_X con `S3b_h.py` (κ = 1, stessi semi) sulle 5 aree di V per tessuto, D con le funzioni di `tools/sonda/S3c_D.py` (stimatore appaiato, bootstrap per coppie stratificato, B = 2000, seme 20261123, metodo (ii) di S3c[^s3c], f = 1,584). M-4 (o M-4c) si toglie per un tessuto solo se la media di V è ≥ 0,99 (regola di S3). Tetto_X = D_X(T); pavimento_X = D_X(B); soglia_X = tetto_X + 1 h.

# A4

Per catena: e_X = D_X(B) − (D_X(T) + 1), errore standard = se bootstrap di D_X(B) − D_X(T) (stessi indici) × f. Esito per catena con la regola dello strumento: regge se e_X + 2 se ≤ 0; non regge se e_X − 2 se > 0; altrimenti incerta. In `results/misura.json` A4 = e_U (catena principale, misure definite come oggi); A4 «regge» solo se regge anche con C. Esiti diversi fra U e C: dubbio scritto per G1.

Controprova CP-chain (la catena distingue il banale dal vero): per catena, D_X(B) − D_X(T) «diverso» (regola di S3c: |Δ| − w > 1 h). Se fallisce, A4 con quella catena è «non misurabile» (livello: affermazione) e non entra nel quadro.

# A2 (catena di oggi, misura principale)

T_max = massimo sulle 12 voci del T_d puntuale di `per_measure.csv` (soglia fissa, come la soglia di A1). Per voce, se_v = deviazione standard bootstrap di B_d (stessi indici di S3c) × f: «oltre» se B_d − 2 se_v > T_max, «dentro» se B_d + 2 se_v ≤ T_max, altrimenti indeterminata. Valore di A2 = voci con B_d puntuale > T_max; se = deviazione standard bootstrap del conteggio × f; esito con la regola dello strumento (serve 0). Descrittivo: lo stesso per U e C.

# A5 col guadagno

Per unità e metodo: guadagno = `score` − `score_shuf` di `results/sonda/S8/per_unit.csv`. T5, S5c, A5 con le stesse divisioni, lo stesso bootstrap (B = 2000, `default_rng(20261201)`) e la stessa regola di esito di S8[^s8], sul guadagno invece del coseno. Codice: funzioni di `S8_A5.py`.

Controprova CP-gain (il guadagno non è la quantità di UMI): tau-b fra la classifica di X1 per guadagno e quella per UMI di G_p assegnati per cellula (`panel_*`, media sulle 10 aree). Se ≥ 0,8, la classifica per guadagno non si distingue da «chi assegna più UMI» e A5 non si legge come fedeltà.

# Cancelli (se uno fallisce i numeri della sua parte non si leggono)

| id | controllo | esito richiesto |
|---|---|---|
| G-MD5 | md5 dei file letti (esportazioni, H5AD, modello) prima e dopo | identici |
| G-REF | `S9_reference.py` con la sorgente HLCA riproduce i file di `S3_sonda/ref`; per LuCA conteggi interi, 13 tipi tenuti con il tipo tumorale, ognuno ≥ 25 cellule | md5 identici; 13/13 |
| G-B2C | uscite di bin2cell di U uguali a quelle di oggi per ogni unità | 50/50 identiche |
| G-CT | CellTypist: rieseguito su un'unità, uscite identiche; 61/61 etichette mappate; sulle cellule del riferimento HLCA di S3 accuratezza di livello 2 ≥ 0,85 (controllo positivo del classificatore sui suoi dati) | 3/3 |
| G-REV | sotto-agente senza storia ricalcola D_X(T), D_X(B), e_X, A2, A5 dalle tabelle | differenze < 1e-9 |

Controlli del metro per catena, come S3c: C-4 (P peggio di T, «diverso») e C-5 (B0 contro B1 non «diverso»). Se uno fallisce, la catena è «non misurabile».

# Attese (scritte prima)

- U: tetto 0,2-1,0 h; D_U(B) 2,5-3,5 h; e_U circa +1,5: A4 non regge. Motivo: il divario di B sta nella disposizione dei tipi (bilancio, 69-71 %), che un riferimento diverso non cambia.
- C: e_C > 0, A4 non regge; meno certo, perché M-4c e M-6c sono misure nuove.
- A2: 9 voci, non regge (già noto). A5 col guadagno: circa −0,6, non regge; CP-gain sotto 0,8.

# Esito e conseguenze

- Cancelli passati: A4, A2 e A5 entrano in `results/misura.json` tramite `tools/sonda/S9_to_misura.py` (copia, nessun ricalcolo); il check smoke la chiama dopo S8. Tetti e soglie per catena in `results/sonda/S9/`.
- Livello: affermazione. Nessun esito di S9 chiude l'obiettivo: decide G1, con le regole del contratto.

[^paper]: docs/PAPER.md
[^s3c]: docs/sonda/REGOLA_S3c.md
[^s8]: docs/sonda/REGOLA_S8.md
[^rc]: docs/riesame/RIESAME_C.md
[^luca]: Salcher et al., High-resolution single-cell atlas reveals diversity and plasticity of tissue-resident neutrophils in non-small cell lung cancer, Cancer Cell 2022, doi 10.1016/j.ccell.2022.10.008
[^hlca]: Sikkema et al., An integrated cell atlas of the lung in health and disease, Nature Medicine 2023, doi 10.1038/s41591-023-02327-2
