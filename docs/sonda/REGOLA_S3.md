---
type: Rule
title: "S3, metro D sul polmone: aree, riferimento, braccio banale, controlli e attese"
description: "Regola scritta e committata prima di leggere l'espressione di H1 e H2: metà tenuta da parte, scelta delle aree di tumore e stroma, riferimento HLCA per RCTD, braccio banale tarato su H2, misura D con pavimento, tetto e soglia, controlli del metro, attese."
status: stable
generated: { by: geo-agent/progetto-master, at: 2026-11-23T09:30:00+01:00 }
sources:
  - id: paper
    path: docs/PAPER.md
  - id: regola_b
    path: docs/riesame/REGOLA_B.md
  - id: fatti_b
    path: docs/riesame/FATTI_B.md
  - id: s1
    path: docs/sonda/REGOLA_S1.md
---
# Perché

Il contratto (`docs/PAPER.md`[^paper]) chiede la misura principale D sul polmone: 12 voci (tumore e stroma per M-1...M-6 dopo bin2cell e RCTD con HLCA), unità l'area (5 per tessuto), bootstrap B = 2000. Il pavimento è D del braccio banale B tarato su H2; il tetto è D di H1 contro H2; la soglia è tetto + δ, δ = 1,0 h. Oggi il simulatore coincide con la catena banale, quindi D(B) è anche il valore di oggi di A1. Questa sessione costruisce metro e catena banale e misura i tre numeri. Deroga: S3 aperta nella stessa chat di S2 (Luca, 2026-11-21).

Stato al momento del commit: H1 e H2 scaricati (`tools/sonda/S3_download.sh`), solo dimensioni e md5 controllati; nessuna espressione di H1 o H2 letta. Letti finora: elenco dei tipi e dei tessuti di HLCA core.

# Dati

- H2 = `Visium_HD_Human_Lung_Cancer_HD_Only_Experiment2` è V, la sezione di taratura. H1 = `..._Experiment1` è T, il tetto. Stesso blocco FFPE di X1 e X2.
- Conteggi: `square_002um/raw_feature_bc_matrix.h5` (tutti i quadratini, come F3.1 sul colon), geni = feature della matrice; nomi dei geni come in Space Ranger (resi unici se servisse). Posizioni: `square_002um/spatial/tissue_positions.parquet` (in_tissue, array_row, array_col, pixel a piena risoluzione); µm per pixel da `square_008um/spatial/scalefactors_json.json`.
- Metà tenuta da parte: quella di `docs/sonda/REGOLA_S1.md`[^s1], invariata. Ogni script lavora solo su quadratini con centro x < m − 50 µm e lo verifica con un'asserzione registrata nel suo JSON (x massimo usato). La matrice si legge solo per le colonne dei quadratini delle aree (lettura per intervalli di `indptr`): l'espressione della metà tenuta da parte non si carica.

# Aree (5 di tumore e 5 di stroma per sezione, regola deterministica)

- Coordinate: x = (array_col + 0,5) × 8 µm, y = (array_row + 0,5) × 8 µm sui quadratini da 8 µm (`square_008um`). Candidati: quadrati di lato 1000 µm (come il colon), angolo su un reticolo da 100 µm, interamente nella metà di sviluppo (x massimo del quadrato < m − 50 µm).
- Copertura = quota dei quadratini da 8 µm del quadrato con in_tissue = 1 (nessuna espressione). Candidati ammessi: copertura ≥ 0,9; se i disgiunti non bastano per 10 aree, 0,8; poi lato 750 µm con 0,9 e 0,8. Ogni rilassamento si annota.
- Punteggio (espressione, `square_008um/filtered_feature_bc_matrix.h5`, quadratini in tissue del quadrato): e = quota di UMI sui geni epiteliali E = {EPCAM, KRT8, KRT18, KRT19, KRT7, CDH1}; s = quota sui geni stromali S = {COL1A1, COL1A2, COL3A1, DCN, LUM, COL6A3}; r = log2((e + 1e-4) / (s + 1e-4)). Un gene assente dal pannello di sonde si toglie e si annota.
- Tumore: si sceglie il candidato con r massimo, si tolgono i candidati che lo intersecano, cinque volte. Stroma: fra i candidati disgiunti da tutte le aree già scelte, r minimo, cinque volte. Parità: x minore, poi y minore. Numerazione 1...5 nell'ordine di scelta.
- Un'area di tumore con e ≤ s, o di stroma con s ≤ e, resta e si segnala.
- Stessa regola su H2 (V) e H1 (T). Uscita: `results/sonda/S3/areas_H1.csv`, `areas_H2.csv`.

# Catena (uguale per tutti i bracci)

1. Esportazione nel formato di F3.1 (bins, X, genes, cells). Nuclei: StarDist `2D_versatile_he` sull'H&E a piena risoluzione (`tissue_image.btf`) ritagliato sull'area più 50 pixel, parametri di F3.3 (normalizzazione 1 e 99,8, `predict_instances_big` a blocchi di 2048, seme 0, determinismo delle operazioni attivo); ogni quadratino da 2 µm prende l'etichetta del nucleo sotto il pixel del suo centro. Tutti i nuclei, senza filtri. in_reg = in_tissue e centro a ≥ 25 µm dal bordo del quadrato.
2. bin2cell (`tools/F3.1_b2c.py`, invariato) e RCTD in modalità doppietto (`tools/F3.1_rctd.R`, invariato, seme 20261118) con il riferimento qui sotto (sigla del riferimento `L`, uguale per tumore e stroma).
3. Misure M-1...M-6 con le etichette trasferite (funzioni di `tools/F3_3_common.py`, invariate, con il riferimento del polmone).

# Riferimento RCTD (HLCA)

- HLCA core, CELLxGENE, dataset 066943a2-fdac-4b29-b348-40cede398e4e (Sikkema et al. 2023, doi 10.1038/s41591-023-02327-2), md5 3d9352bec0733c42ef0ad53b125faba8; conteggi grezzi `raw/X`.
- Etichetta: `ann_level_2` (11 tipi). Motivo: grana vicina a quella del colon (pochi tipi per tessuto) e costo di RCTD contenuto; il livello 3 ha 25 tipi e una classe «None».
- Cellule: tessuto diverso da «nose»; UMI ≥ 20 sui geni comuni con il pannello di sonde; al più 3000 per tipo (estrazione con `numpy.random.default_rng(20261123)`); tipi con meno di 25 cellule esclusi.
- Geni: identificativi Ensembl comuni fra HLCA e le feature di H2, con i nomi di Space Ranger.
- Marcatori per M-3: regola di `tools/F3.1_reference.py` (50 geni per tipo per log2 del rapporto fra CP10k medio nel tipo e negli altri tipi, CP10k medio nel tipo ≥ 1).
- Limite dichiarato ora: HLCA è tessuto sano, le cellule tumorali non vi hanno un tipo; vale per tutti i bracci allo stesso modo.

# Braccio banale B (taratura su H2, metà di sviluppo)

Parametri per tessuto, sommando numeratori e denominatori sulle 5 aree di V dopo la catena:

- densità λ = nuclei della cellula (nucleo con almeno metà dei quadratini in in_reg, definizione di F3.3) / area in_reg;
- area del nucleo = numero medio di quadratini per nucleo, n̄; raggio r = √(4 n̄ / π) µm;
- composizione = quote dei tipi trasferiti fra quelle cellule (tipo non nullo);
- profilo di ogni tipo = somma dei conteggi dei nuclei allargati di bin2cell delle cellule di quel tipo, normalizzata;
- profondità ρ = UMI nei quadratini in_reg / numero di quadratini in_reg.

Simulazione su ciascuna delle 10 aree di V (stessa griglia, stessi in_tissue e in_reg), semi k = 0 e 1, generatore `numpy.random.default_rng([20261125, i, k])`, i = 0...9 nell'ordine tumore 1...5, stroma 1...5:

1. N ~ Poisson(λ · 4 µm² · numero di quadratini in tissue); centri = quadratini in tissue estratti con ripetizione, più uno scarto uniforme di ±1 µm.
2. Territorio = quadratini in tissue assegnati al centro più vicino (Voronoi discreto sulla griglia da 2 µm); cellule senza quadratini tolte e contate.
3. Nucleo = quadratini del territorio a distanza ≤ r dal centro.
4. Tipo estratto dalla composizione, indipendente.
5. Conteggi: n ~ Poisson(ρ · quadratini del territorio); geni ~ Multinomiale(n, profilo del tipo); ogni trascritto in un quadratino del territorio estratto a caso.
6. Fuori tessuto zero conteggi; nessuno spostamento, nessuna striscia, nessun fondo.

Differenza dal colon, con motivo: là B usava la cattura dello Step 2 (strisce, fondo) tarata sul colon. Qui ogni parametro viene da H2 e B non contiene pezzi del simulatore: cattura, spostamento e fondo sono pezzi che dovranno battere B. Controlli per unità (JSON): λ, composizione, ρ e n̄ realizzati contro i bersagli.

# Misura

- D(braccio) = radice con segno della media sulle celle tessuto × misura di (Δ² − Var Δ)/h², Δ = media del braccio − media di V; stimatore, deviazione combinata (V e T), bootstrap parametrico per aree (B = 2000, seme 20261123) e confronti (w = semi-ampiezza al 95 %) di `tools/riesame_D.py`, invariati. Bracci con due semi: prima la media dei semi per area. T = H1, un solo insieme.
- h per tessuto e misura: costruzione di `tools/riesame_B_h.py` sulle 5 aree di V (bootstrap della media, B = 2000, `default_rng([20261115, i, j])`, κ e fattore di asimmetria di R7). Il polmone non ha κ suoi: si usano quelli del colon, tumore = A2 (i = 1), stroma = A3 (i = 2). h di M-4 almeno 0,01.
- M-4 si toglie per un tessuto se la media dei singoletti di V è ≥ 0,99, deciso sulla tabella di V prima di calcolare D di ogni braccio.
- Gruppi: tumore + stroma (principale), tumore, stroma.
- Pavimento = D(B), tetto = D(T), soglia = tetto + 1,0 h; sensibilità con δ = 0,5 e 2. A1 oggi = D(B), errore standard = deviazione standard della sua distribuzione bootstrap.

# Controlli del metro (prima di confrontare i bracci)

| id | controllo | esito richiesto | se fallisce |
|---|---|---|---|
| C-1 | G-3 di REGOLA_B: spostamenti costruiti 0, 1, 2 h con numerosità e deviazioni osservate (V, H1), 200 repliche, seme 20261124 | quadrato medio a 0 entro ±0,15; 1 h stimato entro 1 ± 0,25; 2 h entro 2 ± 0,25 | stimatore distorto: D non si legge |
| C-2 | ricalcolo, codice del metro: il percorso di calcolo di S3 sulla tabella del colon (`results/riesame/B_area_values.csv`, geni comuni) ritrova D di T, S e B di `results/riesame/B_D_common.json` | entro 1e-9 relativo | errore nel metro: si corregge prima di tutto |
| C-3 | ricalcolo, catena: la catena ridotta del check (ritaglio di 250 µm di un'area di V e un B sullo stesso ritaglio) rifatta da capo dà misure identiche | entro 1e-9 | catena non riproducibile: si dichiara e si corregge |
| C-4 | controllo positivo: braccio P = aree di V con i vettori di conteggio dei quadratini in tissue permutati a caso fra i quadratini in tissue (`default_rng([20261126, i])`), nuclei fermi; la disposizione spaziale dell'espressione è distrutta | D(P) − D(T) «diverso», P peggio | il metro non vede un difetto grossolano: D non si usa per decidere |
| C-5 | controllo negativo: D(B, seme 0) − D(B, seme 1), ciascuno con le sue 5 aree per tessuto | non «diverso» (\|Δ\| − w ≤ δ) | il rumore è sottostimato: nessun confronto di B si legge |

# Attese (scritte prima)

- Tetto D(T) circa 2 h (colon 2,23).
- Pavimento D(B) circa 7-8 h (colon 7,77); quindi A1 oggi non regge (D(B) oltre la soglia, attesa circa 3 h).
- Precondizione di G1: pavimento − tetto > 1 h. Se D(B) − D(T) ≤ δ = 1 h, il contratto dice chiudere: lo si scrive e si porta a G1.
- D(P): oltre D(T) di più di δ; nessuna attesa numerica.

# Calcolo

Al più 20 processi a un thread in tutto (OMP, OPENBLAS, MKL, NUMBA a 1), numero letto dal file `WORKERS` del run; uscite in `/mnt/micron/geo_spatialtrans/S3_sonda/`, un file per unità, rilanciabile. Un pilota su un'area di V misura il tempo: se la corsa completa (V, T, B, P: 50 unità) supera circa 4 ore con 20 processi, ci si ferma e si chiede.

# Esito e conseguenze

- I controlli C-1, C-2, C-4, C-5 passano: pavimento, tetto e soglia entrano in `progetto.json` e A1 = D(B) in `results/misura.json`.
- Un controllo fallisce: il numero relativo si riporta come provvisorio, non entra in `progetto.json`, e la sessione lo dichiara.
- Bilancio dell'errore, prima lettura: per voce, termine di B meno termine di T (unità h²), dalla più grande.

[^paper]: docs/PAPER.md
[^s1]: docs/sonda/REGOLA_S1.md
