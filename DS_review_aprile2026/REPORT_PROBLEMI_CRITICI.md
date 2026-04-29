# Problemi Critici e Proposta di Nuova Architettura per il Framework di Simulazione Visium HD

**Data**: Marzo 2025  
**Progetto**: Simulazione Spatial Transcriptomics (Visium HD)

---

## 1. Il Problema Centrale

Il framework attuale genera profili trascrizionali **direttamente per ogni spot** della griglia, trattando ciascuno come entità indipendente. In realtà, nel Visium HD gli spot sono quadrati di 2×2 µm che coprono un tessuto fatto di **cellule** (diametro 15–30 µm). Più spot coprono la stessa cellula, e gli spot ai confini cellulari catturano RNA da cellule diverse.

Manca dunque il livello intermedio tra tessuto e griglia: **le cellule**.

La proposta è riorganizzare la simulazione in **due step sequenziali**, distinguendo chiaramente ciò che **esiste già** nel codebase da ciò che **va implementato ex novo**, e identificando i parametri che potranno essere definiti solo tramite **reverse engineering da dati reali Visium HD**.

---

## 2. Step 1 — Tessuto e Cellule

Questo step parte da un'immagine (reale o sintetica) e produce una **mappa di cellule** individuali con posizione, dimensione, tipo e struttura interna.

### 2.1 Cosa GIÀ ESISTE nel codebase

| Componente | File | Cosa fa |
|---|---|---|
| Caricamento immagine | `03_image_processing.R` | Carica PNG, converte in scala di grigi, applica soglia per identificare regioni di tessuto |
| Generazione tessuto sintetico | `03b_generate_synthetic_tissue.R` | Genera immagini sintetiche con blob gaussiani, patch Voronoi, noise frattale (3 livelli di complessità) |
| Clustering regioni | `04_clustering.R` | K-means spaziale che combina intensità pixel e posizione → assegna label di "tipo cellulare" a ciascun pixel |

Questi componenti coprono la parte **macro**: definire regioni tissutali e assegnare tipi cellulari alle zone dell'immagine. Funzionano e possono essere riutilizzati quasi invariati.

### 2.2 Cosa VA IMPLEMENTATO

#### A. Piazzamento cellule dentro le regioni [NUOVO]

Una volta definite le regioni tissutali dal clustering, bisogna **popolarle con cellule individuali**. Ogni cellula è un'entità con:

- **Centro (x, y)** in µm — posizione nel tessuto
- **Tipo cellulare** — ereditato dalla regione in cui cade, con possibilità di eterogeneità (non tutte le cellule in una regione sono dello stesso tipo)
- **Struttura interna**: nucleo (zona centrale, raggio ~3–5 µm) e citoplasma (corona esterna)
- **Forma e confini**: definiti dall'algoritmo di piazzamento (vedi sotto)

##### Forma delle cellule e packing

La tabella `cell_df` contiene le proprietà di ciascuna cellula, ma **la forma e il grado di addossamento** non sono attributi della tabella: sono il **risultato dell'algoritmo di piazzamento**. L'approccio raccomandato è la **tassellazione di Voronoi con vincoli di densità**:

1. Per ogni regione tissutale, generare N centroidi casuali (dove N = densità × area regione)
2. Applicare tassellazione di Voronoi → ogni centroide genera un **poligono irregolare** che rappresenta una cellula
3. I poligoni riempiono tutto lo spazio senza gap (cellule naturalmente addossate)
4. Il "raggio" equivalente si calcola dall'area del poligono: `radius = sqrt(area/π)`
5. Il nucleo è un cerchio centrato nel centroide del poligono
6. Lo spazio extracellulare si modella come zone dell'immagine senza centroidi (intensità sotto soglia)

**Perché Voronoi è l'approccio giusto:**

- Le cellule epiteliali sono realmente addossate, e Voronoi dà poligoni piccoli e compatti in regioni dense
- Le cellule stromali sono più rade, e Voronoi dà poligoni grandi (parte del poligono è matrice extracellulare)
- I confini tra poligoni adiacenti sono il luogo naturale dove gli spot catturano RNA da due cellule
- La forma irregolare dei poligoni è più realistica di cerchi perfetti

```
CONCETTO: Voronoi con densità variabile

  ┌───────────── Regione "epitelio" (denso) ──────────────┐
  │ ┌──┬──┬──┐                                            │
  │ │●1│●2│●3│  poligoni piccoli, addossati               │
  │ ├──┼──┼──┤  ~4000-6000 cellule/mm²                   │
  │ │●4│●5│●6│  ogni ● = centroide, box = poligono Voronoi│
  │ └──┴──┴──┘                                            │
  ├───────────── Regione "stroma" (rado) ─────────────────┤
  │ ┌──────┬────────┐                                     │
  │ │  ○1  │   ○2   │  poligoni grandi, pochi             │
  │ │      │        │  ~1000-2000 cellule/mm²             │
  │ └──────┴────────┘                                     │
  ├───────────── Confine regioni ─────────────────────────┤
  │ ┌──┬──┬────┬────┐  mix di tipi al confine             │
  │ │●7│●8│ ○3 │ ○4 │                                    │
  │ └──┴──┴────┴────┘                                     │
  └───────────────────────────────────────────────────────┘
```

La **densità cellulare** varia per tipo di regione:
- Epitelio: denso (~4000–6000 cellule/mm²) → poligoni piccoli
- Stroma: rado (~1000–2000 cellule/mm²) → poligoni grandi
- Regioni necrotiche: pochissime cellule intatte
- Spazio extracellulare: nessun centroide (lumi vascolari, matrice)

#### B. Risoluzione dell'immagine di input [CHIARIMENTO]

L'immagine **non deve avere** risoluzione cellulare. Serve solo a definire le regioni macro del tessuto. Anche un'immagine a bassa risoluzione (200×200 pixel, dove 1 pixel = ~30 µm) è sufficiente per distinguere stroma da epitelio.

Le cellule vengono poi piazzate **computazionalmente** all'interno delle regioni, a precisione micrometrica, indipendentemente dalla risoluzione dell'immagine di partenza.

Se in futuro si volessero usare le posizioni reali dei nuclei (da un'immagine H&E ad alta risoluzione), si potrebbe usare un tool di segmentazione (StarDist, Cellpose) per estrarre i centroidi — ma questo è un'estensione futura, non necessaria ora.

### 2.3 Output dello Step 1

Due strutture dati complementari:

**1. Tabella `cell_df`** — una riga per cellula (proprietà scalari):

| cell_id | x | y | eq_radius | cell_type | nucleus_r | region_id | area_um2 |
|---|---|---|---|---|---|---|---|
| 1 | 102.3 | 55.7 | 10.0 | epithelial | 3.5 | region_1 | 314 |
| 2 | 115.1 | 58.2 | 9.5 | epithelial | 3.2 | region_1 | 283 |
| 3 | 245.0 | 300.4 | 15.0 | stromal | 4.5 | region_3 | 707 |
| ... | ... | ... | ... | ... | ... | ... | ... |

Dove `eq_radius = sqrt(area/π)` è il raggio equivalente calcolato dall'area del poligono Voronoi.

**2. Lista `cell_polygons`** — i poligoni Voronoi (geometria):

Ogni cellula ha un poligono associato (lista di vertici o oggetto `sf`/`sp`) che definisce i confini esatti della cellula. Questi poligoni vengono usati nello Step 2 per calcolare le frazioni di overlap spot–cellula.

```r
# Esempio struttura
cell_polygons[[1]]  # poligono Voronoi della cellula 1
# POLYGON((97.1 50.2, 107.5 50.2, 107.5 61.2, 97.1 61.2, 97.1 50.2))
```

---

## 3. Step 2 — Griglia Visium HD e Profili Trascrizionali

Questo step sovrappone la griglia di spot Visium HD al tessuto cellulare, genera i profili trascrizionali per ogni cellula, e poi aggrega l'espressione a livello di spot.

### 3.1 Cosa GIÀ ESISTE nel codebase

| Componente | File | Cosa fa | Riutilizzabile? |
|---|---|---|---|
| Creazione griglia | `05_grid_sampling.R` | Genera griglia regolare con risoluzione e dimensioni configurabili | ✅ Sì, adattando i parametri a 2/4/8 µm |
| Parametri espressione | `06a_expression_params.R` | Inizializza e valida parametri per la generazione | ✅ Sì, invariato |
| Baseline expression | `06b_expression_baseline.R` | Profili baseline per tipo cellulare con marker genes | ✅ Sì, invariato |
| Dispersione NB | `06d_dispersion_params.R` | Parametri dispersione spazialmente variabili | ✅ Sì, invariato |
| Library size | `06e_library_size.R` | Generazione library size con effetti spaziali | ⚠️ Da adattare: i valori vanno ricalibrati per la risoluzione spot (molto più bassi a 8µm) |
| Dropout | `06f_dropout_models.R` | Dropout expression-dependent + ambient RNA | ⚠️ Da adattare: parametri diversi per spot piccoli |
| Moduli genici | `06g_gene_modules.R` | Co-espressione e network regolatori | ✅ Sì, invariato |
| Correlazione spaziale (GRF) | `06h_spatial_correlation.R` | Gaussian Random Field | ⚠️ Da modificare: serve un campo GRF diverso per modulo/gene, non uno solo per tutti |
| Generazione conteggi | `06j_expression_generation.R` | Sampling NB + scaling + validazione | ⚠️ Da adattare: va applicato a livello cellula, non spot |

**In sintesi**: la maggior parte dei moduli statistici (NB, dispersione, dropout, moduli genici) è solida e riutilizzabile. Va cambiato il **livello a cui operano** (cellula anziché spot) e vanno ricalibrati i parametri per la risoluzione Visium HD.

### 3.2 Cosa VA IMPLEMENTATO

#### A. Mapping spot → cellula [NUOVO]

Per ogni spot della griglia, calcolare **quali cellule copre** e in che **percentuale**:

```
CONCETTO: Spot che coprono una singola cellula (cerchio grande)
         e spot al confine tra due cellule

         Cellula A (tipo epiteliale)    Cellula B (tipo stromale)
         raggio 10µm                    raggio 15µm
              ┌──┐
              │  │ ← spot 2×2 µm, 100% dentro cellula A
         ┌────────────┐          ┌──────────────────┐
         │    ┌──┐    │          │                  │
         │    │Nu│cleo│          │                  │
         │    └──┘    │          │                  │
         │         ┌──┼──┐       │                  │
         │         │60│40│ ← spot al confine: 60% A, 40% B
         └─────────┼──┼──┘       │                  │
                   └──┘          │                  │
                                 └──────────────────┘
```

Il risultato è una tabella **many-to-many**:

| spot_id | cell_id | fraction | in_nucleus |
|---|---|---|---|
| spot_1001 | cell_42 | 1.00 | TRUE |
| spot_1002 | cell_42 | 1.00 | FALSE |
| spot_1003 | cell_42 | 0.60 | FALSE |
| spot_1003 | cell_43 | 0.40 | FALSE |
| spot_1004 | — | 0.00 | — |

Dove:
- `fraction` = percentuale dell'area dello spot coperta dalla cellula
- `in_nucleus` = lo spot cade nella zona del nucleo (utile per variazione intra-cellulare)
- `spot_1004` non copre nessuna cellula → spazio extracellulare → riceverà solo ambient RNA

Questa tabella **è il ground truth per la deconvoluzione**. Chi sviluppa algoritmi di deconvoluzione dovrà ricostruire queste frazioni partendo solo dai conteggi degli spot.

#### B. Generazione profili trascrizionali per CELLULA [ADATTAMENTO]

I profili vanno generati a livello di **cellula** (non di spot):

1. Ogni cellula riceve un profilo baseline dal suo tipo cellulare (`06b` — esiste)
2. Si aggiungono effetti spaziali tramite GRF (`06h` — esiste, da modificare per generare campi multipli)
3. Si aggiungono moduli genici co-espressi (`06g` — esiste)
4. Si campiona con NB e dispersione cellula-specifica (`06d`, `06j` — esistono)
5. Si applica library size per cellula (`06e` — esiste, da riscalare)

**Output**: matrice `cell_counts` (geni × cellule) — questo è il **ground truth cellulare**.

#### C. Aggregazione cellula → spot [NUOVO — Cuore del realismo]

Per ogni spot della griglia, l'espressione osservata è una **miscela pesata** delle cellule che copre:

```
spot_counts[spot, gene] = Σ_cellule ( fraction[spot, cellula] × cell_counts[cellula, gene] )
                        + ambient_rna[spot, gene]
                        + noise_tecnico_cattura
```

Tre scenari:

| Scenario | Esempio | Cosa produce |
|---|---|---|
| **Spot puro** | Spot interamente dentro una cellula | Profilo di quella cellula + variazione intra-cellulare |
| **Spot misto** | Spot a cavallo tra cellula A e cellula B | Mix: 60% profilo A + 40% profilo B → **target deconvoluzione** |
| **Spot vuoto** | Spot su spazio extracellulare | Solo ambient RNA (pochissimi UMI) |

#### D. Variazione intra-cellulare: nucleo vs membrana [NUOVO — Da calibrare]

Gli spot sulla stessa cellula non avranno profili identici. Oltre al noise stocastico di cattura, la **posizione dello spot dentro la cellula** influenza quali trascritti vengono catturati:

- **Spot sul nucleo**: arricchiti in pre-mRNA, trascritti nascenti, RNA nucleari (lncRNA, alcuni snoRNA)
- **Spot sul citoplasma**: arricchiti in mRNA maturo in traduzione attiva
- **Spot alla periferia/membrana**: arricchiti in mRNA localizzato (es. trascritti per proteine di membrana, trascritti trasportati attivamente)

#### E. Binning 2µm → 8µm [NUOVO — Semplice da implementare]

Una volta generati i conteggi a livello di spot 2×2 µm, il binning è una semplice somma:

```
bin_8um[i,j] = Σ spot_2um[4i+a, 4j+b]  per a,b ∈ {0,1,2,3}
```

Ogni bin 8×8 µm somma 16 spot nativi. Questo riduce sparsità e aumenta la library size per bin, esattamente come fa il software Space Ranger di 10x.

---

## 4. Cosa Richiede Reverse Engineering da Dati Reali

Alcuni parametri del modello **non possono essere definiti teoricamente** e richiedono analisi di dati reali Visium HD per essere calibrati. Questa è una fase fondamentale che dovrà precedere (o accompagnare) il completamento del simulatore.

### 4.1 Variazione intra-cellulare

**Domanda**: Quanto variano i profili trascrizionali tra spot che coprono la stessa cellula?

**Come rispondere con dati reali:**
1. Prendere un dataset Visium HD reale con H&E ad alta risoluzione
2. Usare un tool di segmentazione nucleare (StarDist, Cellpose) sull'H&E per ottenere i confini cellulari
3. Mappare gli spot (o bin 8µm) alle cellule segmentate
4. Per ogni cellula coperta da ≥ 3 spot: calcolare la **correlazione media** tra i profili dei suoi spot
5. Confrontare con la correlazione tra spot di cellule diverse dello stesso tipo

**Metriche attese:**
- Correlazione intra-cellulare: alta (ρ ~ 0.7–0.95?)
- Correlazione inter-cellulare stesso tipo: media (ρ ~ 0.3–0.6?)
- Correlazione inter-cellulare tipo diverso: bassa (ρ ~ 0.0–0.2?)

**Parametro da calibrare**: `intra_cell_variation` — la deviazione standard del noise aggiunto tra spot della stessa cellula.

### 4.2 Effetto nucleo vs citoplasma

**Domanda**: Quali geni sono arricchiti nel nucleo vs nel citoplasma, e di quanto?

**Come rispondere con dati reali:**
1. Sullo stesso dataset, classificare gli spot come "nucleari" o "citoplasmatici" in base alla segmentazione H&E
2. Per ogni gene, confrontare l'espressione media tra spot nucleari e citoplasmatici della stessa cellula
3. Identificare geni con bias significativo (test di Wilcoxon o simile)

**Risultato atteso**: una lista di geni con il loro **rapporto nucleo/citoplasma** (nuclear enrichment score).

**Nota**: in prima approssimazione, questo effetto può essere trascurato e modellato semplicemente come noise aggiuntivo intra-cellulare. Il modello nucleo/citoplasma è un raffinamento per una fase successiva.

### 4.3 Profilo degli spot di confine (deconvoluzione ground truth)

**Domanda**: Come si comportano realmente gli spot a cavallo tra due cellule?

**Come rispondere con dati reali:**
1. Identificare spot che cadono sul confine tra due cellule (segmentazione H&E)
2. Stimare le frazioni di overlap geometrico
3. Verificare se il profilo trascrizionale dello spot è effettivamente un **mix lineare** dei profili delle due cellule, o se ci sono deviazioni (es. effetti di diffusione RNA, permeabilizzazione asimmetrica)

**Parametri da calibrare:**
- Il mixing è lineare? O ci sono effetti non-lineari?
- C'è diffusione laterale di RNA (spot che "vedono" RNA da cellule che non toccano geometricamente)?
- Quanto è ampio l'effetto di diffusione (in µm)?

### 4.4 Library size a risoluzione 8µm

**Domanda**: Quanti UMI cattura un tipico bin 8×8 µm?

**Come rispondere con dati reali:**
1. Prendere un dataset Visium HD pubblico a risoluzione 8µm
2. Calcolare la distribuzione di UMI per bin
3. Stratificare per tipo di tessuto (denso vs rado)

**Valori attesi** (da letteratura preliminare):
- Bin su tessuto denso: ~500–3000 UMI
- Bin su tessuto rado: ~50–500 UMI
- Bin su spazio extracellulare: ~0–50 UMI (solo ambient)
- Il framework attuale usa 8000 UMI/spot che è realistico per Visium classico (55µm) ma troppo alto per Visium HD (8µm)

### 4.5 Sparsità a risoluzione 8µm

**Domanda**: Quale percentuale di zeri ci si aspetta nei bin 8µm?

**Valore atteso**: >95% di zeri (molto più alto del ~89% attualmente simulato, che è realistico per Visium classico ma non per HD a 8µm).

---

## 5. Riepilogo: Esiste vs Da Implementare vs Da Calibrare

### ✅ ESISTE e funziona

- Caricamento e preprocessing immagini
- Generazione tessuto sintetico (regioni macro)
- Clustering regioni in tipi cellulari
- Creazione griglia di spot (da riparametrizzare per 2/8 µm)
- Modelli statistici di espressione: NB, dispersione variabile, dropout, moduli genici co-espressi, GRF
- Infrastruttura: chunking, parallelizzazione, matrice sparsa, validazione

### 🔨 DA IMPLEMENTARE

- **Piazzamento cellule** dentro le regioni (Voronoi/circle packing, con nucleo e citoplasma)
- **Mapping spot → cellula** (calcolo frazioni di overlap geometrico)
- **Aggregazione cellula → spot** (mix pesato + noise)
- **GRF multipli** (un campo per modulo genico, non uno solo per tutti)
- **Binning 2µm → 8µm** (somma di 4×4 spot)
- **Ricalibrazione parametri** per risoluzione HD (library size, dropout, sparsità)

### 🔬 DA CALIBRARE CON DATI REALI (reverse engineering)

- **Variazione intra-cellulare** → correlazione tra spot della stessa cellula
- **Effetto nucleo vs citoplasma** → enrichment score per gene
- **Linearità del mixing ai confini** → validazione del modello di deconvoluzione
- **Diffusione laterale RNA** → range in µm
- **Library size e sparsità** a risoluzione 8µm → distribuzione empirica
- **Rapporto spot puri / spot misti / spot vuoti** → dipende da densità cellulare e dimensione bin

Questi parametri definiscono la **differenza tra una simulazione plausibile e una simulazione realistica**. La prima può essere costruita con stime ragionevoli dalla letteratura; la seconda richiede confronto diretto con dati Visium HD reali.

---

## 6. Ordine di Implementazione Suggerito

```
FASE 1 — Architettura base (settimane 1-3)
├── 1a. Implementare cell layer (Voronoi/packing)
├── 1b. Implementare spot-cell mapping
├── 1c. Adattare generazione espressione a livello cellula
├── 1d. Implementare aggregazione cellula → spot
└── 1e. Implementare binning 2→8µm
    Output: pipeline funzionante con parametri stimati

FASE 2 — Reverse engineering (settimane 3-5)
├── 2a. Ottenere dataset Visium HD pubblico
├── 2b. Segmentazione cellulare su H&E (StarDist/Cellpose)
├── 2c. Analisi variazione intra-cellulare
├── 2d. Analisi spot di confine
├── 2e. Calibrazione library size e sparsità
└── 2f. Aggiornamento parametri nel simulatore
    Output: parametri calibrati su dati reali

FASE 3 — Validazione e raffinamento (settimane 5-7)
├── 3a. Confronto distribuzioni simulate vs reali
├── 3b. Test di deconvoluzione su dati simulati
├── 3c. GRF multipli per modulo genico
├── 3d. Modello nucleo/citoplasma (se dati lo supportano)
└── 3e. Benchmark con algoritmi di analisi
    Output: framework validato, pronto per paper
```
