# Step 1 — Cell Layer Design Doc

**Data inizio**: 2026-04-29
**Stato**: 🟡 In brainstorming (sezioni 1-3 approvate, 4-5 da discutere)
**Autori**: Luca Vedovelli (project owner), Claude (assist), feedback da Daniele
**Scope**: Solo Step 1 della nuova architettura proposta in `DS_review_aprile2026/REPORT_PROBLEMI_CRITICI.md`. Step 2 (mapping spot↔cellula, aggregazione, binning) sarà oggetto di un design doc separato.

---

## 1. Contesto e motivazione

Il framework attuale `geo_spatialtrans` genera profili trascrizionali **direttamente per ogni spot** della griglia, trattando ciascuno come entità indipendente. Nella realtà di Visium HD gli spot sono quadrati di 2×2 µm che coprono un tessuto fatto di **cellule** (diametro 15-30 µm): più spot coprono la stessa cellula, e gli spot ai confini cellulari catturano RNA da cellule diverse.

Manca quindi il livello intermedio tra tessuto e griglia: **le cellule**.

Step 1 introduce questo livello: un *cell layer* che, partendo dal clustering esistente, genera una mappa di cellule individuali con posizione, tipo, territorio e nucleo. Lo Step 2 (futuro) consumerà il cell layer per generare i profili trascrizionali a livello cellulare e aggregarli sugli spot Visium HD.

### Vincoli del design

- **Coesistenza** con la pipeline esistente: il vecchio `cell_df` (che è in realtà uno *spot grid*) resta funzionante. Step 1 introduce nuove strutture dati come percorso parallelo.
- **Riutilizzo massimo**: il clustering (`04_clustering.R`), la generazione del tessuto (`03b_generate_synthetic_tissue.R`) e il preprocessing immagine (`03_image_processing.R`) sono invariati.
- **Stile codebase**: coerenza con il pattern esistente (snake_case, file numerati, funzioni piccole con responsabilità singola; vedi modulo `06a-06k`).

---

## 2. Decisioni chiave

Tracciate cronologicamente; le decisioni RIVISTE sono il risultato del feedback ricevuto da Daniele dopo le prime schermate.

| # | Decisione | Stato |
|---|---|---|
| 1 | **Modello geometrico** = Voronoi puro con **nucleo ingrandito** (~40-60% di eq_radius). La cellula coincide col territorio Voronoi clipped alla maschera tessutale. Solo 2 zone: nucleo + citoplasma. Niente "body" inscritto separato. Niente ECM intra-regione esplicita (Voronoi non lascia gap). | ✅ |
| 2 | **Eterogeneità** = composizione esplicita per cluster. Catalogo `cell_types` separato dai cluster. Ogni cluster ha una `region_composition` (frazioni di tipi che lo popolano). Tipi infiltranti hanno proprietà PROPRIE (raggio nucleo, densità, forma). | ✅ |
| 3 | **Parametrizzazione** = 2 tabelle layered. (a) `cell_types`: catalogo dei tipi. (b) `region_composition`: mistura per cluster. Defaults dal `tissue_preset`; override fini via `override_cell_types` / `override_composition`. | ✅ |
| 4 | **Forma nucleo** = cerchio default per ogni tipo, ellisse opt-in per tipo. Forma del territorio emerge dal Voronoi+densità. | ✅ |
| 5 | **Architettura** = pipeline a fasi separate (5 funzioni: orchestrator + 4 step), allineata allo stile `06a-06k`. | ✅ |
| 6 | **Stack tecnico** = `deldir` (Voronoi) + `sf` (clipping/aree/spatial ops) + `polylabelr` (Pole of Inaccessibility per posizionare nuclei). Docker base = `rocker/geospatial`. | ✅ |
| 7 | **Coesistenza** con pipeline esistente: il vecchio `cell_df` resta. Riconciliazione delegata a Step 2. | ✅ |
| 8 | **Tissue presets** (livello sopra `cell_types`): 3 preset generici inclusi nel package — `"epithelial"`, `"tumor_microenv"`, `"stromal_rich"`. | ✅ |
| 9 | **Smossatura angoli** territori Voronoi = parametro `corner_smoothing` opzionale (default 0). Algoritmo Chaikin (corner cutting) o buffer-pos+buffer-neg. | ✅ |

**Refinement futuri esplicitamente rinviati a "Step 1.5":**
- Nucleo offset rispetto al centroide del territorio
- Nucleo lievemente irregolare (ellisse anisotropa o blob organico)
- Forma anisotropa del territorio (epitelio colonnare → bias direzionale del Voronoi)
- Eterogeneità ricca: gradienti di transizione + cluster localizzati di infiltrati (era opzione 3 della schermata eterogeneità)

---

## 3. Architettura e file structure

### Nuovi file (R/)

```
03_image_processing.R              [esistente]
03b_generate_synthetic_tissue.R    [esistente]
04_clustering.R                    [esistente]

▼ NUOVO: cell layer si interpone qui ▼
04b_cell_placement.R               [NUOVO]  orchestrator place_cells()
04b1_extract_regions.R             [NUOVO]  pixel cluster → blob connessi (sfc_POLYGON)
04b2_seed_centroids.R              [NUOVO]  Poisson-disk + composizione + densità
04b3_tessellate_voronoi.R          [NUOVO]  Voronoi clipped + smussatura
04b4_derive_cell_geometry.R        [NUOVO]  POI nucleo + raggi + aree
04b5_cell_layer_validation.R       [NUOVO]  sanity checks finali
▲

05_grid_sampling.R                 [invariato — Step 2 lo aggiornerà]
06a-06k_*.R                        [invariati — Step 2 li agganciaerà]
```

### Test e visualizzazione (R/testing/)

```
test_cell_layer.R         [NUOVO] esegue place_cells su tessuto sintetico e valida
visualize_cell_layer.R    [NUOVO] pannello: regioni + centroidi + Voronoi + nuclei
```

### Convenzione spaziale

| Stadio | Coordinate |
|---|---|
| Input `clust` | pixel (output del clustering) |
| Parametro globale | `pixel_size_um` (default `1.0`) |
| Output `cell_layer` | tutto in **µm** (centroidi, raggi, aree) |
| Step 2 | erediterà la stessa convenzione |

> ⚠️ Il codebase attuale ha incoerenze (`pixel_size_um=1` in `create_sampling_grid` ma `=10` in `full_test.R`). Step 1 non tocca questi file: il nuovo `place_cells()` riceve `pixel_size_um` esplicitamente dal chiamante.

---

## 4. Strutture dati e API

### 4.1 Output `cell_layer`

```r
cell_layer <- list(
  cell_df          = <data.frame>,    # una riga per cellula
  region_df        = <data.frame>,    # una riga per regione (blob connesso)
  region_polygons  = <sfc_POLYGON>,   # geometria delle regioni (clip mask)
  cell_territories = <sfc_POLYGON>,   # poligoni Voronoi clipped → la cellula
  cell_nuclei      = <sfc_POLYGON>,   # nuclei (cerchio o ellisse)
  metadata         = list(
    pixel_size_um, n_cells, n_regions,
    tissue_preset_used, corner_smoothing,
    cell_types_resolved, region_composition_resolved,
    random_seed, package_version, timestamp,
    validation                    # output di 04b5_*
  )
)
```

Il join è sempre per **`cell_id`** (1..N nello stesso ordine in `cell_df`, `cell_territories`, `cell_nuclei`).

### 4.2 Schema `cell_df`

| campo | tipo | unità | descrizione |
|---|---|---|---|
| `cell_id` | int | — | univoco |
| `region_id` | int | — | regione (blob connesso) di appartenenza |
| `cluster_id` | int | — | cluster di origine |
| `cell_type` | factor | — | tipo cellulare assegnato (catalogo `cell_types`) |
| `x`, `y` | num | µm | centroide del territorio Voronoi |
| `territory_area` | num | µm² | area del poligono Voronoi clipped |
| `eq_radius` | num | µm | √(territory_area/π) |
| `nucleus_shape` | char | — | `"circle"` (default) o `"ellipse"` |
| `nucleus_radius` | num | µm | raggio nucleo (se circle); NA se ellipse |
| `nucleus_major_r`, `nucleus_minor_r`, `nucleus_orient` | num | µm, rad | NA se circle |
| `nucleus_cx`, `nucleus_cy` | num | µm | centroide nucleo (= POI nel territorio) |
| `nucleus_area` | num | µm² | area del nucleo |
| `cytoplasm_area` | num | µm² | `territory_area − nucleus_area` |
| `is_minority_type` | logi | — | TRUE se `cell_type` ≠ tipo dominante della regione |

### 4.3 Schema `region_df`

| campo | tipo | unità | descrizione |
|---|---|---|---|
| `region_id` | int | — | univoco |
| `cluster_id` | int | — | cluster di origine |
| `dominant_cell_type` | factor | — | tipo con `fraction` massima nella `region_composition` |
| `area_um2` | num | µm² | area del blob |
| `perimeter_um` | num | µm | perimetro |
| `n_cells` | int | — | cellule piazzate nella regione |
| `target_density_weighted` | num | cell/mm² | media pesata `Σ fraction_i × density_i` |
| `achieved_density` | num | cell/mm² | `n_cells / area_um2 × 1e6` |

### 4.4 Catalogo `cell_types`

```r
cell_types <- data.frame(
  cell_type           = c("epithelial_dense", "stromal_loose", "immune_T",
                          "fibroblast", "fibroblast_CAF", "vascular", "extracellular"),
  density             = c(5500, 1500, 8000, 1200, 2000, 3000, 0),    # cell/mm²
  nucleus_shape       = c("circle", "circle", "circle", "circle",
                          "circle", "circle", "circle"),
  nucleus_radius      = c(NA, NA, NA, NA, NA, NA, NA),               # se NA → uso ratio
  nucleus_to_eq_ratio = c(0.45, 0.45, 0.55, 0.40, 0.42, 0.40, 0.0)
)
```

Significato di `nucleus_to_eq_ratio`: `nucleus_radius_actual = ratio × eq_radius`. Si adatta alla densità (cellula compressa → nucleo proporzionalmente più piccolo). Se l'utente specifica `nucleus_radius` assoluto in µm, quello prevale e `ratio` viene ignorato.

### 4.5 Composizione `region_composition`

```r
region_composition <- data.frame(
  cluster_id = c(1, 1, 1, 2, 2),
  cell_type  = c("epithelial_dense", "immune_T", "fibroblast", "stromal_loose", "immune_T"),
  fraction   = c(0.92, 0.05, 0.03, 0.95, 0.05)
)
```

Le frazioni per `cluster_id` devono sommare a 1.0 (validato; il package normalizza con warning se necessario). Default = `fraction = 1.0` per un singolo tipo dominante (omogeneità totale, ground truth pulito); l'utente abilita l'eterogeneità aggiungendo righe.

### 4.6 Tissue presets

| `tissue_preset` | Tipi inclusi | Composizione tipica per cluster |
|---|---|---|
| `"epithelial"` | epithelial_dense, stromal_loose, immune_T, immune_B | Cluster epiteliali: ~88% epi + ~5% stromale + ~5% immune. Cluster stromali: ~95% stromale + ~5% immune. |
| `"tumor_microenv"` | epithelial_tumor, fibroblast_CAF, immune_T, immune_B, stromal_loose | Cluster tumorali: ~70% tumor + ~15% TIL + ~10% CAF + ~5% stromale. Aggressivo su infiltrazione. |
| `"stromal_rich"` | stromal_dense, fibroblast, immune_T, vascular | Cluster stromali: ~75% stromale + ~10% fibroblast + ~10% vascular + ~5% immune. |

Tre preset generici di proposito — il package non vuole curare un atlante. La via "personalizzazione fine" passa da `override_*`.

### 4.7 API rivista

```r
place_cells(
  clust,                                # output di cluster_image()
  tissue_preset       = "epithelial",   # "epithelial" | "tumor_microenv" | "stromal_rich"
  cell_types          = NULL,           # NULL → ereditato dal preset
  region_composition  = NULL,           # NULL → ereditato dal preset
  override_cell_types = NULL,           # patch incrementale sul catalogo
  override_composition = NULL,          # patch incrementale sulla composizione
  cluster_to_dominant = NULL,           # mapping esplicito cluster_id → tipo dominante
                                        # se NULL: euristica intensità media + ordine cluster
  pixel_size_um       = 1.0,
  corner_smoothing    = 0,              # 0 = Voronoi puro, > 0 = Chaikin
  min_region_area_um2 = 100,            # regioni più piccole vengono ignorate
  random_seed         = 42,
  verbose             = TRUE
) -> cell_layer

# step interni (esposti per testabilità):
extract_regions(clust, pixel_size_um, min_region_area_um2)
seed_centroids(regions, cell_types, region_composition, random_seed)
tessellate_voronoi(centroids, region_polygons, corner_smoothing)
derive_cell_geometry(centroids, territories, cell_types)
```

---

## 5. Algoritmi chiave

### 5.1 `extract_regions()`
1. Da `clust` (pixel + `intensity_cluster`) ricostruisco una matrice etichetta `(img_height, img_width)`.
2. Per ogni `cluster_id` → matrice binaria → **componenti connesse** (4-vicinato) tramite `mmand::components()` (o algoritmo nativo BFS se evitiamo `mmand`).
3. Ogni componente connessa → boundary tracing (Moore) → poligono in **pixel** → conversione in µm via `× pixel_size_um` → `sfc_POLYGON` semplificato (`sf::st_simplify` con `dTolerance ≈ 0.5 µm`) per ridurre il numero di vertici.
4. Filtro per `area_um2 ≥ min_region_area_um2` (default 100 µm² ≈ 10×10 µm).
5. Output: `region_df` + `region_polygons` (sfc_POLYGON).

### 5.2 `seed_centroids()`
Per ogni `region_id`:
1. Determino la composizione (filtro `region_composition` per `cluster_id`).
2. Densità target ponderata: `Σ_i fraction_i × density_i`.
3. Numero target cellule: `n_cells_target = density_pondered × area_mm²`.
4. **Allocazione per tipo:** ripartisco `n_cells_target` per tipo secondo `fraction_i`, con arrotondamento *largest-remainder* per somma esatta.
5. **Piazzamento Poisson-disk** (Bridson, O(n)) con `r_min = ⅔ × eq_radius_target` per spaziatura. Esecuzione **per tipo in ordine decrescente di densità** (i tipi rari si infilano nei vuoti residui dei dominanti). Mask di accettazione = `region_polygon`.
6. Output: `centroids` data.frame con `cell_id`, `region_id`, `x`, `y`, `cell_type`, `cluster_id`.

> Razionale Poisson-disk vs jittered grid: distribuzioni naturali (no artifact), spaziatura garantita, tipi rari ottengono territori plausibilmente isolati.

### 5.3 `tessellate_voronoi()`
1. `deldir::deldir()` su tutti i centroidi (globali).
2. Per ogni cellula, estraggo il poligono Voronoi (`deldir::tile.list`).
3. **Clipping** con `sf::st_intersection(voronoi_poly, region_polygons[centroidi$region_id])`. Le cellule sui bordi del tessuto si tagliano correttamente.
4. **Smussatura** (se `corner_smoothing > 0`): Chaikin corner cutting con `n_iter = round(corner_smoothing × 3)` iterazioni (1, 2 o 3).
5. Output: `cell_territories` (sfc_POLYGON) + `territory_area`.

### 5.4 `derive_cell_geometry()`
Per ogni cellula:
1. `eq_radius = sqrt(territory_area / π)`.
2. **Centroide nucleo** = `polylabelr::poi(territory)` (Pole of Inaccessibility) — punto più "interno" del poligono, evita di mettere il nucleo in punte sottili.
3. **Raggio nucleo effettivo:**
   - Se `cell_types$nucleus_radius` è NA → `nucleus_radius_actual = nucleus_to_eq_ratio × eq_radius`
   - Altrimenti → `nucleus_radius_actual = min(nucleus_radius, 0.85 × distanza POI dal bordo)` (clamp di sicurezza: il nucleo non esce dalla cellula)
4. **Forma nucleo:**
   - `circle` → poligono circolare a 32 vertici intorno al POI
   - `ellipse` → ellisse con `nucleus_major_r`, `nucleus_minor_r`, `nucleus_orient`
5. Calcolo `nucleus_area`, `cytoplasm_area = territory_area − nucleus_area`.
6. Output: `cell_df` completo + `cell_nuclei` (sfc_POLYGON).

---

## 6. Validazione e testing

> 🟡 Da discutere e finalizzare nella Sezione 4 del brainstorming.

Bozza di sanity checks (`04b5_cell_layer_validation.R`):

- **Density check**: `|achieved − target| / target < 0.20` per regione (warning).
- **No overlap**: somma `territory_area` per regione ≈ area regione (tolleranza < 1%).
- **No nucleus outside territory**: `st_within(cell_nuclei, cell_territories)` deve essere TRUE per tutte (hard error).
- **Composition fidelity**: frazioni effettive vs target — warning se scarto > 10% per tipo.
- **Reproducibility**: stesso seed → stesso output.

Test scripts:
- `R/testing/test_cell_layer.R` — pipeline end-to-end su tessuto sintetico
- `R/testing/visualize_cell_layer.R` — pannello visivo

---

## 7. Hand-off a Step 2

> 🟡 Da discutere nella Sezione 4 del brainstorming.

Lista preliminare di "agganci" che Step 2 troverà pronti:

- `cell_layer$cell_df` come substrato per generazione profili trascrizionali a livello cellulare
- `cell_layer$cell_territories` + `cell_layer$cell_nuclei` per il mapping spot↔cellula (Sezione 3.2.A del REPORT)
- `cell_type` (catalogo) come chiave per i parametri di espressione cell-type-specific (06b)
- `nucleus_area` / `cytoplasm_area` per la variazione intra-cellulare (Sezione 3.2.D del REPORT)
- `metadata$pixel_size_um` come convenzione spaziale ereditata

---

## 8. Estensioni future (Step 1.5)

- Nucleo offset rispetto al centroide territorio
- Nucleo irregolare (ellisse o blob organico)
- Forma anisotropa del territorio (Voronoi pesato per modellare epitelio colonnare)
- Eterogeneità ricca: gradienti di transizione + cluster localizzati di infiltrati
- Ulteriori tissue presets (richiede curatela e validazione biologica)

---

## 9. Storico revisioni

| Data | Cosa | Da chi |
|---|---|---|
| 2026-04-29 | Decisioni 1-7 prese durante prima sessione brainstorming | Luca + Claude |
| 2026-04-29 | Decisioni 1-3 RIVISTE post-feedback Daniele (Voronoi + nucleo grande, composizione esplicita, catalogo disaccoppiato) | Daniele review |
| 2026-04-29 | Decisioni 8-9 aggiunte (tissue presets, corner smoothing) | Daniele review |
| 2026-04-29 | Sezioni 1-3 approvate e consolidate in questo doc | Luca |
| TBD | Sezione 4-5 (validazione + hand-off Step 2) | brainstorming in corso |
