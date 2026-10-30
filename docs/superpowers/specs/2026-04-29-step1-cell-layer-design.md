# Step 1 — Cell Layer Design Doc

**Data inizio**: 2026-04-29
**Versione**: 1.9 (2026-10-30, S1.6b: solo testo allineato al codice e alle decisioni già prese, nessun comportamento nuovo; la 1.8 è di S1.6; storico completo in §9; la 1.7 è di S1.5c). Versione 1.1 (2026-09-18, sessione S0) — incorpora la revisione di Daniele (commenti C-A, C-B), il check C8 e il pannello di 6 archetipi/preset. v1.0 = commit `7836829`.
**Stato**: 🟢 Design v1.9 — contratto dello Step 1, implementato e consolidato (sessioni S1.1–S1.6b) **[v1.9, S1.6b: la riga diceva «Design v1.1» fino a v1.8]**
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
| 1 | **Modello geometrico** = Voronoi puro con **nucleo ingrandito** (~40-60% di eq_radius). La cellula coincide col territorio Voronoi clipped alla maschera tessutale. Solo 2 zone: nucleo + citoplasma. Niente "body" inscritto separato. Niente ECM intra-regione esplicita (Voronoi non lascia gap). **[v1.9, S1.6b: la taglia del nucleo non è una frazione fissa di eq_radius: viene dal modello di taglia di §5.4 (v1.4), con i parametri ρ, β, σ di ogni preset in `R/04b0_tissue_presets.R`]** | ✅ |
| 2 | **Eterogeneità** = composizione esplicita per cluster. Catalogo `cell_types` separato dai cluster. Ogni cluster ha una `region_composition` (frazioni di tipi che lo popolano). Tipi infiltranti hanno proprietà PROPRIE (raggio nucleo, densità, forma). | ✅ |
| 3 | **Parametrizzazione** = 2 tabelle layered. (a) `cell_types`: catalogo dei tipi. (b) `region_composition`: mistura per cluster. Defaults dal `tissue_preset`; override fini via `override_cell_types` / `override_composition`. | ✅ |
| 4 | **Forma nucleo** = cerchio default per ogni tipo, ellisse opt-in per tipo. Forma del territorio emerge dal Voronoi+densità. | ✅ |
| 5 | **Architettura** = pipeline a fasi separate (1 orchestrator + 5 step interni: extract_regions → seed_centroids → tessellate_voronoi → derive_cell_geometry → validate_cell_layer), allineata allo stile `06a-06k`. | ✅ |
| 6 | **Stack tecnico** = ~~`deldir`~~ **GEOS `sf::st_voronoi` (v1.3, S1.3, 2026-10-08)** (Voronoi) + `sf` (clipping/aree/spatial ops) + `polylabelr` (Pole of Inaccessibility per posizionare nuclei). Docker base = `rocker/geospatial`. | ✅ |
| 7 | **Coesistenza** con pipeline esistente: il vecchio `cell_df` resta. Riconciliazione delegata a Step 2. | ✅ |
| 8 | **Tissue presets** (livello sopra `cell_types`): **6 preset** inclusi nel package, uno per archetipo del pannello A1–A6 (v1.1; erano 3 in v1.0): `"epithelial"` (A1), `"tumor_microenv"` (A2), `"stromal_rich"` (A3), `"lymphoid"` (A4), `"brain_cortex"` (A5), `"columnar_muscle"` (A6, **controllo negativo previsto** per il Voronoi isotropo). I parametri numerici dei preset sono **provvisori** finché non compaiono in `docs/BIO_REFERENCES.md` con fonte (binario R, sessioni R1–R3). | ✅ v1.1 |
| 9 | **Smossatura angoli** territori Voronoi = parametro `corner_smoothing` opzionale (default 0). Algoritmo Chaikin (corner cutting) o buffer-pos+buffer-neg. | ✅ |
| 10 | **Scala px→µm obbligatoria (C-A, Daniele)**: `pixel_size_um` **non ha default** in `place_cells()`. O lo fornisce il chiamante (immagini a risoluzione reale: raccomandato 1 px = 1 µm) oppure viene stimato esplicitamente con `estimate_pixel_size(clust, cell_types, ...)`, che ricava la scala dalla dimensione attesa delle cellule del tipo dominante e la scrive in `metadata$pixel_size_um` con `metadata$pixel_size_source = "estimated"` e un WARN. Nuovo check di plausibilità **C8** sull'`eq_radius` mediano. **v1.6 (S1.5)**: la stima dall'immagine non ha superato la controprova (CP-1: scala trovata a 1 µm/px solo nel linfonodo; numeri sbagliati senza rifiuto) ed esce dall'uso raccomandato: `place_cells()` la usa solo con `allow_estimate = TRUE`. Fonti raccomandate: il valore dell'utente; `scalefactors_json.json` dello Space Ranger controllato contro la griglia dei bin (esatto sui 5 dataset); mai il tag di risoluzione del TIFF (BL-016). | ✅ v1.1, rivista v1.6 |
| 11 | **Diffusione laterale degli mRNA (C-B, Daniele)** è un requisito dello **Step 2** (kernel gaussiano isotropo con parametro `diffusion_sigma_um`, da calibrare su dati reali — binario R5). Per Step 1 non cambia la geometria: `cell_territories` e `cell_nuclei` sono le **sorgenti** del kernel (sez. 7.1). | ✅ v1.1 |

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
06a-06k_*.R                        [invariati — Step 2 li aggancerà]
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
| Parametro globale | `pixel_size_um` — **obbligatorio, senza default** (decisione 10 / C-A). Immagine a risoluzione reale: raccomandato 1 px = 1 µm; altrimenti l'utente lo specifica o lo legge da `scalefactors_json` (v1.6; la stima dall'immagine solo con `allow_estimate = TRUE`) |
| Output `cell_layer` | tutto in **µm** (centroidi, raggi, aree) |
| Step 2 | erediterà la stessa convenzione |

> ⚠️ Il codebase attuale ha incoerenze (`pixel_size_um=1` in `create_sampling_grid` ma `=10` in `full_test.R`). Step 1 non tocca questi file: il nuovo `place_cells()` riceve `pixel_size_um` esplicitamente dal chiamante.

> **C-A (Daniele, revisione v1.0)**: una volta scelti tessuto e tipi cellulari, la dimensione dello spot Visium HD (2×2 o 8×8 µm) va calcolata **in proporzione al tessuto**, dividendo l'immagine in una griglia proporzionale. Se l'immagine ha risoluzione reale, 1 px = 1 µm; altrimenti l'utente **deve** specificare la scala. Conseguenze: (1) `pixel_size_um` obbligatorio (4.7); (2) `metadata$pixel_size_um` esposto da Step 1 ed ereditato da Step 2, che costruisce la griglia 2 µm da lì; (3) check C8 (6.1).

**Helper `estimate_pixel_size(clust, cell_types, cluster_to_dominant = NULL, method = "nuclear_spacing")`** — stima la scala quando è ignota. Idea: la spaziatura media fra nuclei (o l'eq_radius atteso) del tipo dominante è nota in µm dal catalogo (`sqrt(1e6 / (π × density))`); misurando la stessa grandezza in **pixel** sull'immagine (es. distanza al primo vicino dei massimi locali di intensità nel cluster dominante) il rapporto dà µm/px. Restituisce `list(pixel_size_um, method, n_features, ci)`. Non è mai silenzioso: `place_cells()` lo invoca solo se l'utente passa `pixel_size_um = "estimate"`, e il risultato finisce in `metadata` con un WARN. La sua validazione (recupero di scale note 1, 2, 5, 10 µm/px) è la controprova della sessione S1.5. **v1.6 (S1.5)**: implementata con il laplaciano della gaussiana sulla scurezza del cluster dominante, rifiuto se i nuclei non sono risolti o sono meno di 30, controllo di coerenza della taglia (`size_check`); metodo `"metadata"` aggiunto. La controprova l'ha smentita (CP-1, CP-1b, BL-145): non è raccomandata.

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
    pixel_size_um, pixel_size_source,   # "user" | "estimated" (C-A)
    n_cells, n_regions,
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
| `x`, `y` | num | µm | posizione del generatore (centro seminato della cellula; v1.6, S1.5 D-9: prima scritto «centroide del territorio Voronoi») |
| `territory_area` | num | µm² | area del poligono Voronoi clipped |
| `eq_radius` | num | µm | √(territory_area/π) |
| `nucleus_shape` | char | — | `"circle"` (default) o `"ellipse"` |
| `nucleus_radius` | num | µm | raggio nucleo (se circle); NA se ellipse |
| `nucleus_major_r`, `nucleus_minor_r`, `nucleus_orient` | num | µm, rad | NA se circle |
| `nucleus_cx`, `nucleus_cy` | num | µm | centro del nucleo (generatore con «clip», POI del territorio con «scale», v1.5) |
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
| `target_density_weighted` | num | cell/mm² | media **armonica** `1 / Σ (fraction_i / density_i)` (v1.8: testo allineato a S1.2, D3, già applicata in §5.2) |
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

**(v1.4)** Colonne opzionali: `nucleus_area_slope` (β, default 1), `nucleus_area_sdlog` (σ, default 0), `nucleus_aspect` (default 1), `nucleus_orient` (default `"random"`); vedi §5.4. **(v1.5)** `nucleus_nc_max` (φ, tetto del rapporto nucleo/territorio, in (0, 1], default 1; solo modello «clip»). Significato di `nucleus_to_eq_ratio` (con β = 1, σ = 0): `nucleus_radius_actual = ratio × eq_radius`. Si adatta alla densità (cellula compressa → nucleo proporzionalmente più piccolo). Se l'utente specifica `nucleus_radius` assoluto in µm, quello prevale e `ratio` viene ignorato. **(v1.8)** Un tipo con `density = 0` non può avere frazione > 0 nella composizione (`seed_centroids()` dà errore): lo spazio extracellulare non è un tipo cellulare ma l'area fuori dalle regioni o dentro i vuoti della semina (§5.2), quindi la riga `extracellular` dell'esempio sopra è solo illustrativa.

### 4.5 Composizione `region_composition`

```r
region_composition <- data.frame(
  cluster_id = c(1, 1, 1, 2, 2),
  cell_type  = c("epithelial_dense", "immune_T", "fibroblast", "stromal_loose", "immune_T"),
  fraction   = c(0.92, 0.05, 0.03, 0.95, 0.05)
)
```

Le frazioni per `cluster_id` devono sommare a 1.0 (validato; il package normalizza con warning se necessario). Default = `fraction = 1.0` per un singolo tipo dominante (omogeneità totale, ground truth pulito); l'utente abilita l'eterogeneità aggiungendo righe.

### 4.6 Tissue presets — pannello A1–A6 (v1.1)

Il pannello fissato il 2026-09-18 (roadmap §3) entra nel design come **6 preset**, uno per archetipo. Ogni test dello Step 1 gira su tutti e sei (principio di generalità). I tipi elencati e le composizioni "tipiche" sono **ipotesi di partenza**: un preset è dichiarato biologicamente solido solo quando i check B passano contro il dataset reale della sua riga (binario R), e i suoi numeri (densità, `nucleus_to_eq_ratio`) entrano in `docs/BIO_REFERENCES.md` con fonte.

| `tissue_preset` | Archetipo | Tipi inclusi (ipotesi) | Composizione tipica per cluster (ipotesi) | Attesa qualitativa da verificare |
|---|---|---|---|---|
| `"epithelial"` | A1 epitelio semplice denso | epithelial_dense, stromal_loose, immune_T, immune_B | epiteliale: ~88% epi + ~5% stromale + ~5% immune; stromale: ~95% stromale + ~5% immune | densità alta, poligoni compatti e regolari |
| `"tumor_microenv"` | A2 tumore solido con microambiente | epithelial_tumor, fibroblast_CAF, immune_T, immune_B, stromal_loose | tumorale: ~70% tumor + ~15% TIL + ~10% CAF + ~5% stromale | composizione mista, alta variabilità di area |
| `"stromal_rich"` | A3 stroma lasso / connettivo | stromal_dense, fibroblast, immune_T, vascular | stromale: ~75% stromale + ~10% fibroblast + ~10% vascular + ~5% immune | densità bassa, poligoni grandi, N/C basso |
| `"lymphoid"` | A4 tessuto linfoide | immune_B, immune_T, fdc (follicular dendritic), macrophage | follicolo: ~80% B + ~15% T + ~5% altro; zona T: ~70% T + ~25% B + ~5% altro | densità molto alta, N/C ≈ 0.8–0.9, poligoni minuscoli |
| `"brain_cortex"` | A5 cervello (corteccia) | neuron, astrocyte, oligodendrocyte, microglia, vascular | grigia: ~40% neuroni + ~55% glia + ~5% vascolare; bianca: ~90% glia | densità bassa-media, distribuzione **bimodale** delle aree |
| `"columnar_muscle"` | A6 muscolo / epitelio colonnare | myocyte (o epithelial_columnar), fibroblast, vascular | ~85% miociti + ~10% fibroblast + ~5% vascolare | **controllo negativo previsto**: il Voronoi isotropo non riproduce cellule allungate → i check B su forma (eccentricità) devono **fallire**; quantificare quanto e come → motiva Step 1.5 |

Sei preset di proposito, non un atlante: la personalizzazione fine passa da `override_*`. I tre preset nuovi (A4–A6) sono ammessi nel package solo con i parametri marcati "provvisori" nei metadata (`cell_types_resolved$provisional = TRUE`) finché BIO_REFERENCES non li copre. **[v1.9, S1.6b: superato; la condizione è soddisfatta e la colonna non esiste: vedi «Stato attuale dei preset» qui sotto]**

**(v1.8, S1.6)** Stato attuale dei preset. (a) I tipi, le frazioni e i modelli di taglia per tipo vengono dai dati (R6, S1.5b: DEC-013...DEC-019), non dalle ipotesi della tabella sopra; il preset è **una composizione per tessuto applicata a ogni cluster** (`R/04b_cell_placement.R`), non una composizione per tipo di cluster. (b) La condizione della colonna `provisional` è soddisfatta: ogni valore dei sei preset ha una riga con fonte in `docs/BIO_REFERENCES.md` (colonne `src_*` di `R/04b0_tissue_presets.R`, verificate dal test di S1.5); la colonna non esiste e non serve. I valori ancora provvisori sono quelli dichiarati tali nelle decisioni e nelle loro righe (per esempio lo stroma misto, DEC-018, e la composizione per zone del linfonodo, DEC-015). (c) Riferimenti per tipo della taglia di linfonodo e cuore sull'insieme del modello di taglia (B-265...B-269, DEC-030).

### 4.7 API rivista

*(v1.6, S1.5: firma implementata in `R/04b_cell_placement.R`; i semi degli stadi, le regole di `stride`, `small_components` e `size_conversion` sono definiti in `docs/steps/S1.5.md`)*

```r
place_cells(
  clust,                                # output di cluster_image() (colonne x, y, cluster_col)
  tissue_preset       = "epithelial",   # "epithelial" | "tumor_microenv" | "stromal_rich" |
                                        # "lymphoid" | "brain_cortex" | "columnar_muscle"  (A1–A6)
  cell_types          = NULL,           # NULL → ereditato dal preset
  region_composition  = NULL,           # NULL → ereditato dal preset
  override_cell_types = NULL,           # patch incrementale sul catalogo
  override_composition = NULL,          # patch incrementale sulla composizione
  cluster_to_dominant = NULL,           # mapping esplicito cluster_id → tipo dominante
  pixel_size_um,                        # OBBLIGATORIO (C-A): numero (µm/px), list(scalefactors_json = <file>),
                                        # oppure "estimate" (solo con allow_estimate = TRUE, v1.6)
  corner_smoothing    = 0,              # 0 = Voronoi puro, > 0 = Chaikin
  min_region_area_um2 = 100,
  random_seed         = 42,             # da qui tre semi indipendenti: semina, taglia, nuclei pallidi
  cluster_col         = "intensity_cluster",
  stride              = NULL,           # NULL → attributo "stride" di clust; errore se rilevato > 1 e non dichiarato (BL-042)
  simplify_tol_um     = 0.5,
  small_components    = c("merge", "exclude"),        # BL-043
  nucleus_model       = "clip",                       # v1.5
  size_conversion     = c("preserve_area", "none"),   # conversione del modello di taglia alla densita' del preset (S1.5, scelta 2)
  validate = TRUE, strict = TRUE,       # strict: un ERROR dei controlli ferma la funzione
  keep_tiles          = FALSE,          # restituisce i tile di Voronoi (regola delle celle interne)
  allow_estimate      = FALSE,          # v1.6: la stima dall'immagine non e' raccomandata (S1.5, CP-1)
  verbose             = TRUE
) -> cell_layer

# v1.8 (S1.6): argomenti aggiunti da S1.5b-S1.6 (firma completa in R/04b_cell_placement.R):
#   types = getOption("geo.spatialtrans.preset_types", "multi"), size_arm = NULL, composition = "tissue",
#   spatial = NULL ("uniform" | "compartments"), min_dist_mode = c("tissue", "nucleus"), labels = c("field", "blind"),
#   override_compartments = NULL, override_heterogeneity = NULL, type_order = NULL,
#   quotas = getOption("geo.spatialtrans.quotas", FALSE), thinning_bound = getOption("geo.spatialtrans.thinning_bound", "exact"),
#   size_resid = NULL, resid_mode = c("decile", "pooled")
# Opzioni di sessione (S1.6, DEC-030): geo.spatialtrans.preset_types ("single" = preset a un tipo di S1.5),
#   geo.spatialtrans.quotas e geo.spatialtrans.thinning_bound (TRUE e "grid" riproducono S1.5c-S1.5e bit per bit).

# step interni (esposti per testabilità):
estimate_pixel_size(clust, cell_types, cluster_to_dominant = NULL, method = c("nuclear_spacing", "metadata"),
                    scalefactors_json = NULL, value_col = "value", cluster_col = "intensity_cluster",
                    n_boot = 200L, block_px = 64L, size_check = TRUE)
validate_cell_layer(cell_layer, c8_interval = NULL, pixel_size_source = "user", grid_um = 1e-6,
                    sum_rel_tol = 1e-9, min_region_area_um2 = 0, type_order = NULL)   # C1–C9; type_order per C7 (v1.8)
overlay_prec(op, x, y, grid_um = 1e-6, origin = NULL, area = TRUE)                     # unico helper degli overlay (F1, BL-081)
plot_cell_layer(cell_layer, mode = c("voronoi", "types", "density", "regions"))
extract_regions(clust, pixel_size_um, min_region_area_um2)
seed_centroids(regions, cell_types, region_composition, random_seed, dmin_factor = 2/3, max_attempts_factor = 1000,
               max_consecutive_rejections = 10000, type_order = c("sequential", "blind"), spatial = c("uniform", "compartments"),
               compartments = NULL, labels = c("field", "blind"), n_features = 256, heterogeneity = NULL,
               quotas = FALSE, thinning_bound = c("exact", "grid"))                        # v1.8
tessellate_voronoi(centroids, region_polygons, corner_smoothing)
derive_cell_geometry(centroids, territories, cell_types, ...)   # firma completa in §5.4 (piu' verbose, resid_table, resid_mode da S1.5d)
```

---

## 5. Algoritmi chiave

### 5.1 `extract_regions()`
**[v1.9, S1.6b: i passi 2-4 sono il testo di v1.0. Implementato (S1.1, approvato da Luca il 2026-10-07, ROADMAP §6): componenti 4-connesse con union-find sui tratti di riga (`.er_label_runs()`, `.er_union_find()`) invece della BFS; contorno sui lati dei pixel (`.er_component_polygons()`, area esatta, buchi inclusi) invece del tracciamento di Moore sui centri; filtro di area sull'area esatta in pixel. `extract_regions()` esclude ancora le componenti sotto `min_region_area_um2` (passo 4), ma `place_cells()` per predefinito le unisce prima al cluster adiacente (`merge_small_components()`, `small_components = "merge"`, BL-043, v1.6); con `small_components = "exclude"` vale il passo 4.]**
1. Da `clust` (pixel + `intensity_cluster`) ricostruisco una matrice etichetta `(img_height, img_width)`.
2. Per ogni `cluster_id` → matrice binaria → **componenti connesse** (4-vicinato) tramite BFS nativo R (vedi nota implementativa sotto).
3. Ogni componente connessa → boundary tracing (Moore) → poligono in **pixel** → conversione in µm via `× pixel_size_um` → `sfc_POLYGON` semplificato (`sf::st_simplify` con `dTolerance ≈ 0.5 µm`) per ridurre il numero di vertici.
4. Filtro per `area_um2 ≥ min_region_area_um2` (default 100 µm² ≈ 10×10 µm). I pixel che cadevano in regioni escluse vengono **trattati come spazio extracellulare macro** in Step 2 (== zona below threshold del clustering originale): nessun centroide viene piazzato lì e nessun territorio Voronoi le copre.
5. Output: `region_df` + `region_polygons` (sfc_POLYGON).

> **Implementazione componenti connesse**: usiamo BFS nativo R (queue + visited matrix), non `mmand::components()`. Motivo: evitare una nuova dipendenza (mmand non è già nel codebase) e mantenere il pacchetto più leggero. La performance è sufficiente per immagini fino a ~5000×5000 px.

### 5.2 `seed_centroids()`
Per ogni `region_id`:
1. Determino la composizione (filtro `region_composition` per `cluster_id`).
2. Densità target ponderata: `Σ_i fraction_i × density_i`. **[v1.2, S1.2, D3: sostituita da `1 / Σ_i (fraction_i / density_i)` — media armonica; vedi §9]**
3. Numero target cellule: `n_cells_target = density_pondered × area_mm²`.
4. **Allocazione per tipo:** ripartisco `n_cells_target` per tipo secondo `fraction_i`, con arrotondamento *largest-remainder* per somma esatta.
5. **[v1.2, S1.2: RSA a n fissato invece di Bridson, che è massimale; vedi §9]** **Piazzamento Poisson-disk** (Bridson, O(n)) con `r_min = ⅔ × eq_radius_target` per spaziatura **[v1.9, S1.6b: la distanza minima è `min_dist_um` del tipo, nei preset la distanza d* stimata dai tessuti (DEC-009, S1.2); il fattore ⅔ (`dmin_factor`) vale solo per i tipi del catalogo senza `min_dist_um`; ordine dei tipi e semina cieca: «Stato attuale della semina» più sotto]**. Esecuzione **per tipo in ordine decrescente di densità** (i tipi rari si infilano nei vuoti residui dei dominanti). Mask di accettazione = `region_polygon`.
6. Output: `centroids` data.frame con `cell_id`, `region_id`, `x`, `y`, `cell_type`, `cluster_id`.

**(v1.8, S1.6) Stato attuale della semina.** (a) **Preset per tipo**: semina **cieca al tipo** (un solo compito di semina per regione, tipi assegnati dopo con i conteggi ripartiti; DEC-019), così il territorio non dipende dal tipo; i cataloghi dell'utente restano in sequenza salvo `type_order = "blind"`. (b) **Semina a compartimenti** (opzione, `spatial = "compartments"`, v1.7): campi latenti per tipo, vuoti che spostano le cellule (famiglia RV, DEC-022, DEC-024); da S1.6 **senza quote per tipo** per predefinito (DEC-025, regola Q di S1.5e) e con il **limite esatto del diradamento**: le proposte uniformi sono accettate con probabilità w(x) / w_max, w = 1 / Σ_t p_t(x)/ρ_t ≤ max_t ρ_t dei tipi della regione, per il massimo del fattore della densità (1 senza campo della densità; il tetto H con `dens_cap`); con un campo della densità senza tetto il limite è 1.1 volte il massimo sulla griglia, registrato (DEC-030). (c) Il predefinito dei preset resta la semina uniforme (regola B di S1.5c-S1.5e, DEC-020, DEC-022, DEC-024).

> Razionale Poisson-disk vs jittered grid: distribuzioni naturali (no artifact), spaziatura garantita, tipi rari ottengono territori plausibilmente isolati.

### 5.3 `tessellate_voronoi()`
*(v1.3, S1.3, 2026-10-08/09 — testo riscritto; le versioni precedenti sono nella storia git)*
Firma: `tessellate_voronoi(centroids, regions, corner_smoothing = 0, keep_tiles = FALSE, verbose = TRUE)`, con `regions` = output di `extract_regions()`.
1. **Per ogni regione** (D-S1.3.1): diagramma di Voronoi dei soli centroidi della regione con GEOS (`sf::st_voronoi`, `point_order = TRUE`; decisione 6 rivista), nel riquadro della regione allargato di max(1 µm, 10 % del lato); 1 centroide → territorio = regione.
2. **Ritaglio**: solo i tile non contenuti nella regione vengono intersecati con la regione (ritagliata localmente al riquadro del tile).
3. **Frammenti** (D-S1.3.2): se il ritaglio spezza una cella, il pezzo con il generatore resta; ogni pezzo orfano va al territorio della stessa regione con il confine condiviso più lungo (tratto > 1e-6 µm, misurato entro 1e-7 µm; decisioni simultanee per passata, la prima contro i soli pezzi principali; a pari merito il `cell_id` minore; il donatore non è candidato). Un orfano senza confine (solo regioni MULTIPOLYGON) resta al generatore.
4. **Smussatura** (D-S1.3.3, se `corner_smoothing > 0`): Chaikin (taglio 1/4–3/4) con `n_iter = max(1, round(3 × cs))`, cs ∈ (0, 1], su ogni parte del territorio, intersecato con la parte originale; se una parte si spezza resta il pezzo con il generatore. Le lacune sono spazio extracellulare (`region_check$gap_area`).
5. Output: `cell_territories` (ordine di `cell_id`), `territory_df`, `region_check`, `tiles` (opzionale), `info` (frammenti, agganci, isolati, riparazioni). `st_make_valid()` solo come ultima risorsa, contato.
> Avvertenza (BL-069): i vertici condivisi possono differire all'ultima cifra; le operazioni di sovrapposizione GEOS su questi territori (unione, intersezione) possono dare risultati sbagliati — usare predicati o una precisione fissa.

### 5.4 `derive_cell_geometry()`
*(v1.5, S1.4d, 2026-10-20 — modello «clip» predefinito, decisione di Luca del 2026-10-16 su BL-098; il modello «scale» di v1.4 resta, invariato)*
Firma: `derive_cell_geometry(centroids, territories, cell_types, random_seed, nucleus_model = c("clip", "scale"), nucleus_center = NULL, n_vertices = 32, clamp_factor = 0.85, clip_inset_um = 1e-5, snap_grid_um = 1e-6, area_rel_tol = 1e-7, poi_precision_um = 1e-3, strict = TRUE)`.
1. `eq_radius = sqrt(territory_area / π)`; flusso casuale `set.seed(random_seed); z <- rnorm(n); u <- runif(n)` in ordine di `cell_id` (uguale nei due modelli: stesse aree bersaglio con lo stesso seme).
2. **Centro**: `nucleus_center = NULL` → generatore con «clip», POI di `polylabelr` con «scale»; `d` = distanza GEOS dal bordo del territorio (errore se d = 0).
3. **Taglia (D-S1.4.1)**: `log A_N = log(ρ²·A_ref) + β·(log A_T − log A_ref) + σ·z`, `A_ref = 1e6/density` del tipo; `nucleus_radius` assoluto → `log A_N = log(π r²) + σ·z`. Default β = 1, σ = 0. Stima: `fit_nucleus_size_model()`. S1.4: β = 0.11–0.62 sui 6 archetipi (B-091…B-108).
4. **Forma**: poligono di `n_vertices` vertici su ellisse con aspetto `nucleus_aspect` (1 = cerchio), orientazione `nucleus_orient`; semiassi moltiplicati per k = sqrt(2π/(n sin(2π/n))) così che l'area del poligono sia esattamente l'area voluta.
5. **Modello «clip» (predefinito; K-2 di S1.4c)**. Tetto `A* = min(A_N, φ·min(A_T, area(T_in)))`, φ = colonna `nucleus_nc_max` del tipo (in (0, 1], default 1; valori del pannello: q99 del rapporto nucleo/territorio reale: A1 0.757 (B-121), A2 1 (B-122, DEC-005), A3 0.947 (B-123), A4 0.791 (B-143), A5 0.474 (B-125), A6 0.149 (B-154); v1.6: A4 e A6 dalla maschera unica di R3b (DEC-006), prima 0.816 e 0.180 (B-124, B-126)). `T_in` = territorio rientrato di `clip_inset_um` (buffer senza modello di precisione, poi sulla griglia `snap_grid_um`; errore se vuoto: BL-118). Nucleo = poligono(s) ∩ `T_in` (intersezione sulla griglia), con la scala s ≥ 1 tale che area = A*: s = 1 per tutte le celle, poi solo per le celle ritagliate primo passo s = √(A*/area(1)) e falsa posizione di Illinois con ricaduta sulla bisezione; arresto a |area − A*|/A* ≤ `area_rel_tol` o quando lo spostamento dei vertici fra gli estremi è < griglia/4. Nuclei contenuti nel proprio territorio (C3) e disgiunti per costruzione (fessura di 2·`clip_inset_um` fra nuclei a contatto). Colonne: `nucleus_area_target` (A_N), `nucleus_area_cap` (A*), `nucleus_capped`, `nucleus_clipped`, `nucleus_scale` (s), `nucleus_n_iter`, `center_dist`; `nucleus_radius` = raggio del cerchio di pari area; `cell_nuclei` convertiti in `sf` via WKB, `sfc_MULTIPOLYGON` se un nucleo ritagliato è in più parti: allora tutti i nuclei dell'insieme sono di tipo MULTIPOLYGON e il numero di parti si legge da `lengths()` (v1.8, S1.6: nel banco K2 di S1.4e 532 259 nuclei di tipo MULTIPOLYGON, 180 davvero in più parti; BL-173).
6. **Modello «scale» (D-S1.4.2, v1.4, invariato)**: fattore uniforme s ≤ 1; cerchio: raggio circoscritto ≤ `clamp_factor`·d; ellisse: dal limite del cerchio circoscritto, bisezione finché il poligono dilatato di 1/`clamp_factor` resta nel territorio. Colonne `nucleus_area_target`, `nucleus_scale`, `nucleus_clamped`, `center_dist`. Riproduce i risultati fino a S1.4c (S1.4d, cancello G-4b).
7. **Operatore di osservazione per i check B sui nuclei simulati (BL-112, fissato in S1.4c, scritto qui in S1.4d)**: nuclei rasterizzati sulla griglia nativa del ROI (regola del centro del pixel, `tools/S1.4c_raster.py`), centro = centroide dei pixel, territorio osservato = Voronoi dei centroidi osservati con `tessellate_voronoi()` sulla regione del ROI, lo stesso operatore dei nuclei reali; contatto = pixel di bordo entro distanza di Chebyshev ≤ 2 px da un'altra etichetta.
8. **Frazione di nuclei a contatto sulle finestre (S-1; BL-131, S1.4d)**: riferimento = valore manuale di Luca per archetipo (B-127…B-132, decisione del 2026-10-20); simulati = nuclei con centroide raster nella finestra a distanza ≥ d_max dai lati e dalle esclusioni (regola dei focali di R2c); per seme la media delle frazioni di finestra **a peso uguale** (come il riferimento, 15 nuclei per finestra), poi la media sui semi.
9. Output: `cell_df` (§4.2; `is_minority_type` è dell'orchestratore), `cell_nuclei`, `info`.

---

## 6. Validazione e testing

### 6.1 Sanity checks (`04b5_cell_layer_validation.R`)

Eseguiti dopo la pipeline e riportati in `cell_layer$metadata$validation`. Severità con azione:

| # | Check | Severità | Azione |
|---|---|---|---|
| C1 | `\|achieved_density − target_density\| / target < 0.20` per regione | **WARN** | log + scrivi in metadata, non blocca |
| C2 | `Σ territory_area` per regione ≈ area regione (v1.8: tolleranza relativa `sum_rel_tol` = 1e-9 più la perdita dichiarata della smussatura; coppie di territori disgiunte e territori dentro la regione entro g·P + g² sulla griglia g = 1e-6 µm; S1.3b, S1.3d, S1.5) | **ERROR** | hard fail: bug Voronoi/clipping |
| C3 | `st_within(cell_nuclei, cell_territories)` per ogni cellula | **ERROR** | hard fail: bug derive_geometry |
| C4 | Frazione effettiva tipo vs target — v1.8: scarto assoluto `\|osservata − target\|` > 0.10, solo nelle regioni con almeno 10 cellule | **WARN** | log + suggerimento (regione troppo piccola per tipi rari) |
| C5 | `nucleus_area > 0` per ogni cellula | **ERROR** | hard fail |
| C6 | `n_cells > 0` per ogni `region_id` con `area_um2 ≥ min_region_area_um2` | **WARN** | log: parametri sbagliati |
| C7 | `eq_radius` per `cell_type` coerente con il valore atteso | **WARN** | log: 90° percentile dentro `[0.5×, 2.0×]` di `expected_eq_radius` |
| C8 | **Plausibilità della scala (C-A, v1.1)**: `eq_radius` mediano delle cellule piazzate dentro un intervallo plausibile per cellule di mammifero; valore **provvisorio** [5, 25] µm dalla revisione di Daniele, da sostituire con l'intervallo per archetipo di `BIO_REFERENCES.md` (R2–R3) | **WARN** (ERROR se `pixel_size_source = "estimated"`) | log: la scala px→µm è probabilmente sbagliata; riportare `pixel_size_um` usato e mediana osservata. **v1.6 (S1.5)**: intervalli per archetipo nei preset (B-109…B-111, B-142, B-113, B-153). Con la densità del catalogo in cellule/mm² una scala sbagliata cambia il numero di cellule, non il loro `eq_radius` (test E10d): C8 vede una densità in unità sbagliate (test E10c), non una scala sbagliata. Deciso da Luca il 2026-10-21 (DEC-011) |
| C9 | (v1.8, da S1.5) validità delle geometrie, ordine di `cell_id` e legame righe-geometrie fra `cell_df`, `cell_territories`, `cell_nuclei`, `region_df$n_cells` (vincoli di §7.2) | **ERROR** | hard fail |
| TEST | `set.seed → output identico` | **TEST** | unit test, non runtime |

`expected_eq_radius` (per C7) si deriva dalla `density` **del cell_type** (non dalla densità ponderata della regione), perché il check valuta il singolo tipo a prescindere dal suo contesto regionale:
```
expected_eq_radius_um[cell_type] = sqrt(1e6 / (π × density[cell_type]))
# es. density=5500 → expected ≈ 7.6 µm
#     density=1500 → expected ≈ 14.6 µm
```
**Nota**: in regioni miste i territori Voronoi reali dipendono dalla composizione locale (un linfocita raro in epitelio denso avrà territorio più grande della sua "naturale" densità isolata). Il check C7 è quindi una **WARN soft** — sforare significa che il tipo si trova in un contesto sbagliato di densità, non necessariamente un bug.

**(v1.8, S1.6, BL-156; decisione di Luca, [DEC-031](/decisions/DEC-031.md))** Con la semina **cieca al tipo** (preset per tipo con semina uniforme) il territorio non dipende dal tipo e la densità del tipo puro dei preset (ρ_t = D·m̄/m_t) è solo un'area di riferimento del modello di taglia: C7 usa allora, cella per cella, la densità della miscela della regione (`target_density_weighted`), e il controllo è il 90° percentile di `eq_radius / expected_eq_radius_miscela` in [0.5, 2]. Con la semina in sequenza e con i compartimenti resta la densità del tipo, come sopra. Nel banco di S1.6: 0 WARN su 300 unità contro 2 (cervello) con la regola di prima.

### 6.2 Test scripts

**(v1.8, S1.6)** I test sono per step: `R/testing/test_<step>.R` con il loro runner `run_<step>.sh` (`test`, `mutants`; mutanti come copie modificate del codice, `tools/mutate.py`), e la batteria di S1.6 (`bash R/testing/run_S1.6.sh test`) esegue tutte le precedenti. I casi limite della tabella sotto sono coperti da `test_S1.1.R`...`test_S1.6.R` (in particolare `test_S1.5.R`); i tre script elencati qui non esistono.

```
R/testing/test_cell_layer.R              [PRIMARY] integrazione end-to-end
R/testing/test_cell_layer_edge_cases.R   [robustezza]
R/testing/visualize_cell_layer.R         [visivo]
```

**`test_cell_layer.R`** — esegue la pipeline su tessuto sintetico e valida:

```r
syn <- generate_synthetic_tissue(600, 600, complexity = 2, seed = 42)
clust <- cluster_image(syn$img_df_thresh, k = 4,
                       random_seed = 42, spatial_weight = 0.6)

# v1.1: il test gira su TUTTI i 6 preset (A1–A6); qui il caso A1
cell_layer <- place_cells(
  clust,
  tissue_preset = "epithelial",
  pixel_size_um = 1.0,          # esplicito (C-A): nessun default
  random_seed   = 42
)

# Asserzioni:
#   - n_cells > 1000
#   - n_regions >= 4
#   - metadata$validation: tutti C1-C8 PASS o solo WARN
#   - C2, C3, C5: ERROR = 0 (hard fail conditions)
#   - Run con stesso seed → stesso n_cells (riproducibilità)
#   - Run con seed diverso → n_cells in ±5% (stabilità densità)
```

**`test_cell_layer_edge_cases.R`** — robustezza:

| Caso | Comportamento atteso |
|---|---|
| Cluster con `area_um2 < min_region_area_um2` | regione esclusa, log info |
| `region_composition` con frazioni che non sommano 1.0 | normalizzazione + warning |
| `cell_type` in composition ma assente dal catalogo | hard error con messaggio puntuale |
| Tessuto enorme (4000×4000 px) | non crash, < 5 min, RAM < 8 GB (v1.6, S1.5 P-8: rispettato in 5 tessuti; linfonodo 310 s e 8.50 GB, BL-144) |
| Tessuto con un solo cluster | funziona, una sola regione |
| `corner_smoothing = 1` con territori piccoli | nessun poligono degenerato (area > 0) |
| Composizione con `fraction = 0` per un tipo | quel tipo non viene piazzato |
| `pixel_size_um` omesso | hard error con messaggio puntuale (C-A) |
| `pixel_size_um = "estimate"` su immagine con scala nota (1, 2, 5, 10 µm/px) | scala recuperata entro ±20%; `pixel_size_source = "estimated"` + WARN (controprova S1.5; v1.6: smentita, la stima richiede `allow_estimate = TRUE`) |
| `pixel_size_um` sbagliato di 10× | C8 scatta (WARN/ERROR). **v1.6 (S1.5): non realizzabile** con densità assolute (vedi C8, §6.1); errore della pre-registrazione (DEC-011, Luca, 2026-10-21); la scala viene dall'utente o dai metadati |

### 6.3 Visualizzazione (`R/testing/visualize_cell_layer.R`)

Pannello a 4 sub-plot, sfondo bianco (vincolo CLAUDE.md):

1. **Regioni del tessuto** — `region_polygons` colorate per `cluster_id`
2. **Centroidi e tipi cellulari** — punti colorati per `cell_type`, `region_polygons` di sfondo trasparente
3. **Tassellazione Voronoi** — `cell_territories` con bordi neri, colorati per `cell_type`; nuclei sovrapposti in grigio scuro
4. **Densità raggiunta** — choropleth: `region_polygons` colorato per `achieved_density / target_density × 100%` (validazione visiva di C1)

Output: `R/testing/cell_layer_panel.png` (1600×1200, `bg = "white"` in `ggsave()`).

Utility funzione esposta:
```r
plot_cell_layer(cell_layer, mode = c("voronoi", "types", "density", "regions"))
```

---

## 7. Hand-off a Step 2

### 7.1 Strutture dati consumate da Step 2

| Data | Uso in Step 2 | Sezione del REPORT |
|---|---|---|
| `cell_layer$cell_df` | substrato per generazione profili a livello cellula | 3.2.B |
| `cell_layer$cell_territories` | input geometrico per `st_intersection(spot_grid, ...)` → frazioni di overlap | 3.2.A |
| `cell_layer$cell_nuclei` | distinzione spot intra-nucleo vs intra-citoplasma → profilo nucleare diverso | 3.2.D |
| `cell_layer$region_polygons` | maschera per gli "spot extracellulari macro" (solo ambient RNA) | 3.2.C |
| `cell_df$nucleus_area` / `cytoplasm_area` | variazione intra-cellulare ponderata per area | 3.2.D |
| `cell_layer$metadata$pixel_size_um` (+ `pixel_size_source`) | convenzione spaziale per costruire la spot grid 2 µm **in proporzione al tessuto** (C-A) | 3.2 |
| `cell_layer$cell_territories` + `cell_nuclei` | **sorgenti del kernel di diffusione laterale degli mRNA** (C-B, v1.1): in Step 2 i conteggi cellulari vengono diffusi con un kernel gaussiano isotropo (`diffusion_sigma_um`, ordine di grandezza atteso di pochi µm, **da stimare** in R5) prima della somma sugli spot | 4.3 |

### 7.2 Vincoli su Step 1 (perché Step 2 funzioni)

- `cell_territories` devono essere **VALID** (`sf::st_is_valid()`) e senza self-intersection — `tessellate_voronoi()` chiamerà `sf::st_make_valid()` come passo finale di sicurezza.
- `cell_nuclei` devono essere **CONTAINED** nei rispettivi territori (già garantito da check C3 di 6.1).
- `metadata$pixel_size_um` deve essere sempre presente e numerico (C-A): Step 2 non ha un proprio default di scala.
- Ordine di `cell_id` consistente fra `cell_df`, `cell_territories`, `cell_nuclei` (1..N stesso ordine — già garantito dall'orchestrator).
- I poligoni `cell_territories` saranno restituiti come `sfc` con bbox calcolato (default di `sf`), così Step 2 può sfruttare l'R-tree di GEOS per `st_intersects` in O(log n) sulla spot grid HD.

### 7.3 Pre-validazioni offerte gratis a Step 2

- **Copertura del tessuto**: % area tessuto coperta da `cell_territories` (≈100% sopra threshold, 0% sotto threshold) — diagnostica per il rate di "spot vuoti" attesi
- **Distribuzione `eq_radius` per `cell_type`** — per calibrare i profili NB cell-type-specific (06d) e i parametri di library size cell-type-specific (06e)
- **Conteggi cellule per `region_id`** — per testare ergonomia del mapping `region_id → spot_ids` in Step 2

---

### 7.4 Limiti dichiarati del livello cellulare (v1.8, S1.6)

Ciò che il livello cellulare non riproduce, misurato contro i riferimenti veri, e da tenere presente quando si testano i metodi di analisi sui dati simulati:

- **Nuclei a contatto in eccesso** nei tessuti densi: tumore del colon 0.69 contro al più 0.57 del conteggio manuale, linfonodo 0.83 contro al più 0.75, al limite nello stroma del tumore 0.36 contro al più 0.28 (DEC-029; nessun disegno del nucleo provato li porta nel vero; nel tumore del colon il riferimento manuale è il meno sicuro).
- **Territori troppo uniformi** con la semina uniforme dei preset (variabilità dell'area dei territori e numero di vicini sotto il vero in intestino, tumore del colon, linfonodo; S1.5b-S1.5e): la semina a compartimenti con i vuoti è un'opzione che ne corregge una parte ma non è il predefinito (regola B).
- **Cuore**: territori isotropi, cellule vere allungate (controllo negativo previsto, §4.6, R3).
- **Densità dei preset** riferita al tessuto pieno di cellule del conteggio manuale, non al tessuto con i suoi vuoti (BL-167, da decidere in uno step dedicato).

## 8. Estensioni future (Step 1.5)
**[v1.9, S1.6b: elenco di v1.1 (in S0 è cambiata solo la voce del territorio anisotropo). Implementati poi: nucleo non centrato sul baricentro del territorio (centro nel generatore con il modello «clip», `nucleus_center`, §5.4), nucleo a ellisse come opzione (§5.4, v1.4), semina a compartimenti con vuoti e densità eterogenea come opzione (§5.2, v1.8); il territorio anisotropo non è implementato e il cuore resta il controllo negativo (§7.4).]**

- Nucleo offset rispetto al centroide territorio
- Nucleo irregolare (ellisse o blob organico)
- Forma anisotropa del territorio (Voronoi pesato per modellare epitelio colonnare / muscolo) — motivato quantitativamente dal fallimento atteso dei check B su A6 (`"columnar_muscle"`)
- Eterogeneità ricca: gradienti di transizione + cluster localizzati di infiltrati
- Ulteriori tissue presets (richiede curatela e validazione biologica)

---

## 9. Storico revisioni

| Data | Cosa | Da chi |
|---|---|---|
| 2026-04-29 | Decisioni 1-7 prese durante prima sessione brainstorming | Luca + Claude |
| 2026-04-29 | Decisioni 1-3 RIVISTE post-feedback Daniele (Voronoi + nucleo grande, composizione esplicita, catalogo disaccoppiato) | Daniele review |
| 2026-04-29 | Decisioni 8-9 aggiunte (tissue presets, corner smoothing) | Daniele review |
| 2026-04-29 | Sezioni 1-5 approvate e consolidate in questo doc | Luca |
| 2026-04-29 | Sezioni 6-7 finalizzate (validazione, testing, visualizzazione, hand-off Step 2) | Luca |
| 2026-04-29 | v1.0 inviata a Daniele per revisione (commit `7836829`) | Luca |
| 2026-09-18 | **v1.1 (S0)** — C-A: `pixel_size_um` obbligatorio, helper `estimate_pixel_size()`, `metadata$pixel_size_source`, check C8, casi limite sulla scala (decisione 10). C-B: diffusione laterale mRNA come requisito di Step 2, riga in 7.1, vincolo in 7.2 (decisione 11). Pannello A1–A6 → 6 preset in 4.6/4.7 (decisione 8 rivista); test su tutti i preset. Corretto refuso nella riga C6 (nome campo `region_id`). | Daniele (commenti) + Luca + Claude |
| 2026-10-08 | **v1.2 (S1.2, parziale)** — §5.2: densità di miscela armonica (D3); RSA a n fissato invece di Bridson (Bridson satura a ≈ 5× la densità target con d = ⅔·eq_radius); arrotondamento stocastico; distanza minima globale d_ij = (d_i+d_j)/2; `min_dist_um` opzionale. §5.1 passo 4 (BL-043) e stride (BL-042) decisi, da implementare | Luca + Claude |
| 2026-10-08 | **v1.3 (S1.3)** — §5.3: Voronoi **per regione** (D-S1.3.1; il testo «globale + ritaglio sulla regione del centroide» lascia scoperta in mediana il 9.4 % delle regioni piccole, CP-1); frammenti orfani al territorio con il confine condiviso più lungo (D-S1.3.2); Chaikin con `n_iter = max(1, round(3·cs))` intersecato con il territorio originale (D-S1.3.3); firma `tessellate_voronoi(centroids, regions, corner_smoothing, keep_tiles, verbose)`; decisione 6: motore GEOS invece di deldir (confronto C-10, `results/S1.3/`) | Luca + Claude |
| 2026-10-15 | **v1.4 (S1.4)** — §5.4 riscritta: modello di taglia a due parametri (β, σ) che contiene il design (D-S1.4.1); clamp in tutti i modi sul raggio circoscritto, bisezione per l'ellisse (D-S1.4.2); centro POI o generatore; §4.4 colonne opzionali; §4.2 `nucleus_cx/cy`. Contenimento in conflitto con i dati nei tessuti densi (PH-1, BL-098) | Luca (scelte in chat) + S1.4 |
| 2026-10-20 | **v1.5 (S1.4d)** — §5.4: modello «clip» (nucleo ritagliato nel proprio territorio e portato all'area del modello di taglia, K-2 di S1.4c, decisione di Luca del 2026-10-16) predefinito, «scale» (v1.4) invariato come opzione; colonna `nucleus_nc_max` (φ) in §4.4; operatore di osservazione dei check B (BL-112) e aggregazione di S-1 sulle finestre a peso uguale (BL-131) | Luca (decisioni S1.4c, R2c) + S1.4d |
| 2026-10-21 | **v1.6 (S1.5)** — §4.7 firma implementata di `place_cells()` (`stride`, `small_components`, `size_conversion`, `keep_tiles`, `allow_estimate`), `validate_cell_layer()` (C1–C9), `overlay_prec()`, `plot_cell_layer()`; §4.2 `x`, `y` = generatori (D-9); §5.4 tetti φ di A4 e A6 da B-143 e B-154 (revisione FND-11); decisione 10 e §3: stima della scala dall'immagine fuori dall'uso raccomandato (CP-1), fonti ammesse (BL-016); §6.1 C8 e §6.2: una scala sbagliata non è rilevabile da C8 con densità assolute (D-3, DEC-011); prestazioni del linfonodo oltre il limite (BL-144) | S1.5 (agente; C8 e caso limite decisi da Luca, DEC-011) |
| 2026-10-23 | **v1.7 (S1.5c)**: §5.2 semina a compartimenti per tipo come opzione (`spatial = "compartments"`: campi latenti per tipo, densità del tipo puro, quote; distanza minima per nucleo); il predefinito resta la semina uniforme (DEC-020); §8 «eterogeneità ricca» in parte implementata, spazi vuoti dentro il tessuto ancora da fare (BL-157) **[v1.9, S1.6b: il commit di v1.7 (`915d044`) aggiunge solo questa riga; il testo di §5.2 sulla semina a compartimenti è entrato con v1.8, §8 non è più cambiato dopo v1.1]** | geo-agent (S1.5c) |
| 2026-10-29 | **v1.8 (S1.6)**: §5.2 stato della semina (cieca al tipo per i preset, DEC-019; compartimenti senza quote, DEC-025; limite esatto del diradamento, DEC-030); §4.7 firme e opzioni di sessione; §4.3 media armonica; §4.4 tipi a densità 0; §4.6 preset dai dati, composizione unica per cluster, condizione di `provisional` soddisfatta; §5.4 tipo MULTIPOLYGON (BL-173); §6.1 C2, C4, C9 come implementati, C7 con la miscela per la semina cieca (BL-156, decisione di Luca, DEC-031); §6.2 test per step; §7.4 limiti dichiarati (DEC-029). Testo allineato al codice dopo la revisione avversariale di S1.6 (RC-06...RC-14, RC-21) | geo-agent (S1.6); C7: Luca |
| 2026-10-30 | **v1.9 (S1.6b)**: solo testo, nessun comportamento nuovo. Intestazione e riga di stato; §2 decisione 1 (taglia del nucleo dal modello di §5.4); §5.1 come implementato in S1.1 (union-find, contorni sui lati dei pixel) e componenti piccole unite per predefinito da `place_cells()` (BL-043); §4.6 rimando sulla colonna `provisional`; §5.2 passo 5 (distanza minima dei preset = d*, DEC-009); §8 (cosa è stato implementato); nota sulla riga di v1.7. Discrepanze trovate dalla mappa design-codice per il dossier di Daniele | geo-agent (S1.6b) |
