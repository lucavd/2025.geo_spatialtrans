# Guida Completa alla Validazione Biologica per Simulazioni Spatial Transcriptomics

## Indice
1. [Panoramica del Sistema](#panoramica-del-sistema)
2. [Script di Validazione Biologica](#script-di-validazione-biologica)
3. [Parametri Chiave e loro Effetti](#parametri-chiave-e-loro-effetti)
4. [Risoluzione Problemi Comuni](#risoluzione-problemi-comuni)
5. [Ottimizzazione per Realismo Biologico](#ottimizzazione-per-realismo-biologico)
6. [Metriche di Validazione](#metriche-di-validazione)
7. [Troubleshooting Avanzato](#troubleshooting-avanzato)

---

## Panoramica del Sistema

### Architettura della Simulazione
Il sistema è composto da moduli interconnessi che generano dati spatial transcriptomics realistici:

```
Input Image → Clustering → Grid Sampling → Expression Generation → Validation
     ↓              ↓            ↓                ↓                  ↓
  Tissue Map → Cell Types → Spatial Coords → Gene Expression → Quality Report
```

### File Principali
- **`R/run_full_simulation_simple.R`** - Script principale di simulazione
- **`R/biological_validation_report.R`** - Sistema di validazione biologica
- **`R/functions/06*_expression*.R`** - Moduli di generazione espressione
- **`R/functions/07*_simulation*.R`** - Pipeline di simulazione
- **`R/functions/08_validation_plots.R`** - Plot di validazione

---

## Script di Validazione Biologica

### Utilizzo Base

```bash
# Esegui simulazione
Rscript R/run_full_simulation_simple.R

# Esegui validazione biologica
Rscript R/biological_validation_report.R
```

```r
# In R direttamente
source("R/biological_validation_report.R")
validation_results <- generate_biological_validation_report("results/simple_simulation.rds")
```

### Output della Validazione

Lo script genera automaticamente:

1. **`plots/biological_validation/biological_validation_report.md`**
   - Report testuale dettagliato
   - Sommario esecutivo con ✓/✗ per ogni metrica
   - Raccomandazioni specifiche

2. **`plots/biological_validation/biological_validation_summary.png`**
   - Grafico a barre dei punteggi di validazione
   - Soglia di accettabilità (80%)
   - Categorizzazione per tipo di metrica

3. **`plots/biological_validation/biological_validation_umi_distribution.png`**
   - Istogramma distribuzione UMI per cella
   - Linee di riferimento per range biologico (1000-15000)

4. **`plots/biological_validation/biological_validation_detailed_results.rds`**
   - Risultati completi per analisi avanzate

### Metriche Validate

| Metrica | Range Target | Descrizione |
|---------|--------------|-------------|
| **UMI Totali per Cella** | 1,000 - 15,000 | Realistic sequencing depth |
| **Espressione Massima** | < 10,000 | Evita geni irrealisticamente alti |
| **Mediana Espressione** | 5 - 100 | Baseline expression appropriato |
| **Coerenza Spaziale** | Cluster compatti | Spatial organization realistico |
| **Specificità Marker** | Fold-change > 1.5 | Marker genes distinguibili |

---

## Parametri Chiave e loro Effetti

### 1. Library Size Parameters (`library_size_params`)

```r
library_size_params = list(
  mean_library_size = 500,           # Media UMI per cella
  library_size_cv = 0.15,            # Coefficiente di variazione
  spatial_effect_on_library = 0.05,  # Effetto spaziale (moltiplicativo)
  cell_type_effect = FALSE           # Effetto cell-type (moltiplicativo)
)
```

**Effetti:**
- `mean_library_size`: Controlla UMI totali base
- `library_size_cv`: Variabilità tra celle
- `spatial_effect_on_library`: Gradiente spaziale (usa `exp()`, pericoloso!)
- `cell_type_effect`: Differenze tra tipi cellulari

**Valori Storici Problematici:**
- Default originale: `mean_library_size = 10000` → UMI ~400k (troppo alto)
- Con effetti moltiplicativi attivi: esplosione dei valori

### 2. Baseline Expression (`06b_expression_baseline.R`)

```r
mu <- rep(-1.2, n_genes)  # baseline log-expression
```

**Effetti:**
- `mu = 2`: exp(2) ≈ 7.4 → UMI ~400k per simulazione
- `mu = -1`: exp(-1) ≈ 0.37 → UMI ~20k per simulazione  
- `mu = -1.2`: exp(-1.2) ≈ 0.30 → UMI ~16k per simulazione
- `mu = -1.5`: exp(-1.5) ≈ 0.22 → UMI ~12k per simulazione

**Problema Mediana Bassa:**
Baseline uniforme → molti geni con espressione simile → mediana compressa

### 3. Dropout Parameters (`difficulty_level = "medium"`)

```r
dropout_params = list(
  dropout_range = c(0.2, 0.5),              # 20-50% dropout rate
  dispersion_range = c(2.0, 1.0),           # Negative binomial dispersion
  expression_dependent_dropout = TRUE,       # Dropout correlato a espressione
  dropout_curve_midpoint = 0.5,
  dropout_curve_steepness = 5
)
```

**Problemi Identificati:**
- Dropout 20-50% troppo aggressivo per spatial transcriptomics
- Dispersione bassa (2.0, 1.0) → più zeri dalla negative binomial
- Expression-dependent dropout amplifica il problema

### 4. Marker Parameters

```r
marker_params = list(
  marker_genes_per_type = 7,          # Geni marker per tipo cellulare
  marker_expression_fold = 1.2,       # Fold-change per marker (medium)
  marker_overlap_fold = 0.2           # Overlap tra tipi adiacenti
)
```

---

## Risoluzione Problemi Comuni

### Problema 1: UMI Troppo Alti (>50k per cella)

**Cause:**
1. `mean_library_size` troppo alto
2. Effetti moltiplicativi (spatial, cell-type) troppo forti
3. Baseline expression troppo alto

**Soluzioni:**
```r
# 1. Ridurre library size
library_size_params$mean_library_size <- 300  # da 500

# 2. Disabilitare effetti moltiplicativi
library_size_params$spatial_effect_on_library <- 0
library_size_params$cell_type_effect <- FALSE

# 3. Ridurre baseline
mu <- rep(-1.8, n_genes)  # da -1.2
```

### Problema 2: Mediana Espressione Troppo Bassa (<1)

**Cause:**
1. Baseline troppo basso
2. Dropout troppo aggressivo
3. Dispersione negative binomial troppo bassa
4. Tutti i geni partono dallo stesso valore

**Soluzioni:**

#### A. Diversificare Baseline Expression
```r
# In 06b_expression_baseline.R
# Invece di: mu <- rep(-1.2, n_genes)

# Distribuizione più realistica
mu <- rnorm(n_genes, mean = -1.0, sd = 0.3)
mu <- pmax(mu, -2.5)  # Floor minimo
mu <- pmin(mu, 0.5)   # Ceiling massimo
```

#### B. Ridurre Dropout
```r
# In 07b_difficulty_setup.R per "medium"
dropout_range = c(0.1, 0.3)  # da c(0.2, 0.5)
```

#### C. Aumentare Dispersione NB
```r
# In 07b_difficulty_setup.R per "medium"  
dispersion_range = c(5.0, 3.0)  # da c(2.0, 1.0)
```

### Problema 3: Marker Genes Non Identificabili

**Cause:**
1. `marker_expression_fold` troppo basso
2. Overlap tra tipi cellulari troppo alto
3. Rumore spaziale troppo forte

**Soluzioni:**
```r
# Aumentare fold-change marker
marker_expression_fold = 2.0  # da 1.2

# Ridurre overlap  
marker_overlap_fold = 0.1  # da 0.2

# Ridurre rumore spaziale
spatial_noise_intensity = 0.5  # da 1.0
```

---

## Ottimizzazione per Realismo Biologico

### Configurazione Raccomandata

```r
# In run_full_simulation_simple.R
simulate_spatial_transcriptomics(
  # ... altri parametri ...
  difficulty_level = "custom",  # Per controllo completo
  
  library_size_params = list(
    mean_library_size = 400,
    library_size_cv = 0.12,
    spatial_effect_on_library = 0.02,
    cell_type_effect = FALSE
  ),
  
  dropout_params = list(
    dropout_range = c(0.05, 0.25),
    dispersion_range = c(8.0, 5.0),
    expression_dependent_dropout = TRUE,
    dropout_curve_midpoint = 0.3,
    dropout_curve_steepness = 3
  ),
  
  marker_params = list(
    marker_genes_per_type = 8,
    marker_expression_fold = 1.8,
    marker_overlap_fold = 0.1
  )
)
```

```r
# In 06b_expression_baseline.R - baseline diversificato
set.seed(random_seed + k)  # Per riproducibilità per tipo
mu_base <- -1.1
mu_sd <- 0.25

# Distribuizione log-normale troncata per baseline
mu <- rlnorm(n_genes, meanlog = mu_base, sdlog = mu_sd)
mu <- log(mu)  # Torna a log-space
mu <- pmax(mu, -2.0)  # Floor
mu <- pmin(mu, -0.5)  # Ceiling
```

### Pipeline di Ottimizzazione

1. **Test Baseline**: Parti da parametri conservativi
2. **Valida UMI**: Assicurati che UMI totali siano 5-20k
3. **Valida Mediana**: Target 3-15 per mediana espressione
4. **Valida Marker**: Fold-change marker > 1.5
5. **Valida Spaziale**: Coherenza cluster appropriata

### Script di Test Rapido

```r
# Test rapido parametri
test_parameters <- function(mean_lib, baseline, dropout_max) {
  # Modifica parametri
  # Esegui simulazione ridotta (n_genes=100, piccola griglia)
  # Valida solo UMI e mediana
  # Ritorna TRUE se nel range target
}

# Grid search
best_params <- NULL
for (lib in c(300, 400, 500)) {
  for (base in c(-1.3, -1.1, -0.9)) {
    for (drop in c(0.2, 0.3, 0.4)) {
      if (test_parameters(lib, base, drop)) {
        best_params <- list(lib=lib, base=base, drop=drop)
        break
      }
    }
  }
}
```

---

## Metriche di Validazione

### Distribuzione UMI Ideale

```
Spatial Transcriptomics (es. Visium HD):
- Mediana: 8,000-12,000 UMI
- Range: 3,000-25,000 UMI  
- 95% celle: 1,000-15,000 UMI

Single Cell (riferimento):
- Mediana: 3,000-8,000 UMI
- Range: 500-15,000 UMI
```

### Gene Expression Guidelines

```
Per Gene Expression Level:
- Housekeeping genes: 10-100 counts
- Marker genes: 50-500 counts  
- Highly expressed: 100-2000 counts
- Maximum realistic: <10,000 counts

Per Cell:
- Median non-zero expression: 5-50
- Fraction zeros: 60-90% (spatial transcriptomics)
- Genes detected per cell: 1000-4000
```

### Spatial Coherence Metrics

```
Cluster Spatial Properties:
- Intra-cluster distance: <50% of tissue diameter
- Inter-cluster separation: >10% of tissue diameter
- Cluster compactness: Area/Perimeter ratio > 0.3
- Spatial autocorrelation: Moran's I > 0.3 for marker genes
```

---

## Troubleshooting Avanzato

### Debug Workflow

1. **Abilita Debug Mode**:
```r
# In 06e_library_size.R (già aggiunto)
cat("DEBUG - Library size stats: mean =", round(mean(library_size)), "\n")

# Aggiungi in 06j_expression_generation.R
cat("DEBUG - Gene", g, "raw_counts range:", range(raw_counts), "\n")
cat("DEBUG - Gene", g, "final_counts range:", range(scaled_counts), "\n")
```

2. **Analisi Step-by-Step**:
```r
# Carica risultati e analizza componenti
sim_results <- readRDS("results/simple_simulation.rds")

# UMI distribution
umi_per_cell <- colSums(sim_results$expression)
summary(umi_per_cell)

# Gene expression stats  
gene_means <- rowMeans(sim_results$expression)
gene_medians <- apply(sim_results$expression, 1, median)
summary(gene_means)
summary(gene_medians)

# Zero fraction
zero_fraction <- mean(sim_results$expression == 0)
cat("Zero fraction:", zero_fraction, "\n")
```

3. **Component Isolation**:
```r
# Test solo library size generation
test_lib_sizes <- generate_library_sizes(
  cell_df, 
  library_size_params = list(mean_library_size = 500, ...)
)
summary(test_lib_sizes)

# Test solo baseline expression
test_baseline <- generate_baseline_expression(100, 5, marker_params)
sapply(test_baseline, summary)
```

### Performance Optimization

```r
# Per dataset grandi, riduci per test
simulate_spatial_transcriptomics(
  n_genes = 200,           # da 2000
  grid_resolution = 50,    # da 30  
  # ... test rapido
)

# Usa subset per validazione
validation_subset <- sim_results
validation_subset$expression <- sim_results$expression[1:100, 1:1000]
validation_subset$coordinates <- sim_results$coordinates[1:1000, ]
validation_subset$intensity_cluster <- sim_results$intensity_cluster[1:1000]
```

### File Locations per Debug

```
File di configurazione:
- R/functions/06a_expression_params.R (default parameters)
- R/functions/07b_difficulty_setup.R (difficulty levels)

File di generazione:
- R/functions/06b_expression_baseline.R (baseline μ)
- R/functions/06e_library_size.R (library size generation)
- R/functions/06j_expression_generation.R (final count generation)

File di output:
- results/simple_simulation.rds (simulation results)
- plots/biological_validation/ (validation reports)
- plots/validazione_simple/ (basic validation plots)
```

---

## Checklist Finale

### Pre-Simulazione
- [ ] Parametri library_size appropriati (<1000 mean)
- [ ] Baseline expression bilanciato (-1.5 a -0.8)
- [ ] Dropout non eccessivo (<40% max)
- [ ] Effetti moltiplicativi controllati

### Post-Simulazione  
- [ ] UMI mediana 5-20k
- [ ] Espressione massima <10k
- [ ] Mediana espressione 3-50
- [ ] >80% celle in range realistico
- [ ] Marker genes identificabili (fold-change >1.5)

### Validazione
- [ ] Report biologico generato senza errori
- [ ] Grafici di distribuzione ragionevoli  
- [ ] Coerenza spaziale appropriata
- [ ] Specificità marker sufficiente

**Target Finale: Tutti i parametri ✓ REALISTICO nel report di validazione biologica.**