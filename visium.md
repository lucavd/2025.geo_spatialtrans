# Implementazione di una Simulazione Avanzata di Dati Visium HD

## Introduzione

Questo documento descrive le modifiche apportate al framework di simulazione per renderlo più aderente alle caratteristiche dei dati di trascrittomica spaziale Visium HD di 10x Genomics. Visium HD rappresenta un'evoluzione significativa rispetto alla tecnologia Visium originale, con una risoluzione notevolmente superiore (2μm x 2μm per spot) che consente di ottenere dati più vicini alla risoluzione di singola cellula.

Le modifiche implementate mirano a:
1. Riprodurre accuratamente la struttura spaziale e la risoluzione di Visium HD
2. Modellare realisticamente i pattern spaziali di espressione genica
3. Simulare caratteristiche tecniche come dropout e dimensione della libreria
4. Fornire metriche di validazione per confrontare i dati simulati con quelli reali

## 1. Struttura Spaziale e Risoluzione

### Implementazione della Griglia Regolare ad Alta Risoluzione

**Problema affrontato:** La versione precedente del simulatore campionava casualmente punti dall'immagine, perdendo la struttura regolare e la risoluzione specifica delle tecnologie Visium HD.

**Soluzione implementata:**
- Introduzione della modalità `grid_mode` che crea una griglia regolare
- Parametrizzazione attraverso `grid_resolution` (2μm per Visium HD) e `grid_spacing`
- Proiezione dell'immagine di input sulla griglia per mantenere la corrispondenza con la morfologia tissutale

```r
# Calcola le dimensioni della griglia in base all'immagine
img_width_um <- w
img_height_um <- h

# Calcola il numero di bin necessari
n_bins_x <- ceiling(img_width_um / grid_resolution)
n_bins_y <- ceiling(img_height_um / grid_resolution)

# Crea le coordinate della griglia
x_coords <- seq(1, img_width_um, by = grid_resolution + grid_spacing)
y_coords <- seq(1, img_height_um, by = grid_resolution + grid_spacing)

# Crea un dataframe con tutte le coordinate possibili
grid_points <- expand.grid(x = x_coords, y = y_coords)
```

**Fondamento teorico:** Visium HD utilizza una griglia regolare di spots adiacenti a 2μm x 2μm, consentendo la cattura dell'RNA con una risoluzione prossima a quella cellulare. Questa struttura regolare, a differenza del campionamento casuale, permette di:
1. Modellare correttamente effetti di autocorrelazione spaziale
2. Riflettere la continuità tissutale presente nei dati reali
3. Facilitare l'integrazione con altre modalità di imaging

## 2. Modellazione dei Pattern Spaziali

### 2.1 Implementazione di Gradienti tra Regioni

**Problema affrontato:** Il modello originale creava confini netti tra cluster, mentre nei tessuti reali esistono spesso zone di transizione graduale.

**Soluzione implementata:**
- Parametro `gradient_regions` che attiva gradienti ai confini tra cluster
- Calcolo della distanza dal confine per ogni punto della griglia
- Interpolazione dell'espressione genica basata sulla distanza dal confine

```r
if (spatial_params$gradient_regions && exists("boundary_dist", where = cell_df) && 
    exists("in_gradient", where = cell_df)) {
  # Usa l'approccio gradiente per i confini
  for (i in 1:N) {
    if (cell_df$in_gradient[i]) {
      # Se questo punto è in zona di gradiente
      orig_cluster <- as.integer(cluster_labels[i])
      
      # Trova cluster vicini
      nearby_idx <- which(dist_mat[i, ] < spatial_params$gradient_width * grid_resolution)
      nearby_clusters <- unique(as.integer(cluster_labels[nearby_idx]))
      nearby_clusters <- nearby_clusters[nearby_clusters != orig_cluster]
      
      if (length(nearby_clusters) > 0) {
        # Peso basato sulla distanza dal bordo
        gradient_weight <- (1 - cell_df$boundary_dist[i])^2
        
        # Mix dell'espressione
        other_expr <- mean(sapply(nearby_clusters, function(k) mean_expression_list[[k]][g]))
        base_expr[i] <- base_expr[i] * (1 - gradient_weight) + other_expr * gradient_weight
      }
    }
  }
}
```

**Fondamento teorico:** I gradienti di espressione sono un fenomeno biologico comune nei tessuti reali, dove:
1. Zone di transizione tra tipi cellulari mostrano profili di espressione intermedi
2. Gradienti morfogenetici influenzano l'espressione genica in modo continuo
3. Segnali paracrini creano pattern di espressione che decadono con la distanza

Studi come quelli di Edsgärd et al. (2018) e Svensson et al. (2018) hanno dimostrato che i gradienti sono un elemento cruciale per la modellazione realistica dei dati spaziali.

### 2.2 Modellazione dell'Autocorrelazione Spaziale

**Problema affrontato:** Il modello originale usava un solo approccio per la correlazione spaziale, senza quantificare o validare il grado di autocorrelazione.

**Soluzione implementata:**
- Supporto per due metodi di simulazione: Gaussian Random Fields (GRF) e Conditional Autoregressive (CAR)
- Calcolo dell'indice di Moran I per quantificare l'autocorrelazione spaziale
- Visualizzazione dei valori di autocorrelazione per validazione

```r
# GRF implementation
if (correlation_method == "grf") {
  # Converti cell_df in oggetto spatial
  sp_df <- cell_df
  coordinates(sp_df) <- ~ x + y
  
  # Crea un modello gstat per il GP con parametri personalizzati
  gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                beta = 0, model = vgm(psill=spatial_params$spatial_noise_intensity, 
                                     range=spatial_params$spatial_range, 
                                     model="Exp"), 
                nmax=10)
  
  # Genera il noise spaziale
  gp_noise <- predict(gp_sim, newdata = sp_df, nsim = 1)$sim1
  
  # Normalizza il noise
  gp_noise <- scale(gp_noise)
}

# CAR implementation
else if (correlation_method == "car") {
  # Crea una griglia di vicinanza
  coords_matrix <- as.matrix(cell_df[, c("x", "y")])
  nb <- dnearneigh(coords_matrix, 0, spatial_params$spatial_range)
  
  # Creiamo la matrice dei pesi spaziali
  lw <- nb2listw(nb, style = "W")
  
  # Simula un processo spaziale autocorrelato
  gp_noise <- spam.spGauss(lw, sigma2 = spatial_params$spatial_noise_intensity, n = 1)
}

# Calcolo dell'indice di Moran I per validazione
if (requireNamespace("spdep", quietly = TRUE)) {
  coords_matrix <- as.matrix(cell_df[, c("x", "y")])
  nb <- spdep::dnearneigh(coords_matrix, 0, spatial_params$spatial_range)
  lw <- spdep::nb2listw(nb, style = "W")
  
  # Calcola Moran's I per i primi 10 geni
  moran_genes <- min(10, ncol(expression_data))
  moran_i <- numeric(moran_genes)
  
  for (g in 1:moran_genes) {
    if (!all(expression_data[, g] == 0)) {
      moran_i[g] <- spdep::moran.test(expression_data[, g], lw)$estimate[1]
    }
  }
}
```

**Fondamento teorico:** 
1. **Gaussian Random Fields (GRF):** Modellano fenomeni spaziali attraverso funzioni di covarianza esponenziale che rispecchiano il decadimento della correlazione con la distanza, similmente ai processi di diffusione molecolare nei tessuti.

2. **Modelli Conditional Autoregressive (CAR):** Particolarmente adatti per dati su griglia regolare, modellano l'espressione in un punto come dipendente dai valori nei punti circostanti. Matematicamente, per ogni posizione i:
   ```
   X_i | X_{-i} ~ N(ρ * Σ_j w_{ij} X_j, τ²)
   ```
   dove X_{-i} sono i valori a tutte le altre posizioni, w_{ij} sono pesi spaziali.

3. **Indice di Moran I:** Metrica fondamentale in statistica spaziale che quantifica il grado di autocorrelazione:
   ```
   I = (n/S0) * Σ_i Σ_j w_{ij}(x_i - x̄)(x_j - x̄) / Σ_i(x_i - x̄)²
   ```
   I valori variano da -1 (dispersione perfetta) a +1 (clustering perfetto), con 0 che indica casualità spaziale.

## 3. Simulazione dell'Espressione Genica

### 3.1 Modellazione della Dimensione della Libreria

**Problema affrontato:** Il modello originale non considerava la variabilità nella dimensione della libreria (profondità di sequenziamento) tra gli spots, un fattore critico nei dati reali.

**Soluzione implementata:**
- Simulazione delle dimensioni delle librerie con distribuzione log-normale
- Effetto spaziale sulla dimensione della libreria
- Effetto del tipo cellulare sulla dimensione della libreria
- Scaling dei conteggi di espressione in base alla dimensione della libreria

```r
# Genera dimensioni libreria con distribuzione log-normale
library_size_sd <- library_size_params$mean_library_size * library_size_params$library_size_cv
log_mean <- log(library_size_params$mean_library_size^2 / 
               sqrt(library_size_params$mean_library_size^2 + library_size_sd^2))
log_sd <- sqrt(log(1 + (library_size_sd^2 / library_size_params$mean_library_size^2)))

library_size <- rlnorm(N, meanlog = log_mean, sdlog = log_sd)

# Aggiungi effetto spaziale sulla dimensione libreria
if (library_size_params$spatial_effect_on_library > 0) {
  # Crea un GP con correlazione spaziale a raggio ampio
  sp_df_lib <- cell_df
  coordinates(sp_df_lib) <- ~ x + y
  
  lib_gp <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                 beta = 0, model = vgm(psill = 1.0, 
                                      range = spatial_params$spatial_range * 1.5, 
                                      model = "Exp"), 
                 nmax = 15)
  
  lib_noise <- predict(lib_gp, newdata = sp_df_lib, nsim = 1)$sim1
  
  # Applica alla dimensione libreria (effetto moltiplicativo)
  lib_noise <- scale(lib_noise)
  lib_effect <- library_size_params$spatial_effect_on_library * lib_noise
  library_size <- library_size * exp(lib_effect)
}

# Aggiungi effetto del tipo cellulare sulla dimensione della libreria
if (library_size_params$cell_type_effect) {
  # Diversi tipi cellulari hanno diversi contenuti di RNA
  cell_type_effect <- numeric(N)
  
  # Crea effetti diversi per diversi tipi cellulari
  type_effects <- rnorm(k_cell_types, mean = 0, sd = 0.2)
  
  # Assegna effetto in base al tipo cellulare
  for (k in 1:k_cell_types) {
    cell_type_effect[cluster_labels == k] <- type_effects[k]
  }
  
  # Applica effetto moltiplicativo
  library_size <- library_size * exp(cell_type_effect)
}

# Applica l'effetto della dimensione della libreria ai conteggi
scaled_counts <- raw_counts * (library_size / mean(library_size))
expression_data[, g] <- round(scaled_counts)
```

**Fondamento teorico:** La distribuzione della dimensione della libreria nei dati di trascrittomica segue tipicamente una distribuzione log-normale (Lytal et al. 2020). Nei dati Visium reali, la dimensione della libreria mostra pattern complessi dovuti a:
1. Variazione dell'efficienza di cattura RNA in diverse regioni tissutali
2. Differenze nella densità cellulare tra regioni
3. Contenuto di RNA intrinsecamente diverso tra tipi cellulari (es. cellule secretorie hanno più RNA di cellule di supporto)
4. Effetti tecnici durante la preparazione delle librerie

La modellazione realistica di questa variazione, inclusa la dipendenza dal tipo cellulare, è fondamentale per simulare correttamente la variabilità tecnica e biologica nei dati.

### 3.2 Dropout Dipendente dall'Espressione

**Problema affrontato:** Il modello originale utilizzava una probabilità di dropout uniforme o dipendente solo dalla posizione spaziale, ma nei dati reali il dropout è fortemente correlato al livello di espressione.

**Soluzione implementata:**
- Modello di dropout dipendente dal livello di espressione del gene
- Funzione logistica per modellare la relazione espressione-dropout
- Combinazione con il modello di dropout spaziale

```r
if (dropout_params$expression_dependent_dropout) {
  # Normalizza tra 0 e 1
  norm_expr <- scale01(expression_data[, g])
  
  # Calcola la probabilità di dropout con una funzione logistica
  # Formula: p(dropout) = 1 / (1 + exp((expr - midpoint) * steepness))
  dropout_prob_expr <- 1 / (1 + exp((norm_expr - dropout_params$dropout_curve_midpoint) * 
                                   dropout_params$dropout_curve_steepness))
  
  # Combina con il dropout spaziale base (media pesata)
  dropout_prob <- 0.7 * dropout_prob_expr + 0.3 * base_dropout
  
  # Applica dropout
  zero_idx <- rbinom(N, 1, dropout_prob) == 1
  expression_data[zero_idx, g] <- 0
}
```

**Fondamento teorico:** Numerosi studi (Kharchenko et al. 2014, Pierson & Yau 2015) hanno dimostrato che nei dati di trascrittomica, la probabilità di osservare un falso zero (dropout) è inversamente correlata al livello reale di espressione del gene. Questa relazione segue tipicamente una funzione sigmoide o logistica:

p(dropout) = 1 / (1 + exp((x - x₀) * k))

dove:
- x è il livello di espressione normalizzato
- x₀ è il punto di flesso (dove p(dropout) = 0.5)
- k è il parametro di pendenza che controlla la ripidità della transizione

Questa relazione si basa sul fatto che geni altamente espressi hanno maggiori probabilità di essere catturati anche in presenza di inefficienze tecniche.

### 3.3 Dispersione Variabile per Tipo Cellulare e Posizione

**Problema affrontato:** La variabilità biologica (dispersione) è spesso eterogenea in base sia alla posizione nel tessuto che al tipo cellulare, un fenomeno non adeguatamente catturato dal modello originale.

**Soluzione implementata:**
- Dispersione variabile basata sulla distanza dal confine tra cluster
- Effetto del tipo cellulare sulla dispersione 
- Parametrizzazione avanzata della dispersione nei dati simulati

```r
# Dispersione basata sulla posizione (distanza dal confine)
if (spatial_params$gradient_regions && exists("boundary_dist", where = cell_df)) {
  # Più vicino al confine = più variabilità (parametro dispersione più basso)
  dispersion_param <- dropout_params$dispersion_range[2] + 
    cell_df$boundary_dist * (dropout_params$dispersion_range[1] - dropout_params$dispersion_range[2])
} else {
  # Metodo originale basato sulla distanza media
  dispersion_param <- rescale(mean_dist, to = dropout_params$dispersion_range)
}

# Aggiungi effetto del tipo cellulare sulla dispersione
# Alcuni tipi cellulari sono intrinsecamente più variabili di altri
set.seed(random_seed + 4)
# Crea effetti casuali per ogni tipo cellulare
type_dispersion_effects <- runif(k_cell_types, 
                               min = 1 - dropout_params$cell_type_dispersion_effect,
                               max = 1 + dropout_params$cell_type_dispersion_effect)

# Applica effetto moltiplicativo per tipo cellulare
for (k in 1:k_cell_types) {
  dispersion_param[cluster_labels == k] <- dispersion_param[cluster_labels == k] * type_dispersion_effects[k]
}
```

**Fondamento teorico:** Nella distribuzione Binomiale Negativa, il parametro di dispersione (size) è inversamente correlato alla variabilità extra-Poissoniana. Un valore basso di questo parametro comporta alta variabilità.

L'implementazione si basa sull'osservazione che:

1. Le aree di transizione tra regioni tissutali mostrano maggiore eterogeneità trascrizionale
2. Diversi tipi cellulari hanno intrinsecamente diversi livelli di variabilità trascrizionale (es. cellule staminali vs. cellule differenziate)
3. La variabilità dell'espressione genica è un fenomeno multifattoriale che combina effetti spaziali, di tipo cellulare e di stato cellulare

Studi come Edsgärd et al. (2018) e Tirosh et al. (2016) documentano queste differenze di variabilità sia a livello spaziale che a livello di tipo cellulare, evidenziando l'importanza di una modellazione realistica della dispersione.

### 3.4 Implementazione di Moduli Genici Co-espressi

**Problema affrontato:** Il modello originale trattava ogni gene come indipendente, mentre nei dati reali esistono gruppi di geni con espressione coordinata che formano moduli o programmi trascrizionali.

**Soluzione implementata:**
- Creazione di moduli di geni co-espressi
- Implementazione di rumore correlato all'interno di ciascun modulo
- Parametrizzazione del livello di correlazione tra geni dello stesso modulo

```r
# Crea moduli di geni co-espressi
if (cell_specific_params$use_gene_modules) {
  # Calcola il numero di geni per modulo
  genes_per_module <- ceiling(n_genes / cell_specific_params$n_gene_modules)
  
  # Assegna geni ai moduli
  gene_modules <- list()
  for (m in 1:cell_specific_params$n_gene_modules) {
    start_idx <- (m-1) * genes_per_module + 1
    end_idx <- min(m * genes_per_module, n_genes)
    gene_modules[[m]] <- start_idx:end_idx
  }
  
  # Crea rumore correlato per ogni modulo
  module_noise <- matrix(0, nrow = N, ncol = n_genes)
  
  # Genera rumore base per ogni modulo
  base_module_noise <- matrix(rnorm(cell_specific_params$n_gene_modules * N), 
                             nrow = N, ncol = cell_specific_params$n_gene_modules)
  
  # Applica il rumore del modulo a ciascun gene appartenente al modulo
  for (m in 1:length(gene_modules)) {
    module_genes <- gene_modules[[m]]
    # Assegna lo stesso rumore base a tutti i geni del modulo, scalato per la correlazione
    module_noise[, module_genes] <- base_module_noise[, m] * cell_specific_params$module_correlation
  }
  
  # Aggiunta al modello di espressione
  if (cell_specific_params$use_gene_modules && !is.null(module_noise)) {
    # Aggiungi il rumore correlato del modulo genico a cui appartiene questo gene
    mu_vals <- mu_vals + module_noise[, g]
  }
}
```

**Fondamento teorico:** Nei sistemi biologici reali, i geni non funzionano in modo indipendente ma sono organizzati in programmi trascrizionali coordinati che regolano processi cellulari specifici. Questi programmi creano pattern di co-espressione dove gruppi di geni mostrano correlazione positiva.

Studi come quelli di Stuart et al. (2003), Zhang & Horvath (2005) e Buettner et al. (2017) hanno dimostrato che:

1. I moduli di co-espressione sono una caratteristica fondamentale dei dati di trascrittomica
2. Questi moduli spesso corrispondono a pathway biologici o programmi funzionali
3. La correlazione all'interno dei moduli varia in base al tipo di tessuto e allo stato cellulare

La simulazione di questi pattern di co-espressione è essenziale per creare dati sintetici che riflettano realisticamente la complessità dei sistemi biologici reali.

## 4. Visualizzazione e Validazione

### 4.1 Visualizzazione Avanzata

**Problema affrontato:** La visualizzazione originale non mostrava adeguatamente la struttura a griglia e non forniva informazioni su caratteristiche chiave come dimensione della libreria e autocorrelazione.

**Soluzione implementata:**
- Visualizzazione specifica per griglia utilizzando geom_tile
- Plot della distribuzione spaziale della dimensione della libreria
- Visualizzazione dell'indice di Moran I per i geni selezionati

```r
# Visualizzazione griglia
if (grid_mode) {
  p <- ggplot(cell_df, aes(x = x, y = y, fill = intensity_cluster)) +
    geom_tile(width = grid_resolution, height = grid_resolution) +
    scale_y_reverse() +
    coord_fixed()
}

# Plot della dimensione libreria
p_lib <- ggplot(data.frame(x = cell_df$x, y = cell_df$y, lib_size = library_size), 
               aes(x = x, y = y, fill = lib_size)) +
  geom_tile(width = grid_resolution, height = grid_resolution) +
  scale_fill_viridis_c(option = "plasma")

# Plot dell'autocorrelazione
p_auto <- ggplot(data.frame(gene = result$spatial_autocorrelation$genes_tested,
                           moran_i = result$spatial_autocorrelation$moran_i), 
                aes(x = gene, y = moran_i)) +
  geom_bar(stat = "identity", fill = "steelblue")
```

**Fondamento teorico:** La visualizzazione efficace dei dati simulati è cruciale per:
1. Validare visivamente la somiglianza con dati reali
2. Comprendere la distribuzione spaziale di caratteristiche chiave
3. Verificare il comportamento dei modelli di autocorrelazione implementati

Le visualizzazioni implementate seguono le prassi comuni nell'analisi di dati Visium, dove la dimensione della libreria e l'autocorrelazione spaziale sono diagnostici fondamentali per valutare la qualità del dataset.

### 4.2 Output Arricchito per Validazione

**Problema affrontato:** L'output originale mancava di informazioni cruciali per validare la simulazione rispetto ai dati reali.

**Soluzione implementata:**
- Aggiunta della dimensione della libreria all'output
- Calcolo e inclusione di metriche di autocorrelazione spaziale
- Informazioni aggiuntive sui parametri utilizzati

```r
result <- list(
  coordinates       = cell_df[, c("x", "y")],
  intensity_cluster = cell_df$intensity_cluster,
  expression        = expression_data,
  threshold_used    = threshold_value,
  library_size      = library_size,  # Dimensione libreria per spot
  parameters        = list(...),     # Tutti i parametri
  spatial_autocorrelation = list(    # Metriche di autocorrelazione
    moran_i = moran_i,
    genes_tested = 1:moran_genes
  )
)
```

**Fondamento teorico:** La validazione quantitativa è essenziale per garantire che i dati simulati rappresentino accuratamente le caratteristiche statistiche dei dati reali. Le metriche di autocorrelazione spaziale, in particolare l'indice di Moran I, sono ampiamente utilizzate nella letteratura per quantificare i pattern spaziali di espressione genica (Svensson et al. 2018, Hu et al. 2021).

## 5. Miglioramenti Recenti e Direzioni Future

### 5.1 Miglioramenti Implementati

Le recenti modifiche hanno trasformato significativamente il framework di simulazione, rendendolo molto più aderente alla realtà biologica dei dati Visium HD:

1. **Miglioramento dei gradienti tra regioni:** 
   - Implementazione di una funzione di transizione non lineare (esponente personalizzabile)
   - Considerazione di molteplici tipi cellulari vicini anziché solo coppie
   - Transizioni più graduali che riflettono l'eterogeneità dei confini tissutali

2. **Moduli genici co-espressi:**
   - Simulazione di programmi trascrizionali coordinati
   - Correlazione parametrizzabile tra geni dello stesso modulo
   - Riproduzione di pathway biologici realistici

3. **Effetti specifici del tipo cellulare:**
   - Variazione della dimensione della libreria per tipo cellulare
   - Dispersione genica specifica per tipo cellulare
   - Riproduzione dell'eterogeneità trascrizionale osservata in diversi tipi di cellule

4. **Autocorrelazione spaziale migliorata:**
   - Parametrizzazione avanzata dei modelli GRF e CAR
   - Pattern spaziali più marcati con valori di Moran I più realistici
   - Maggiore variabilità di struttura spaziale tra geni

Questi miglioramenti hanno portato a simulazioni con caratteristiche molto più simili ai dati Visium HD reali, con pattern spaziali complessi, correlazioni tra geni biologicamente plausibili, e variabilità realistica sia a livello tecnico che biologico.

### 5.2 Direzioni Future

Nonostante i significativi progressi, diverse aree rimangono aperte per ulteriori sviluppi:

1. **Dinamiche temporali:** Implementare la simulazione di processi dinamici come sviluppo, differenziamento o risposte a stimoli, che potrebbero essere modellati come serie temporali di dati spaziali.

2. **Interazioni ligando-recettore:** Modellare esplicitamente le interazioni cellula-cellula attraverso coppie ligando-recettore, creando pattern di espressione spazialmente dipendenti basati su meccanismi di comunicazione intercellulare.

3. **Inferenza di parametri da dati reali:** Sviluppare metodi per stimare automaticamente i parametri di simulazione a partire da dati Visium HD reali, per creare simulazioni "gemelle digitali" di dataset specifici.

4. **Modelli gerarchici multi-scala:** Implementare pattern spaziali a diverse scale di risoluzione (organismi, tessuti, regioni, singole cellule) che meglio riflettono l'organizzazione gerarchica dei sistemi biologici.

5. **Integrazione multi-omica:** Estendere il framework per simulare simultaneamente più livelli di dati omici (trascrittoma, proteoma, epigenoma) con le loro reciproche dipendenze.

Il framework continuerà ad evolversi per rispondere alle crescenti esigenze di benchmark avanzati per metodi di analisi di dati di trascrittomica spaziale ad alta risoluzione.

## Bibliografia

- Edsgärd, D., Johnsson, P., & Sandberg, R. (2018). Identification of spatial expression trends in single-cell gene expression data. Nature Methods, 15(5), 339-342.
- Hu, J., Li, X., Coleman, K., et al. (2021). SpaGCN: Integrating gene expression, spatial location and histology to identify spatial domains and spatially variable genes by graph convolutional network. Nature Methods, 18, 1342–1351.
- Kharchenko, P. V., Silberstein, L., & Scadden, D. T. (2014). Bayesian approach to single-cell differential expression analysis. Nature Methods, 11(7), 740-742.
- Lytal, N., Ran, D., & An, L. (2020). Normalization methods on single-cell RNA-seq data: An empirical survey. Frontiers in Genetics, 11, 41.
- Pierson, E., & Yau, C. (2015). ZIFA: Dimensionality reduction for zero-inflated single-cell gene expression analysis. Genome Biology, 16(1), 1-10.
- Svensson, V., Teichmann, S. A., & Stegle, O. (2018). SpatialDE: identification of spatially variable genes. Nature Methods, 15(5), 343-346.