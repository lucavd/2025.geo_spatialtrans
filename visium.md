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
  # Crea un GP per l'effetto spaziale
  sp_df_lib <- cell_df
  coordinates(sp_df_lib) <- ~ x + y
  
  lib_gp <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
                 beta = 0, model = vgm(psill = 1.0, 
                                      range = spatial_params$spatial_range * 1.5, 
                                      model = "Exp"), 
                 nmax = 10)
  
  lib_noise <- predict(lib_gp, newdata = sp_df_lib, nsim = 1)$sim1
  
  # Applica alla dimensione libreria (effetto moltiplicativo)
  lib_noise <- scale(lib_noise)
  lib_effect <- library_size_params$spatial_effect_on_library * lib_noise
  library_size <- library_size * exp(lib_effect)
}

# Applica l'effetto della dimensione della libreria ai conteggi
scaled_counts <- raw_counts * (library_size / mean(library_size))
expression_data[, g] <- round(scaled_counts)
```

**Fondamento teorico:** La distribuzione della dimensione della libreria nei dati di trascrittomica segue tipicamente una distribuzione log-normale (Lytal et al. 2020). Nei dati Visium reali, la dimensione della libreria mostra anche pattern spaziali dovuti a:
1. Variazione dell'efficienza di cattura RNA in diverse regioni tissutali
2. Differenze nella densità cellulare tra regioni
3. Effetti tecnici durante la preparazione delle librerie

La modellazione realistica di questa variazione è fondamentale per simulare correttamente la variabilità tecnica nei dati.

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

### 3.3 Dispersione Variabile Spazialmente

**Problema affrontato:** La variabilità biologica (dispersione) è spesso maggiore ai confini tra regioni tissutali, un fenomeno non adeguatamente catturato dal modello originale.

**Soluzione implementata:**
- Dispersione variabile basata sulla distanza dal confine tra cluster
- Parametrizzazione della dispersione nei dati simulati

```r
if (spatial_params$gradient_regions && exists("boundary_dist", where = cell_df)) {
  # Più vicino al confine = più variabilità (parametro dispersione più basso)
  dispersion_param <- dropout_params$dispersion_range[2] + 
    cell_df$boundary_dist * (dropout_params$dispersion_range[1] - dropout_params$dispersion_range[2])
} else {
  # Metodo originale basato sulla distanza media
  dispersion_param <- rescale(mean_dist, to = dropout_params$dispersion_range)
}
```

**Fondamento teorico:** Nella distribuzione Binomiale Negativa, il parametro di dispersione (size) è inversamente correlato alla variabilità extra-Poissoniana. Un valore basso di questo parametro comporta alta variabilità.

L'implementazione si basa sull'osservazione che le aree di transizione tra regioni tissutali mostrano maggiore eterogeneità trascrizionale dovuta a:

1. Cellule in stati intermedi o transitori
2. Misture di popolazioni cellulari diverse
3. Risposte a segnali contrastanti dalle regioni adiacenti

Questo fenomeno è documentato in studi come Edsgärd et al. (2018) che mostrano come confini tissutali e gradienti funzionali siano caratterizzati da maggiore variabilità di espressione.

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

## 5. Conclusione e Sviluppi Futuri

Le modifiche implementate hanno trasformato il framework di simulazione in uno strumento molto più accurato per replicare le caratteristiche dei dati Visium HD. In particolare:

1. La modalità griglia riproduce fedelmente la struttura spaziale di Visium HD
2. I gradienti e i modelli di correlazione spaziale catturano pattern biologici realistici
3. La modellazione della dimensione della libreria e del dropout riflette le sfide tecniche dei dati reali

### Potenziali Sviluppi Futuri

1. **Composizione cellulare mista:** Implementare la simulazione esplicita di composizioni cellulari miste in ogni spot, particolarmente rilevante per bin a 2μm che possono comunque contenere frammenti di più cellule.

2. **Correlazione tra geni:** Aggiungere la possibilità di simulare strutture di correlazione complesse tra geni, riflettendo pathway e programmi trascrizionali coordinati.

3. **Struttura di covarianza stimata da dati reali:** Permettere l'importazione di matrici di covarianza stimate da dati Visium HD reali per guidare la simulazione.

4. **Integrazione con dati di imaging multipli:** Estendere il framework per simulare dati multi-modali che integrano trascrittomica spaziale con imaging di proteine o altre modalità.

## Bibliografia

- Edsgärd, D., Johnsson, P., & Sandberg, R. (2018). Identification of spatial expression trends in single-cell gene expression data. Nature Methods, 15(5), 339-342.
- Hu, J., Li, X., Coleman, K., et al. (2021). SpaGCN: Integrating gene expression, spatial location and histology to identify spatial domains and spatially variable genes by graph convolutional network. Nature Methods, 18, 1342–1351.
- Kharchenko, P. V., Silberstein, L., & Scadden, D. T. (2014). Bayesian approach to single-cell differential expression analysis. Nature Methods, 11(7), 740-742.
- Lytal, N., Ran, D., & An, L. (2020). Normalization methods on single-cell RNA-seq data: An empirical survey. Frontiers in Genetics, 11, 41.
- Pierson, E., & Yau, C. (2015). ZIFA: Dimensionality reduction for zero-inflated single-cell gene expression analysis. Genome Biology, 16(1), 1-10.
- Svensson, V., Teichmann, S. A., & Stegle, O. (2018). SpatialDE: identification of spatially variable genes. Nature Methods, 15(5), 343-346.