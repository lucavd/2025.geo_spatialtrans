# Hard Clustering Challenge

## Strategie per rendere più difficile il riconoscimento dei cluster nei dati simulati

Il codice originale `image_simulations_2_claude.R` creava dati con cluster troppo ben definiti, rendendo facile per algoritmi come Seurat identificarli. Le seguenti modifiche sono state implementate per creare una versione più challenging:

### 1. Riduzione del differenziale di espressione dei geni marker
- Ridotto il numero di geni marker specifici per cellula da 10 a 5
- Diminuito il differenziale di espressione da +2 a +0.8
- Aggiunto overlapping parziale (0.4) tra cluster adiacenti
```r
# Riduzione del numero e dell'intensità dei marker
start_idx <- (k - 1) * 5 + 1  # Ridotto da 10 a 5 geni marker per tipo
end_idx   <- min(k * 5, n_genes)
if (start_idx <= end_idx) {
  mu[start_idx:end_idx] <- mu[start_idx:end_idx] + 0.8  # Ridotto da +2 a +0.8
}

# Aggiungi espressione parziale nei cluster adiacenti (overlapping)
if (k > 1) {
  prev_markers <- ((k-2) * 5 + 1):min((k-1) * 5, n_genes)
  if (length(prev_markers) > 0) {
    mu[prev_markers] <- mu[prev_markers] + 0.4  # Espressione parziale del cluster precedente
  }
}
```

### 2. Aumento del rumore spaziale
- Incrementato l'intensità del rumore spaziale da 0.5 a 1.5
- Ridotto il range di correlazione spaziale da 30 a 15
- Aggiunto noise random per ulteriore variabilità cellula-specifica
```r
# Crea un modello gstat per il GP con range più corto (più irregolarità spaziali)
gp_sim <- gstat(formula = z ~ 1, locations = ~x+y, dummy = TRUE,
              beta = 0, model = vgm(psill=1.5, range=15, model="Exp"), nmax=10)

# Aggiungi correlazione spaziale se richiesta (più intensa)
if (use_spatial_correlation) {
  mu_vals <- mu_vals + 1.5 * gp_noise  # Aumentato da 0.5 a 1.5
  
  # Aggiungi noise casuale addizionale per confondere i pattern
  random_noise <- rnorm(length(mu_vals), 0, 0.4)
  mu_vals <- mu_vals + random_noise
}
```

### 3. Aumento del dropout e variabilità
- Aumentato il dropout in tutto il tessuto (da 0.1-0.5 a 0.3-0.7)
- Ridotta la differenza di dispersione tra centro e bordi dei cluster
```r
# Dispersione molto elevata ovunque (con poca differenza centro/bordi)
dispersion_param <- rescale(mean_dist, to = c(1.5, 0.8))  # Valori più bassi = più variabilità

# Dropout probability elevata ovunque
dropout_prob <- rescale(mean_dist, to = c(0.3, 0.7))  # Maggior dropout ovunque
```

### 4. Creazione di cellule ibride ai confini tra cluster
- Identificazione delle coppie di cellule vicine appartenenti a cluster diversi
- Creazione di una matrice di ibridazione che rappresenta quanto ogni cellula ha caratteristiche di altri cluster
- Applicazione dell'effetto ibrido all'espressione genica
```r
# Crea cellule ibride ai confini tra cluster
hybrid_pairs <- list()
for (i in 1:N) {
  # Trova cellule vicine (tra le 20 più vicine)
  neighbors <- order(dist_mat[i,])[2:20]
  diff_cluster_neighbors <- neighbors[cluster_labels[neighbors] != cluster_labels[i]]
  
  # Se ci sono vicini di cluster diversi, aggiungi alla lista
  if (length(diff_cluster_neighbors) > 0) {
    hybrid_pairs[[length(hybrid_pairs) + 1]] <- c(i, diff_cluster_neighbors[1])
  }
}

# Applica effetto di ibridazione tra cluster (cellule ai confini)
for (k in 1:k_cell_types) {
  # Trova cellule che hanno componente del cluster k
  hybrid_cells <- which(hybrid_matrix[, k] > 0)
  if (length(hybrid_cells) > 0) {
    # Applica l'effetto ibrido usando il profilo del cluster k
    hybrid_effect <- mean_expression_list[[k]][g] * hybrid_matrix[hybrid_cells, k]
    base_expr[hybrid_cells] <- base_expr[hybrid_cells] * (1 - hybrid_matrix[hybrid_cells, k]) + hybrid_effect
  }
}
```

### 5. Aggiunta di variabilità cellula-specifica
- Ogni cellula riceve un effetto casuale indipendente dal suo cluster
```r
# Aggiungi variabilità cellula-specifica indipendente dal cluster
cell_specific_effect <- rnorm(N, 0, 0.3)
```

## Effetti previsti
Queste modifiche dovrebbero rendere il riconoscimento dei cluster un compito molto più difficile per metodi come Seurat, perché:

1. I confini tra cluster diventano sfumati con cellule ibride e overlapping di espressione
2. Il rumore spaziale e la variabilità cellula-specifica mascherano i pattern di espressione
3. L'aumento del dropout crea ulteriore incertezza nel segnale
4. La riduzione del numero e dell'intensità dei geni marker rende più difficile distinguere i tipi cellulari

Queste sfide sono più rappresentative dei dati di trascrittomica spaziale reali, dove i confini tra tipi cellulari sono spesso ambigui e i gradienti di espressione creano transizioni graduali.