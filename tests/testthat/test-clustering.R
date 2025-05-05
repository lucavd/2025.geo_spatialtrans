test_that("cluster_image funziona correttamente con dati semplici", {
  # Creiamo un dataset di test semplice
  test_df <- data.frame(
    x = rep(1:10, each = 10),
    y = rep(1:10, times = 10),
    value = c(rep(0.1, 50), rep(0.5, 50))
  )
  
  # Applica la funzione
  result <- cluster_image(test_df, k_cell_types = 2, random_seed = 123)
  
  # Verifica il risultato
  expect_s3_class(result, "data.frame")
  expect_true("intensity_cluster" %in% names(result))
  expect_s3_class(result$intensity_cluster, "factor")
  expect_equal(nlevels(result$intensity_cluster), 2)
  
  # Verifica che i cluster siano coerenti con i valori
  # (assumendo che i valori bassi siano in un cluster, quelli alti in un altro)
  cluster_low <- as.integer(result$intensity_cluster[result$value == 0.1][1])
  cluster_high <- as.integer(result$intensity_cluster[result$value == 0.5][1])
  
  expect_true(cluster_low != cluster_high)
  expect_true(all(result$intensity_cluster[result$value == 0.1] == result$intensity_cluster[result$value == 0.1][1]))
  expect_true(all(result$intensity_cluster[result$value == 0.5] == result$intensity_cluster[result$value == 0.5][1]))
})

test_that("spatial_kmeans è il metodo predefinito e retrocompatibilità funziona", {
  # Dataset di test
  test_df <- data.frame(
    x = rep(1:10, each = 10),
    y = rep(1:10, times = 10),
    value = c(rep(0.1, 50), rep(0.5, 50))
  )
  
  # Testa che spatial_kmeans sia il metodo predefinito
  result_default <- cluster_image(test_df, k_cell_types = 2, random_seed = 123)
  
  # Testa retrocompatibilità con kmeans++
  result_kmeans <- cluster_image(test_df, k_cell_types = 2, random_seed = 123, 
                                clustering_method = "kmeans++")
  
  # Verifica che entrambi producano risultati validi
  expect_s3_class(result_default, "data.frame")
  expect_s3_class(result_kmeans, "data.frame")
  expect_true("intensity_cluster" %in% names(result_default))
  expect_true("intensity_cluster" %in% names(result_kmeans))
  
  # I risultati non devono necessariamente essere uguali
  # ma entrambi dovrebbero avere cluster significativi
  expect_equal(nlevels(result_default$intensity_cluster), 2)
  expect_equal(nlevels(result_kmeans$intensity_cluster), 2)
})

test_that("spatial_kmeans genera risultati diversi con diversi pesi spaziali", {
  # Dataset di test con pattern spaziali chiari
  test_df <- data.frame(
    x = rep(1:20, each = 20),
    y = rep(1:20, times = 20),
    value = rep(c(0.1, 0.5), each = 200)
  )
  
  # Applica spatial_kmeans con diversi pesi
  result_w0 <- cluster_image(test_df, k_cell_types = 2, random_seed = 123,
                           clustering_method = "spatial_kmeans", spatial_weight = 0)
  result_w2 <- cluster_image(test_df, k_cell_types = 2, random_seed = 123,
                           clustering_method = "spatial_kmeans", spatial_weight = 2)
  
  # Con peso 0, il risultato dovrebbe essere basato solo sui valori (simile a kmeans++)
  # Con peso alto, il risultato dovrebbe essere più influenzato dalla posizione spaziale
  
  # Calcola la frammentazione spaziale dei cluster
  # (più cluster adiacenti con etichette diverse = più frammentazione)
  count_transitions <- function(df) {
    transitions <- 0
    for (i in 1:(nrow(df)-1)) {
      if (df$x[i] == df$x[i+1] && df$y[i]+1 == df$y[i+1]) {
        if (df$intensity_cluster[i] != df$intensity_cluster[i+1]) {
          transitions <- transitions + 1
        }
      }
    }
    return(transitions)
  }
  
  # Con peso spaziale alto, ci aspettiamo meno transizioni (cluster più coesi)
  transitions_w0 <- count_transitions(result_w0)
  transitions_w2 <- count_transitions(result_w2)
  
  expect_true(transitions_w2 <= transitions_w0)
})

test_that("stima automatica del numero di cluster funziona", {
  # Dataset con 3 cluster chiari
  test_df <- data.frame(
    x = rep(1:10, each = 10),
    y = rep(1:10, times = 10),
    value = c(rep(0.1, 30), rep(0.5, 30), rep(0.9, 40))
  )
  
  # Test con stima automatica (silhouette)
  result_auto_sil <- cluster_image(test_df, k_cell_types = 2, random_seed = 123,
                                 estimate_k = TRUE, k_estimation_method = "silhouette")
  
  # Test con stima automatica (elbow)
  result_auto_elbow <- cluster_image(test_df, k_cell_types = 2, random_seed = 123,
                                   estimate_k = TRUE, k_estimation_method = "elbow")
  
  # Verifica che entrambi i metodi producano risultati validi
  expect_s3_class(result_auto_sil, "data.frame")
  expect_s3_class(result_auto_elbow, "data.frame")
  
  # I metodi potrebbero stimare 2 o 3 cluster, entrambi accettabili per questo dataset
  expect_true(nlevels(result_auto_sil$intensity_cluster) >= 2)
  expect_true(nlevels(result_auto_elbow$intensity_cluster) >= 2)
})

# Test con SLIC solo se supercells è disponibile
if (requireNamespace("supercells", quietly = TRUE)) {
  test_that("SLIC clustering funziona correttamente", {
    # Dataset di test con pattern spaziali
    test_df <- data.frame(
      x = rep(1:15, each = 15),
      y = rep(1:15, times = 15),
      value = c(rep(0.1, 112), rep(0.5, 113))
    )
    
    # Applica SLIC clustering
    result_slic <- cluster_image(test_df, k_cell_types = 2, random_seed = 123,
                               clustering_method = "slic")
    
    # Verifica il risultato
    expect_s3_class(result_slic, "data.frame")
    expect_true("intensity_cluster" %in% names(result_slic))
    expect_equal(nlevels(result_slic$intensity_cluster), 2)
  })
}