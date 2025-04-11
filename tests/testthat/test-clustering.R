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