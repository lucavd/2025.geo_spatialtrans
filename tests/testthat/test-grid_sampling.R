test_that("create_sampling_grid crea una griglia regolare", {
  # Creiamo un dataset di test semplice e un'immagine di test
  test_df <- data.frame(
    x = rep(1:10, each = 10),
    y = rep(1:10, times = 10),
    value = c(rep(0.1, 50), rep(0.5, 50)),
    intensity_cluster = factor(c(rep(1, 50), rep(2, 50)))
  )
  
  # Immagine di test
  test_img_array <- matrix(0.2, nrow = 10, ncol = 10)
  
  # Test con modalità griglia standard
  result <- create_sampling_grid(
    img_df_thresh = test_df,
    img_array = test_img_array,
    img_width = 10,
    img_height = 10,
    grid_mode = TRUE,
    grid_resolution = 2,
    grid_spacing = 0,
    use_fixed_grid = FALSE,
    random_seed = 123
  )
  
  # Verifica del risultato
  expect_s3_class(result, "data.frame")
  expect_true("intensity_cluster" %in% names(result))
  expect_equal(unique(diff(sort(unique(result$x)))), 2) # Spaziatura x corretta
  expect_equal(unique(diff(sort(unique(result$y)))), 2) # Spaziatura y corretta
})

test_that("create_sampling_grid funziona in modalità campionamento casuale", {
  # Preparazione dati di test
  test_df <- data.frame(
    x = rep(1:10, each = 10),
    y = rep(1:10, times = 10),
    value = c(rep(0.1, 50), rep(0.5, 50)),
    intensity_cluster = factor(c(rep(1, 50), rep(2, 50)))
  )
  
  # Immagine di test
  test_img_array <- matrix(0.2, nrow = 10, ncol = 10)
  
  # Test con modalità campionamento
  result <- create_sampling_grid(
    img_df_thresh = test_df,
    img_array = test_img_array,
    img_width = 10,
    img_height = 10,
    grid_mode = FALSE,
    n_cells = 20,
    k_cell_types = 2,
    random_seed = 123
  )
  
  # Verifica del risultato
  expect_s3_class(result, "data.frame")
  expect_true("intensity_cluster" %in% names(result))
  expect_equal(nrow(result), 20)
  
  # Verifica che ci siano cellule di entrambi i cluster
  expect_length(unique(result$intensity_cluster), 2)
})