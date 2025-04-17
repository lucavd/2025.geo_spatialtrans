test_that("prepare_image gestisce correttamente gli errori", {
  # Test con image_path NULL
  expect_error(
    prepare_image(image_path = NULL),
    "È necessario specificare un percorso di immagine valido"
  )
  
  # Test con path non esistente - testiamo solo che generi un errore, senza verificare il messaggio specifico
  expect_error(
    prepare_image(image_path = "path/non/esistente.png")
  )
})

# Nota: per test più completi sarebbe necessario un file immagine di esempio
# che potrebbe essere incluso nel pacchetto come fixture per i test

test_that("prepare_image funziona con immagine sintetica", {
  skip_if_not_installed("png")
  # Creazione di un'immagine sintetica 2x3
  mat <- matrix(c(0.1, 0.6, 0.9,
                  0.2, 0.7, 0.3),
                nrow = 2, byrow = TRUE)
  tmp <- tempfile(fileext = ".png")
  png::writePNG(mat, tmp)
  res <- prepare_image(tmp, threshold_value = 0.5)
  expect_type(res, "list")
  expect_equal(res$width, 3L)
  expect_equal(res$height, 2L)
  expect_equal(dim(res$img_array), c(3L, 2L))
  df <- res$img_df_thresh
  expect_true(all(c("x", "y", "value") %in% colnames(df)))
  expect_equal(nrow(df), 3L)
  expect_equal(sort(df$value), c(0.1, 0.2, 0.3), tolerance = 1e-2)
})