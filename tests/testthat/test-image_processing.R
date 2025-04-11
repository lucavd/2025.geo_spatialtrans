test_that("prepare_image gestisce correttamente gli errori", {
  # Test con image_path NULL
  expect_error(
    prepare_image(image_path = NULL),
    "È necessario specificare il percorso dell'immagine"
  )
  
  # Test con path non esistente - testiamo solo che generi un errore, senza verificare il messaggio specifico
  expect_error(
    prepare_image(image_path = "path/non/esistente.png")
  )
})

# Nota: per test più completi sarebbe necessario un file immagine di esempio
# che potrebbe essere incluso nel pacchetto come fixture per i test