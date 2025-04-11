test_that("scale01 normalizza correttamente", {
  # Test di input normali
  expect_equal(scale01(c(1, 2, 3, 4, 5)), c(0, 0.25, 0.5, 0.75, 1))
  
  # Test con valori negativi
  expect_equal(scale01(c(-10, 0, 10)), c(0, 0.5, 1))
  
  # Test con valori tutti uguali
  expect_equal(scale01(rep(5, 3)), c(0.5, 0.5, 0.5))
  
  # Test con valori NA
  expect_error(scale01(c(1, 2, NA)))
})

test_that("scale01_vec funziona come scale01 ma è vettorizzato", {
  # Test con dati normali
  x <- c(1, 2, 3, 4, 5)
  expect_equal(scale01_vec(x), scale01(x))
  
  # Test con valori tutti uguali
  x <- rep(7, 5)
  expect_equal(scale01_vec(x), scale01(x))
})