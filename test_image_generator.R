library(png)

# Create a simple test image
img <- matrix(0, nrow = 50, ncol = 50)

# Add some patterns
for (i in 1:50) {
  for (j in 1:50) {
    if (i < 25 && j < 25) {
      img[i, j] <- 0.9
    } else if (i >= 25 && j >= 25) {
      img[i, j] <- 0.7
    } else if (i < 25 && j >= 25) {
      img[i, j] <- 0.5
    } else {
      img[i, j] <- 0.3
    }
  }
}

# Create directory if it doesn't exist
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

# Save the test image
writePNG(img, target = "inst/extdata/test_image.png")

cat("Test image saved to inst/extdata/test_image.png\n")