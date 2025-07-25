#' Carica e prepara un'immagine per la simulazione
#'
#' @param image_path Path dell'immagine da caricare
#' @param threshold_value Soglia per il thresholding dell'immagine (default: 0.7)
#' @return Lista contenente l'immagine preparata e i dati derivati
#' @importFrom png readPNG
#' @importFrom dplyr select filter
#' @export
prepare_image <- function(image_path, threshold_value = 0.7) {
  if (is.null(image_path) || !file.exists(image_path)) {
    stop("È necessario specificare un percorso di immagine valido")
  }
  if (!requireNamespace("png", quietly = TRUE)) {
    stop("Package 'png' is required; please install it to read images")
  }
  # Read image (returns [height x width x channels] or [height x width])
  img_raw <- png::readPNG(image_path)
  # Convert to grayscale by averaging RGB channels if needed
  if (length(dim(img_raw)) == 3) {
    # assume at least 3 channels (RGB); ignore alpha if present
    gray <- (img_raw[,,1] + img_raw[,,2] + img_raw[,,3]) / 3
  } else {
    gray <- img_raw
  }
  # Transpose so that array dims match original semantics (width x height)
  img_array <- t(gray)
  w <- nrow(img_array)
  h <- ncol(img_array)
  # Build dataframe of pixel values
  grid <- expand.grid(x = seq_len(w), y = seq_len(h))
  grid$value <- as.vector(img_array)
  # Filter by threshold
  img_df_thresh <- grid[grid$value < threshold_value, , drop = FALSE]
  # Return prepared image data
  return(list(
    img_array     = img_array,
    img_df_thresh = img_df_thresh,
    width         = w,
    height        = h
  ))
}