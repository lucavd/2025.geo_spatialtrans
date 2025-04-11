#' Carica e prepara un'immagine per la simulazione
#'
#' @param image_path Path dell'immagine da caricare
#' @param threshold_value Soglia per il thresholding dell'immagine (default: 0.7)
#' @return Lista contenente l'immagine preparata e i dati derivati
#' @importFrom imager load.image rm.alpha grayscale squeeze as.cimg
#' @importFrom dplyr select filter
#' @export
prepare_image <- function(image_path, threshold_value = 0.7) {
  if (is.null(image_path)) {
    stop("È necessario specificare il percorso dell'immagine (image_path)")
  }
  
  # Caricamento immagine
  img <- load.image(image_path)
  
  # Gestione canali alpha e conversione a scala di grigi
  if (spectrum(img) == 4) {
    img <- rm.alpha(img)
  }
  if (spectrum(img) == 3) {
    img <- grayscale(img)
  }
  
  img <- squeeze(img)
  
  # Converto in array e poi in cimg (x,y,cc,t)
  img_array <- as.array(img)
  w <- nrow(img_array)
  h <- ncol(img_array)
  img_cimg <- as.cimg(img_array, dims = c(w, h, 1, 1))
  
  # Converto in data frame (x, y, value)
  img_df <- as.data.frame(img_cimg) %>%
    dplyr::select(x, y, value)
  
  # Applico soglia
  img_df_thresh <- img_df %>%
    filter(value < threshold_value)
  
  # Ritorno i dati generati
  return(list(
    img = img,
    img_cimg = img_cimg,
    img_array = img_array,
    img_df = img_df,
    img_df_thresh = img_df_thresh,
    width = w,
    height = h
  ))
}