#' Crea una griglia di campionamento per la simulazione
#'
#' Genera una griglia di campionamento spaziale regolare o un campionamento
#' casuale in base ai parametri specificati.
#'
#' @param img_df_thresh Dataframe con i pixel filtrati con intensity_cluster
#' @param img_array Array dell'immagine originale
#' @param img_width Larghezza immagine in pixel
#' @param img_height Altezza immagine in pixel
#' @param grid_mode Se TRUE utilizza una griglia regolare, altrimenti sampling casuale
#' @param n_cells Numero di celle da campionare (per grid_mode=FALSE)
#' @param k_cell_types Numero di tipi cellulari
#' @param grid_resolution Risoluzione della griglia in μm
#' @param grid_spacing Spazio tra bin della griglia
#' @param use_fixed_grid Utilizzo di una griglia di dimensioni fisse
#' @param fixed_grid_width_mm Larghezza della griglia fissa in mm
#' @param fixed_grid_height_mm Altezza della griglia fissa in mm
#' @param pixel_size_um Dimensione di ogni pixel in μm
#' @param threshold_value Soglia per il thresholding
#' @param random_seed Seed per riproducibilità
#' @return Dataframe con le celle/bin campionate
#' @importFrom dplyr filter select mutate rowwise ungroup bind_rows
#' @importFrom tidyr expand_grid
#' @export
create_sampling_grid <- function(img_df_thresh, img_array, img_width, img_height,
                                grid_mode = TRUE, n_cells = 20000, k_cell_types = 5,
                                grid_resolution = 2, grid_spacing = 0,
                                use_fixed_grid = FALSE, fixed_grid_width_mm = 6.5,
                                fixed_grid_height_mm = 6.5, pixel_size_um = 1,
                                threshold_value = 0.7, random_seed = 123) {
  
  # Se grid_mode è TRUE, crea una griglia regolare
  if (grid_mode) {
    # Calcola le dimensioni dell'immagine in μm
    img_width_um <- img_width * pixel_size_um
    img_height_um <- img_height * pixel_size_um
    
    # Determina le dimensioni della griglia e le coordinate
    if (use_fixed_grid) {
      # Usa una griglia fissa con dimensioni standard (6.5mm x 6.5mm)
      # Converti da mm a μm
      grid_width_um <- fixed_grid_width_mm * 1000  # 6.5mm = 6500μm
      grid_height_um <- fixed_grid_height_mm * 1000  # 6.5mm = 6500μm
      
      # Calcola il numero di bin necessari
      n_bins_x <- ceiling(grid_width_um / grid_resolution)
      n_bins_y <- ceiling(grid_height_um / grid_resolution)
      
      # Crea le coordinate della griglia fissa in μm
      x_coords <- seq(0, grid_width_um - grid_resolution, by = grid_resolution + grid_spacing)
      y_coords <- seq(0, grid_height_um - grid_resolution, by = grid_resolution + grid_spacing)
    } else {
      # Usa una griglia che si adatta all'immagine
      grid_width_um <- img_width_um
      grid_height_um <- img_height_um
      
      # Calcola il numero di bin necessari
      n_bins_x <- ceiling(img_width_um / grid_resolution)
      n_bins_y <- ceiling(img_height_um / grid_resolution)
      
      # Crea le coordinate della griglia adattata all'immagine in μm
      x_coords <- seq(0, img_width_um - grid_resolution, by = grid_resolution + grid_spacing)
      y_coords <- seq(0, img_height_um - grid_resolution, by = grid_resolution + grid_spacing)
    }
    
    # Crea un dataframe con tutte le coordinate possibili
    grid_points <- expand.grid(x = x_coords, y = y_coords)
    
    # Per la griglia fissa, assicuriamoci che i punti corrispondano all'immagine
    if (use_fixed_grid) {
      # Calcola l'offset per centrare l'immagine nella griglia, se necessario
      if (grid_width_um > img_width_um) {
        offset_x <- (grid_width_um - img_width_um) / 2
      } else {
        offset_x <- 0
      }
      
      if (grid_height_um > img_height_um) {
        offset_y <- (grid_height_um - img_height_um) / 2
      } else {
        offset_y <- 0
      }
      
      # Assegna a ciascun punto della griglia il valore dell'immagine, 
      # se il punto è all'interno dell'immagine
      grid_df <- grid_points %>%
        rowwise() %>%
        mutate(
          # Calcola le coordinate relative all'immagine, considerando l'offset
          rel_x = x - offset_x,
          rel_y = y - offset_y,
          
          # Controlla se il punto è all'interno dell'immagine
          is_in_image = (rel_x >= 0 && rel_x < img_width_um && rel_y >= 0 && rel_y < img_height_um),
          
          # Se il punto è fuori dall'immagine, assegna un valore superiore alla soglia
          # altrimenti prendi il valore dall'immagine
          value = if (is_in_image) {
            # Converti da coordinate μm a indici di pixel
            img_x = min(max(round(rel_x / pixel_size_um), 1), img_width)
            img_y = min(max(round(rel_y / pixel_size_um), 1), img_height)
            img_array[img_x, img_y]
          } else {
            1.0  # Valore superiore alla soglia, sarà filtrato
          }
        ) %>%
        filter(value < threshold_value) %>%  # Applica la soglia
        ungroup()
    } else {
      # Comportamento originale per la griglia adattata all'immagine
      grid_df <- grid_points %>%
        rowwise() %>%
        mutate(
          # Converti da μm a indici di pixel nell'immagine originale
          img_x = min(max(round(x / pixel_size_um), 1), img_width),
          img_y = min(max(round(y / pixel_size_um), 1), img_height),
          value = img_array[img_x, img_y]
        ) %>%
        filter(value < threshold_value) %>%  # Applica la stessa soglia
        ungroup()
    }
    
    # Assegna cluster usando i centroidi k-means dell'immagine originale
    # Trova i centroidi dei cluster
    cluster_centroids <- sapply(levels(img_df_thresh$intensity_cluster), function(cl) {
      mean(img_df_thresh$value[img_df_thresh$intensity_cluster == cl])
    })
    
    # Assegna ogni punto griglia al cluster più vicino
    grid_df <- grid_df %>%
      mutate(intensity_cluster = factor(apply(outer(value, cluster_centroids, 
                                               FUN = function(x, y) abs(x - y)),
                                        1, which.min)))
    
    # Questa è la nostra "cell_df" finale
    cell_df <- grid_df
    
  } else {
    # Modalità campionamento casuale
    # Calcola le frequenze dei cluster
    clust_counts <- table(img_df_thresh$intensity_cluster)
    clust_freq   <- clust_counts / sum(clust_counts)
    cells_per_cluster <- round(clust_freq * n_cells)
    
    # Campiona le cellule proporzionalmente da ciascun cluster
    cell_list <- vector("list", length(levels(img_df_thresh$intensity_cluster)))
    set.seed(random_seed)
    
    for (i in seq_along(levels(img_df_thresh$intensity_cluster))) {
      clust_name <- levels(img_df_thresh$intensity_cluster)[i]
      n_sub      <- cells_per_cluster[clust_name]
      df_sub     <- img_df_thresh %>% filter(intensity_cluster == clust_name)
      idx_sub    <- sample(seq_len(nrow(df_sub)), min(n_sub, nrow(df_sub)), replace = FALSE)
      cell_list[[i]] <- df_sub[idx_sub, ]
    }
    
    cell_df <- bind_rows(cell_list)
  }
  
  return(cell_df)
}