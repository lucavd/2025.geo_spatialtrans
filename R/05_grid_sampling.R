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
#' @param grid_resolution Risoluzione della griglia in um
#' @param grid_spacing Spazio tra bin della griglia
#' @param use_fixed_grid Utilizzo di una griglia di dimensioni fisse
#' @param fixed_grid_width_mm Larghezza della griglia fissa in mm
#' @param fixed_grid_height_mm Altezza della griglia fissa in mm
#' @param pixel_size_um Dimensione di ogni pixel in um
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
    # Calcola le dimensioni dell'immagine in um
    img_width_um <- img_width * pixel_size_um
    img_height_um <- img_height * pixel_size_um
    
    # Determina le dimensioni della griglia e le coordinate
    if (use_fixed_grid) {
      # Usa una griglia fissa con dimensioni standard (6.5mm x 6.5mm)
      # Converti da mm a um
      grid_width_um <- fixed_grid_width_mm * 1000  # 6.5mm = 6500um
      grid_height_um <- fixed_grid_height_mm * 1000  # 6.5mm = 6500um
      
      # Calcola il numero di bin necessari
      n_bins_x <- ceiling(grid_width_um / grid_resolution)
      n_bins_y <- ceiling(grid_height_um / grid_resolution)
      
      # Crea le coordinate della griglia fissa in um
      # Aggiungi controlli per evitare parametri invalidi
      if (grid_width_um <= grid_resolution) {
        grid_width_um <- grid_resolution * 2  # Forza almeno 2 punti
      }
      if (grid_height_um <= grid_resolution) {
        grid_height_um <- grid_resolution * 2  # Forza almeno 2 punti
      }
      
      # Assicurati che il passo sia positivo
      step <- max(0.1, grid_resolution + grid_spacing)
      
      x_coords <- seq(0, grid_width_um - grid_resolution, by = step)
      y_coords <- seq(0, grid_height_um - grid_resolution, by = step)
    } else {
      # Usa una griglia che si adatta all'immagine
      grid_width_um <- img_width_um
      grid_height_um <- img_height_um
      
      # Calcola il numero di bin necessari
      n_bins_x <- ceiling(img_width_um / grid_resolution)
      n_bins_y <- ceiling(img_height_um / grid_resolution)
      
      # Crea le coordinate della griglia adattata all'immagine in μm
      # Aggiungi controlli per evitare parametri invalidi
      if (img_width_um <= grid_resolution) {
        img_width_um <- grid_resolution * 2  # Forza almeno 2 punti
      }
      if (img_height_um <= grid_resolution) {
        img_height_um <- grid_resolution * 2  # Forza almeno 2 punti
      }
      
      # Assicurati che il passo sia positivo
      step <- max(0.1, grid_resolution + grid_spacing)
      
      x_coords <- seq(0, img_width_um - grid_resolution, by = step)
      y_coords <- seq(0, img_height_um - grid_resolution, by = step)
    }
    
    # Crea un dataframe con tutte le coordinate possibili
    grid_points <- expand.grid(x = x_coords, y = y_coords)
    
    # Per la griglia fissa, assicuriamoci che i punti corrispondano all'immagine
    if (use_fixed_grid) {
      # Calcola l'offset per centrare l'immagine nella griglia
      # (se l'immagine è più grande della griglia, centriamo; altrimenti offset a 0)
      offset_x <- max((img_width_um  - grid_width_um ) / 2, 0)
      offset_y <- max((img_height_um - grid_height_um) / 2, 0)
      
      # Convert grid coordinates to image coordinates and assign image values
      rel_x <- grid_points$x + offset_x
      rel_y <- grid_points$y + offset_y
      # Mask for points inside the image
      in_img <- rel_x >= 0 & rel_x < img_width_um & rel_y >= 0 & rel_y < img_height_um
      # Initialize values (outside -> 1.0 to be filtered)
      value <- rep(1.0, length(rel_x))
      if (any(in_img)) {
        ix <- floor(rel_x[in_img] / pixel_size_um) + 1
        iy <- floor(rel_y[in_img] / pixel_size_um) + 1
        ix <- pmin(pmax(ix, 1), img_width)
        iy <- pmin(pmax(iy, 1), img_height)
        value[in_img] <- img_array[cbind(ix, iy)]
      }
      # Filter by threshold
      keep <- value < threshold_value
      grid_df <- data.frame(
        x = grid_points$x[keep],
        y = grid_points$y[keep],
        value = value[keep]
      )
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
    
      # Assign cluster using k-means centroids of the original image
      cluster_levels    <- levels(img_df_thresh$intensity_cluster)
      cluster_centroids <- sapply(cluster_levels, function(cl) {
        mean(img_df_thresh$value[img_df_thresh$intensity_cluster == cl])
      })
      # Find nearest centroid per point (min abs difference)
      n_pts <- nrow(grid_df)
      best_idx  <- integer(n_pts)
      best_dist <- abs(grid_df$value - cluster_centroids[1])
      best_idx[] <- 1
      for (j in seq_along(cluster_centroids)[-1]) {
        d_j <- abs(grid_df$value - cluster_centroids[j])
        mask <- d_j < best_dist
        if (any(mask)) {
          best_dist[mask] <- d_j[mask]
          best_idx[mask]  <- j
        }
      }
      grid_df$intensity_cluster <- factor(cluster_levels[best_idx], levels = cluster_levels)
      # Final cell_df
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