#' Generate Synthetic Tissue Image
#'
#' Quickly creates a grayscale synthetic tissue image with varying complexity
#' to be used as input for spatial transcriptomics simulation. The function
#' outputs the image matrix and a data.frame of pixel intensities, and can
#' optionally save a PNG file on disk.
#'
#' @param width_px Integer. Image width in pixels.
#' @param height_px Integer. Image height in pixels.
#' @param complexity Integer in 1:4. 1 = Gaussian blobs, 2 = Voronoi patches,
#'   3 = mix + fractal noise, 4 = complexity 3 + biological structures.
#' @param seed Integer. Random seed for reproducibility.
#' @param output_path Character or NA. If not NA, the PNG image is written to
#'   this path.
#' @return List with `img_matrix` (height×width numeric) and `img_df`
#'   (data.frame with x, y, intensity).
#' @examples
#' img <- generate_synthetic_tissue(512, 512, complexity = 2, seed = 42,
#'                                  output_path = "synthetic.png")
#' @export
generate_synthetic_tissue <- function(width_px = 6800, height_px = 6500,
                                      complexity = 2, seed = 123,
                                      output_path = NA) {
  start_time <- Sys.time()
  stopifnot(complexity %in% 1:4)
  cat(sprintf("\n[synthetic_tissue] Generating synthetic tissue (%dx%d, complexity=%d)\n",
              width_px, height_px, complexity))
  set.seed(seed)

  # Base blank image
  img_mat <- matrix(255, nrow = height_px, ncol = width_px) # initialize white

  if (complexity == 1) {
    ## Gaussian blobs ---------------------------------------------------------
    n_blobs <- 8
    cat("  - Placing", n_blobs, "Gaussian blobs...\n")
    for (i in seq_len(n_blobs)) {
      cx <- sample(width_px, 1)
      cy <- sample(height_px, 1)
      sigma <- sample(seq(150, 500, by = 50), 1)
      # Draw blob via outer product of Gaussians
      x_vec <- seq_len(width_px)
      y_vec <- seq_len(height_px)
      gx <- exp(-((x_vec - cx)^2) / (2 * sigma^2))
      gy <- exp(-((y_vec - cy)^2) / (2 * sigma^2))
      blob <- 200 * (gy %o% gx) # intensity scale 0-200
      img_mat <- pmin(img_mat, 255 - blob) # darker blobs on white background
    }
  } else if (complexity == 2) {
    ## Voronoi patches --------------------------------------------------------
    cat("  - Creating Voronoi patches...\n")
    n_centers <- 30
    centers <- cbind(x = runif(n_centers, 1, width_px),
                     y = runif(n_centers, 1, height_px))

    # Precompute coordinate grids (reuse memory)
    xs <- matrix(rep(seq_len(width_px), each = height_px), nrow = height_px)
    ys <- matrix(rep(seq_len(height_px), width_px), nrow = height_px)

    best_d  <- matrix(Inf, nrow = height_px, ncol = width_px)
    patch_id <- matrix(NA_integer_, nrow = height_px, ncol = width_px)

    for (k in seq_len(n_centers)) {
      d <- abs(xs - centers[k, "x"]) + abs(ys - centers[k, "y"]) # Manhattan distance as proxy
      better <- d < best_d
      patch_id[better] <- k
      best_d[better]  <- d[better]
    }

    intensities <- sample(seq(50, 200, by = 5), n_centers, replace = TRUE)
    img_mat <- matrix(intensities[patch_id], nrow = height_px)
  } else if (complexity == 3) {
    ## Mix: blobs + Voronoi + fractal noise ----------------------------------
    cat("  - Generating mixed pattern (blobs + patches + noise)...\n")
    tmp <- generate_synthetic_tissue(width_px, height_px, complexity = 1,
                                     seed = seed * 7 + 123, output_path = NA)
    img_mat <- tmp$img_matrix
    tmp2 <- generate_synthetic_tissue(width_px, height_px, complexity = 2,
                                      seed = seed * 13 + 456, output_path = NA)
    img_mat <- pmin(img_mat, tmp2$img_matrix)
    # Add simple fractal-like noise via Perlin approximation (random clouds)
    set.seed(seed * 17 + 789)  # Set different seed for noise generation
    n_noise <- 6
    cat("  - Adding", n_noise, "noise layers...\n")
    for (i in seq_len(n_noise)) {
      scale <- 2^(i + 1)
      noise <- matrix(runif(ceiling(height_px/scale) * ceiling(width_px/scale)),
                      nrow = ceiling(height_px/scale))
      noise <- noise[rep(seq_len(nrow(noise)), each = scale, length.out = height_px),
                     rep(seq_len(ncol(noise)), each = scale, length.out = width_px)]
      img_mat <- pmin(img_mat, 255 - 30 * noise)
    }
  } else if (complexity == 4) {
    ## Mix + biological structures (vessels/fibers) ---------------------------
    cat("  - Generating mixed pattern with biological structures...\n")
    
    # Start with complexity 3 base (with different seed for variation)
    tmp <- generate_synthetic_tissue(width_px, height_px, complexity = 3,
                                     seed = seed * 19 + 321, output_path = NA)
    img_mat <- tmp$img_matrix
    
    # Add few realistic biological structures (vessels/fibers)
    set.seed(seed * 23 + 999)  # Different seed for biological structures
    n_structures <- 3  # Few realistic biological structures
    cat("  - Adding", n_structures, "biological structures (vessels/fibers)...\n")
    xs <- matrix(rep(seq_len(width_px), each = height_px), nrow = height_px)
    ys <- matrix(rep(seq_len(height_px), width_px), nrow = height_px)
    
    for (i in seq_len(n_structures)) {
      # Define 3 control points for gentle Bezier curves (biological vessels)
      p0 <- c(runif(1, width_px * 0.1, width_px * 0.9), 
              runif(1, height_px * 0.1, height_px * 0.9))
      p1 <- c(runif(1, width_px * 0.2, width_px * 0.8), 
              runif(1, height_px * 0.2, height_px * 0.8))
      p2 <- c(runif(1, width_px * 0.1, width_px * 0.9), 
              runif(1, height_px * 0.1, height_px * 0.9))
      
      # Parameter t from 0 to 1 (fewer points for efficiency)
      t <- seq(0, 1, length.out = 200)
      
      # Bezier curve formula
      x_curve <- (1 - t)^2 * p0[1] + 2 * (1 - t) * t * p1[1] + t^2 * p2[1]
      y_curve <- (1 - t)^2 * p0[2] + 2 * (1 - t) * t * p1[2] + t^2 * p2[2]
      
      # Biological vessel/fiber parameters
      structure_thickness <- sample(15:35, 1)  # Realistic vessel thickness
      intensity_val <- runif(1, 0.3, 0.6)     # Moderate darkness
      
      for (j in seq_along(x_curve)) {
        dist_sq <- (xs - x_curve[j])^2 + (ys - y_curve[j])^2
        mask <- dist_sq < structure_thickness^2
        # Gaussian falloff for natural vessel appearance
        falloff <- exp(-dist_sq[mask] / (2 * structure_thickness^2))
        img_mat[mask] <- pmin(img_mat[mask], 255 * (1 - intensity_val * falloff))
      }
    }
  }

  # Clip and convert to 0-1 for PNG
  img_norm <- (img_mat - min(img_mat)) / (max(img_mat) - min(img_mat))

  if (!is.na(output_path)) {
    cat("  - Writing PNG to", output_path, "...\n")
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    png::writePNG(img_norm, target = output_path)
  }

  # Build data.frame (sample every pixel or stride 2 to save RAM)
  stride <- max(1, floor(sqrt((width_px * height_px)/1e6))) # keep <=1M rows
  idx_x <- seq_len(width_px)[(seq_len(width_px) - 1) %% stride == 0]
  idx_y <- seq_len(height_px)[(seq_len(height_px) - 1) %% stride == 0]
  df <- expand.grid(x = idx_x, y = idx_y)
  df$intensity <- img_mat[df$y + (df$x - 1) * height_px]

  end_time <- Sys.time()
  cat(sprintf("[synthetic_tissue] Done in %.1fs. Output df rows: %d\n",
              as.numeric(difftime(end_time, start_time, units = "secs")), nrow(df)))

  list(img_matrix = img_mat, img_df = df)
}
