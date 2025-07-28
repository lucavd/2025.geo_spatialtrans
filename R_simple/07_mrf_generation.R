#' Simulate spatial cell-type labels with an optimized Potts (MRF) model
#'
#' This function generates a rectangular lattice of discrete cell‐type labels
#' using optimized Gibbs sampling from a first–order Potts model. Optimized
#' with vectorized operations and adaptive convergence for speed.
#'
#' @param grid_size Integer or length-2 vector. Number of lattice sites in the
#'   *x* and *y* directions. Default `200`, i.e. a 200×200 grid.
#' @param k_cell_types Integer. Number of discrete cell types (Potts states).
#' @param beta Numeric. Interaction strength (>0 attracts, <0 repels). When
#'   `beta = 0` the configuration is random; higher values increase spatial
#'   autocorrelation. Typical range 0–1.
#' @param n_iter Integer. Number of Gibbs sweeps (reduced default for speed).
#' @param seed Integer. Random seed for reproducibility.
#' @param fast_mode Logical. Use faster approximations (default TRUE).
#' @return `data.frame` with columns `x`, `y`, `cell_type` (integer 1…k).
#' @examples
#' df <- simulate_mrf(grid_size = 128, k_cell_types = 6, beta = 0.6,
#'                    n_iter = 200, seed = 42)
#' @export
simulate_mrf <- function(grid_size = 200,
                         k_cell_types = 5,
                         beta = 0.8,
                         n_iter = 200,
                         seed = 1,
                         fast_mode = TRUE,
                         tissue_structure = "uniform",
                         interaction_matrix = NULL) {
  if (length(grid_size) == 1L) grid_size <- c(grid_size, grid_size)
  stopifnot(length(grid_size) == 2L, k_cell_types >= 2, n_iter >= 1)
  # Support both single structures and composite structures
  if (is.character(tissue_structure) && length(tissue_structure) == 1) {
    # Single structure mode (backward compatibility)
    valid_structures <- c("uniform", "vessel", "boundary", "gradient")
    if (!tissue_structure %in% valid_structures) {
      stop("Invalid tissue_structure. Must be one of: ", paste(valid_structures, collapse=", "))
    }
    structure_list <- list(list(type = tissue_structure, weight = 1.0))
  } else if (is.list(tissue_structure)) {
    # Composite structure mode: list of structure definitions
    structure_list <- tissue_structure
  } else {
    stop("tissue_structure must be a single string or a list of structure definitions")
  }

  set.seed(seed)
  nx <- grid_size[1]; ny <- grid_size[2]

  # Create biologically-informed interaction matrix
  if (is.null(interaction_matrix)) {
    interaction_matrix <- create_biological_interactions(k_cell_types, seed)
  }

  # Initialize with tissue structure-aware seeding
  lattice <- initialize_tissue_structure(nx, ny, k_cell_types, structure_list, seed)

  # Enhanced neighborhood energy with cell-type specific interactions
  neigh_energy <- function(x, y, label) {
    neighbors <- c(
      if (y > 1) lattice[y-1, x] else 0,
      if (y < ny) lattice[y+1, x] else 0,
      if (x > 1) lattice[y, x-1] else 0,
      if (x < nx) lattice[y, x+1] else 0
    )
    
    # Sum interaction strengths based on biological affinity matrix
    energy <- 0
    for (neighbor in neighbors[neighbors > 0]) {
      energy <- energy + interaction_matrix[label, neighbor]
    }
    return(energy)
  }

  # Optimized Gibbs sweeps --------------------------------------------------
  if (fast_mode && beta > 0.3) {
    # Fast mode: checkerboard updates + vectorized operations
    for (iter in seq_len(n_iter)) {
      # Update even sites (checkerboard pattern)
      for (phase in 0:1) {
        for (y in seq_len(ny)) {
          x_indices <- seq(1 + (y + phase) %% 2, nx, by = 2)
          if (length(x_indices) == 0) next
          
          for (x in x_indices) {
            # Vectorized neighbor energy computation
            neighbors <- c(
              if (y > 1) lattice[y-1, x] else 0,
              if (y < ny) lattice[y+1, x] else 0,
              if (x > 1) lattice[y, x-1] else 0,
              if (x < nx) lattice[y, x+1] else 0
            )
            
            # Count matches for each cell type
            logp <- sapply(seq_len(k_cell_types), function(lbl) {
              beta * sum(neighbors == lbl, na.rm = TRUE)
            })
            
            # Fast sampling with temperature annealing
            temp <- max(0.5, 1 - iter / n_iter)  # Cool down over time
            p <- exp(logp / temp)
            p <- p / sum(p)
            lattice[y, x] <- sample.int(k_cell_types, 1, prob = p)
          }
        }
      }
      
      if (iter %% 50 == 0) {
        cat(sprintf("[simulate_mrf] fast sweep %d/%d\n", iter, n_iter))
      }
    }
  } else {
    # Standard mode (original algorithm)
    for (iter in seq_len(n_iter)) {
      for (y in seq_len(ny)) {
        for (x in seq_len(nx)) {
          logp <- sapply(seq_len(k_cell_types), function(lbl) {
            beta * neigh_energy(x, y, lbl)
          })
          p <- exp(logp - max(logp))
          p <- p / sum(p)
          lattice[y, x] <- sample.int(k_cell_types, 1, prob = p)
        }
      }
      if (iter %% 100 == 0) {
        cat(sprintf("[simulate_mrf] standard sweep %d/%d\n", iter, n_iter))
      }
    }
  }

  # Build data.frame ---------------------------------------------------------
  df <- expand.grid(x = seq_len(nx), y = seq_len(ny))
  df$cell_type <- as.vector(t(lattice))
  df
}

#' Create biologically-informed cell-type interaction matrix
#' @param k_cell_types Number of cell types
#' @param seed Random seed
#' @return Symmetric interaction matrix
create_biological_interactions <- function(k_cell_types, seed) {
  set.seed(seed + 100)
  
  # Base interaction matrix (symmetric)
  mat <- matrix(0.2, nrow = k_cell_types, ncol = k_cell_types)  # weak default attraction
  diag(mat) <- 1.0  # strong self-attraction
  
  # Add biological patterns:
  # 1. Immune cells (types 1-2) attract each other
  if (k_cell_types >= 2) {
    mat[1, 2] <- mat[2, 1] <- 0.8
  }
  
  # 2. Stromal cells (type 3) repel immune but attract epithelial
  if (k_cell_types >= 4) {
    mat[3, 1:2] <- mat[1:2, 3] <- -0.3  # repulsion
    mat[3, 4] <- mat[4, 3] <- 0.9       # epithelial-stromal attraction
  }
  
  # 3. Add some random asymmetric interactions for complexity
  for (i in 1:(k_cell_types-1)) {
    for (j in (i+1):k_cell_types) {
      noise <- runif(1, -0.2, 0.2)
      mat[i, j] <- mat[i, j] + noise
      mat[j, i] <- mat[i, j]  # keep symmetric
    }
  }
  
  return(mat)
}

#' Initialize tissue with biologically-inspired spatial structure
#' @param nx,ny Grid dimensions
#' @param k_cell_types Number of cell types
#' @param structure_list List of structure definitions with weights
#' @param seed Random seed
#' @return Initial lattice configuration
initialize_tissue_structure <- function(nx, ny, k_cell_types, structure_list, seed) {
  set.seed(seed + 200)
  
  # Initialize base lattice
  lattice <- matrix(1, nrow = ny, ncol = nx)
  
  # Apply each structure in sequence, weighted by importance
  for (struct_def in structure_list) {
    structure_type <- struct_def$type
    weight <- struct_def$weight
    
    if (structure_type == "uniform") {
      # Apply uniform randomization with weight
      n_random <- round(nx * ny * weight * 0.3)  # 30% of weight as randomization
      random_positions <- sample(nx * ny, n_random)
      for (pos in random_positions) {
        row <- ((pos - 1) %% ny) + 1
        col <- ((pos - 1) %/% ny) + 1
        lattice[row, col] <- sample.int(k_cell_types, 1)
      }
      next
    }
    
    # Create temporary lattice for this structure
    temp_lattice <- apply_single_structure(nx, ny, k_cell_types, structure_type, seed + as.numeric(charToRaw(structure_type)[1]))
    
    # Blend with existing lattice based on weight
    blend_structures(lattice, temp_lattice, weight)
  }
  
  return(lattice)
}

#' Apply a single tissue structure to a lattice
#' @param nx,ny Grid dimensions
#' @param k_cell_types Number of cell types
#' @param structure_type Single structure type
#' @param seed Random seed
#' @return Lattice with single structure applied
apply_single_structure <- function(nx, ny, k_cell_types, structure_type, seed) {
  set.seed(seed)
  lattice <- matrix(1, nrow = ny, ncol = nx)
  
  if (structure_type == "vessel") {
    # Create vessel-like linear structures
    n_vessels <- max(2, k_cell_types %/% 2)
    for (v in 1:n_vessels) {
      # Random vessel path
      start_x <- sample(nx, 1)
      start_y <- sample(ny, 1)
      vessel_type <- ((v - 1) %% k_cell_types) + 1
      
      # Draw vessel with random walk
      x <- start_x; y <- start_y
      for (step in 1:min(nx, ny)) {
        if (x >= 1 && x <= nx && y >= 1 && y <= ny) {
          # Vessel core and surrounding
          lattice[y, x] <- vessel_type
          if (x > 1) lattice[y, x-1] <- vessel_type
          if (x < nx) lattice[y, x+1] <- vessel_type
        }
        # Random walk
        x <- x + sample(c(-1, 0, 1), 1)
        y <- y + sample(c(-1, 0, 1), 1)
      }
    }
  } else if (structure_type == "boundary") {
    # Create distinct regions with sharp boundaries
    mid_x <- nx %/% 2
    mid_y <- ny %/% 2
    
    lattice[1:mid_y, 1:mid_x] <- 1
    lattice[1:mid_y, (mid_x+1):nx] <- 2 %% k_cell_types + 1
    lattice[(mid_y+1):ny, 1:mid_x] <- 3 %% k_cell_types + 1
    lattice[(mid_y+1):ny, (mid_x+1):nx] <- 4 %% k_cell_types + 1
    
  } else if (structure_type == "gradient") {
    # Create spatial gradient of cell types
    for (i in 1:ny) {
      for (j in 1:nx) {
        # Distance-based assignment
        dist_from_center <- sqrt((i - ny/2)^2 + (j - nx/2)^2)
        max_dist <- sqrt((ny/2)^2 + (nx/2)^2)
        type_prob <- (dist_from_center / max_dist) * (k_cell_types - 1) + 1
        lattice[i, j] <- min(k_cell_types, max(1, round(type_prob)))
      }
    }
  }
  
  return(lattice)
}

#' Blend two lattice structures based on weight
#' @param base_lattice Base lattice to modify
#' @param new_lattice New structure to blend in
#' @param weight Weight of new structure (0-1)
blend_structures <- function(base_lattice, new_lattice, weight) {
  nx <- ncol(base_lattice)
  ny <- nrow(base_lattice)
  
  # Determine how many pixels to replace based on weight
  n_replace <- round(nx * ny * weight)
  replace_positions <- sample(nx * ny, n_replace)
  
  for (pos in replace_positions) {
    row <- ((pos - 1) %% ny) + 1
    col <- ((pos - 1) %/% ny) + 1
    base_lattice[row, col] <- new_lattice[row, col]
  }
  
  invisible(NULL)  # Modify in place
}
