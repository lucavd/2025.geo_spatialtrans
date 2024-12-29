#!/usr/bin/env Rscript

# Required libraries with error checking
for(pkg in c("torch", "igraph", "future", "future.apply")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

#' Create spatial graph from cell coordinates
create_spatial_graph <- function(coordinates, expression, distance_threshold = 50) {
  coord_matrix <- as.matrix(coordinates)
  dist_matrix <- as.matrix(dist(coord_matrix))
  adjacency_matrix <- (dist_matrix < distance_threshold) * 1

  graph <- graph_from_adjacency_matrix(
    adjacency_matrix,
    mode = "undirected",
    weighted = TRUE
  )

  node_features <- as.matrix(expression)

  list(
    graph = graph,
    features = node_features,
    adjacency = adjacency_matrix
  )
}

#' Basic GNN model
GNNModel <- nn_module(
  initialize = function(in_features, hidden_dim = 64, out_features = 2) {
    self$conv1 <- nn_linear(in_features, hidden_dim)
    self$conv2 <- nn_linear(hidden_dim, hidden_dim)
    self$out <- nn_linear(hidden_dim, out_features)
    self$dropout <- nn_dropout(0.5)
  },

  forward = function(x, adj) {
    # First layer
    x <- torch_matmul(adj, x)
    x <- self$conv1(x)
    x <- nnf_relu(x)
    x <- self$dropout(x)

    # Second layer
    x <- torch_matmul(adj, x)
    x <- self$conv2(x)
    x <- nnf_relu(x)
    x <- self$dropout(x)

    # Output layer
    x <- self$out(x)
    x
  }
)

#' Process single file
process_file <- function(data_path) {
  cat("Processing", basename(data_path), "\n")

  # Load data
  data <- readRDS(data_path)

  # Create graph structure
  graph_data <- create_spatial_graph(
    coordinates = data$coordinates,
    expression = data$expression
  )

  # Prepare tensors
  adj_tensor <- torch_tensor(graph_data$adjacency)
  feature_tensor <- torch_tensor(graph_data$features)
  labels <- torch_tensor(as.integer(data$true_hotspots))

  # Initialize model
  model <- GNNModel(
    in_features = ncol(graph_data$features),
    hidden_dim = 64,
    out_features = length(unique(data$true_hotspots))
  )

  # Training setup
  optimizer <- optim_adam(model$parameters, lr = 0.01)

  # Training loop with early stopping
  patience <- 10  # Number of epochs to wait for improvement
  min_delta <- 1e-4  # Minimum change in loss to qualify as an improvement
  best_loss <- Inf
  patience_counter <- 0
  loss_history <- numeric()
  max_epochs <- 500  # Maximum number of epochs to prevent infinite loops

  for(epoch in 1:max_epochs) {
    model$train()
    optimizer$zero_grad()

    output <- model(feature_tensor, adj_tensor)
    loss <- nn_cross_entropy_loss()(output, labels)
    current_loss <- as.numeric(loss)

    loss$backward()
    optimizer$step()

    # Store loss history
    loss_history <- c(loss_history, current_loss)

    # Early stopping check
    if (current_loss < best_loss - min_delta) {
      best_loss <- current_loss
      patience_counter <- 0
    } else {
      patience_counter <- patience_counter + 1
    }

    # Print progress
    if(epoch %% 5 == 0) {
      cat("Epoch", epoch, "Loss:", round(current_loss, 6),
          "Best:", round(best_loss, 6), "\n")
    }

    # Check if we should stop
    if (patience_counter >= patience) {
      cat("\nEarly stopping at epoch", epoch,
          "- Loss hasn't improved for", patience, "epochs\n")
      break
    }

    # Additional stability check - stop if loss is stable enough
    if (epoch > 20) {
      recent_losses <- tail(loss_history, 10)
      loss_std <- sd(recent_losses)
      if (loss_std < min_delta) {
        cat("\nStopping at epoch", epoch,
            "- Loss has stabilized (std:", round(loss_std, 6), ")\n")
        break
      }
    }
  }

  # Plot final loss curve
  plot(1:length(loss_history), loss_history,
       type = "l",
       xlab = "Epoch",
       ylab = "Loss",
       main = paste("Training Loss -", basename(data_path)))

  # Generate predictions
  model$eval()
  predictions <- torch_argmax(model(feature_tensor, adj_tensor), dim = 2)
  predictions <- as.integer(as.array(predictions))

  # Save results
  result <- list(
    hotspots = predictions,
    model = model
  )

  output_path <- file.path("results", paste0("gnn_results_", basename(data_path)))
  saveRDS(result, file = output_path)

  cat("Saved results to:", output_path, "\n")
  result
}

# Main execution function
run_analysis <- function() {
  # Create results directory if needed
  if (!dir.exists("results")) {
    dir.create("results")
    cat("Created results directory\n")
  }

  # Get data files
  data_paths <- list.files(
    path = "data",
    pattern = "simulated_.*_correlation.rds",
    full.names = TRUE
  )

  if(length(data_paths) == 0) {
    stop("No data files found in data/ directory")
  }

  cat("Found", length(data_paths), "files to process\n")

  # Process each file
  for(data_path in data_paths) {
    tryCatch({
      process_file(data_path)
    }, error = function(e) {
      cat("Error processing", basename(data_path), ":", e$message, "\n")
    })
  }
}

# Execute analysis
run_analysis()
