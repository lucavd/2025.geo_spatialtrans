#!/usr/bin/env Rscript
# Test rapido del pipeline completo con clustering basato su espressione
# Versione ridotta per validare l'integrazione

# Modifica temporanea del file principale per test rapido
main_file <- "R/run_full_size_optimized.R"

# Leggi il file
lines <- readLines(main_file)

# Trova e modifica i parametri per test rapido
lines <- gsub("n_genes\\s*=\\s*20000", "n_genes = 1000", lines)
lines <- gsub("grid_resolution\\s*=\\s*30", "grid_resolution = 80", lines)
lines <- gsub("fixed_grid_width_mm\\s*=\\s*6\\.5", "fixed_grid_width_mm = 3.0", lines)
lines <- gsub("fixed_grid_height_mm\\s*=\\s*6\\.5", "fixed_grid_height_mm = 3.0", lines)
lines <- gsub("width_px\\s*=\\s*6800", "width_px = 1500", lines)
lines <- gsub("height_px\\s*=\\s*6500", "height_px = 1500", lines)
lines <- gsub("chunk_size <- 2000", "chunk_size <- 500", lines)

# Salva il file modificato temporaneamente
temp_file <- "R/test_temp_full_size.R"
writeLines(lines, temp_file)

cat("=== TEST PIPELINE COMPLETO CON CLUSTERING ESPRESSIONE ===\n")
cat("File temporaneo creato:", temp_file, "\n")
cat("Parametri ridotti per test rapido\n")
cat("Eseguendo simulazione...\n\n")

# Esegui il test
system(paste("Rscript", temp_file))

# Pulisci
file.remove(temp_file)
cat("\n=== TEST COMPLETATO ===\n")