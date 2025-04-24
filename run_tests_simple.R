# Script semplificato per testare le versioni semplificate dei moduli

# Carica tutte le librerie necessarie
library(tictoc)
library(dplyr)
library(tidyr)
library(ggplot2)
library(future)
library(future.apply)
library(imager)
library(MASS)
library(cluster)
library(ClusterR)
library(gstat)
library(sp)
library(scales)
library(spatstat)
library(spdep)
library(fields)
library(testthat)

# Configura il futuro per la parallelizzazione
future::plan(future::multisession, workers = 2)
options(future.globals.maxSize = 100 * 1024^2) # 100 GB

# Carica le funzioni originali (escludi i moduli che testiamo con versioni semplificate)
tic("Caricamento moduli originali")
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
exclude_patterns <- c("06l_ligand_receptor_interactions\\.R", 
                     "06m_temporal_dynamics\\.R", 
                     "06n_alternative_splicing\\.R", 
                     "06o_anisotropic_patterns\\.R", 
                     "06p_3d_microenvironment\\.R")

for (file in sort(files)) {
  # Salta i file che testiamo con versioni semplificate
  if (!any(sapply(exclude_patterns, function(p) grepl(p, file)))) {
    cat("Caricamento:", basename(file), "\n")
    source(file)
  }
}
toc()

# Carica le versioni semplificate per i test
tic("Caricamento moduli semplificati")
simplified_files <- c(
  "R/functions/06l_ligand_receptor_interactions_simple.R",
  "R/functions/06m_temporal_dynamics_simple.R",
  "R/functions/06n_alternative_splicing_simple.R",
  "R/functions/06o_anisotropic_patterns_simple.R",
  "R/functions/06p_3d_microenvironment_simple.R"
)

for (file in simplified_files) {
  cat("Caricamento:", basename(file), "\n")
  source(file)
}
toc()

# Esecuzione dei test specifici per le funzioni semplificate
tic("Esecuzione tests semplificati")

# Esegui i test per i 5 moduli semplificati
test_files <- c(
  "tests/testthat/test-ligand_receptor_interactions.R",
  "tests/testthat/test-temporal_dynamics.R",
  "tests/testthat/test-alternative_splicing.R",
  "tests/testthat/test-anisotropic_patterns.R",
  "tests/testthat/test-3d_microenvironment.R"
)

for (test_file in test_files) {
  cat("Esecuzione test:", basename(test_file), "\n")
  testthat::test_file(test_file)
}

toc()