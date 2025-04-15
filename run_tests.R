# Script semplificato per caricare i moduli ed eseguire i test

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

# Carica tutti i file dei moduli
tic("Caricamento moduli")
files <- list.files("R/functions", full.names = TRUE, pattern = "\\.R$")
for (file in sort(files)) {
  cat("Caricamento:", basename(file), "\n")
  source(file)
}
toc()

# Esecuzione dei test
tic("Esecuzione tests")
testthat::test_dir("tests/testthat/")
toc()
