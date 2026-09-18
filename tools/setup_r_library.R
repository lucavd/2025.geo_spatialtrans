#!/usr/bin/env Rscript
# tools/setup_r_library.R — S0 (2026-09-18)
# Costruisce la libreria R di progetto da binari Posit PPM (Ubuntu noble, R 4.6).
# Uso: Rscript --vanilla tools/setup_r_library.R
# Idempotente: reinstalla solo i pacchetti mancanti.
lib <- file.path(getwd(), "renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(
  repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"),
  HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(),
                          paste(getRversion(), R.version$platform, R.version$arch, R.version$os)),
  Ncpus = 40,
  install.packages.compile.from.source = "never"
)
# Set A: richiesti dalla roadmap per Step 1 (S0)
pk_step1 <- c("sf", "deldir", "polylabelr", "imager", "testthat", "quarto")
# Set B: dipendenze della pipeline attuale (R/*.R, R/testing/*.R, DESCRIPTION Imports/Suggests)
pk_pipeline <- c("Matrix", "ClusterR", "dplyr", "sp", "gstat", "png", "pbapply",
                 "future", "future.apply", "scales", "RColorBrewer", "ggplot2", "tidyr",
                 "viridis", "reshape2", "hexbin", "gridExtra", "pheatmap",
                 "spatstat", "spdep", "fields", "MASS", "cluster", "tictoc",
                 "tidyverse", "rmarkdown", "knitr", "withr", "renv")
pk <- unique(c(pk_step1, pk_pipeline))
t0 <- Sys.time()
missing <- setdiff(pk, rownames(installed.packages(lib.loc = lib)))
cat("Libreria:", lib, "\nDa installare:", length(missing), "pacchetti\n")
if (length(missing)) install.packages(missing, lib = lib, dependencies = c("Depends", "Imports", "LinkingTo"))
cat("Tempo installazione:", round(as.numeric(difftime(Sys.time(), t0, units = "secs"))), "s\n")

# --- Passo 2 (S0): libreria autosufficiente. Carica tutti i namespace richiesti e installa nella libreria
#     di progetto quelli ancora risolti dalla site-library di sistema (esclusi i pacchetti base di R).
for (p in pk) suppressWarnings(suppressMessages(loadNamespace(p)))
base_pk <- rownames(installed.packages(priority = "base"))
ns <- setdiff(loadedNamespaces(), base_pk)
ns_path <- vapply(ns, function(p) dirname(getNamespaceInfo(p, "path")), character(1))
from_sys <- ns[!startsWith(normalizePath(ns_path), normalizePath(lib))]
cat("Namespace risolti fuori dalla libreria di progetto:", length(from_sys), "\n")
if (length(from_sys)) install.packages(from_sys, lib = lib, dependencies = c("Depends", "Imports", "LinkingTo"))
cat("Pacchetti nella libreria di progetto:", length(rownames(installed.packages(lib.loc = lib))), "\n")

# --- Passo 3 (S0): chiusura ricorsiva delle dipendenze (Depends/Imports/LinkingTo) dal repository PPM,
#     così la libreria è autosufficiente anche per i namespace caricati pigramente a runtime (es. labeling).
ap <- available.packages()
base_pk <- rownames(installed.packages(priority = "base"))
closure <- unique(unlist(tools::package_dependencies(pk, db = ap, which = c("Depends", "Imports", "LinkingTo"), recursive = TRUE)))
closure <- setdiff(unique(c(pk, closure)), base_pk)
missing3 <- setdiff(closure, rownames(installed.packages(lib.loc = lib)))
cat("Chiusura ricorsiva:", length(closure), "pacchetti; mancanti nella libreria di progetto:", length(missing3), "\n")
if (length(missing3)) install.packages(missing3, lib = lib, dependencies = FALSE)
cat("Pacchetti nella libreria di progetto dopo il passo 3:", length(rownames(installed.packages(lib.loc = lib))), "\n")

# Verifica di caricamento (in questo processo; la verifica indipendente e' in R/testing/test_S0_library.R)
res <- vapply(pk, function(p) suppressWarnings(suppressMessages(requireNamespace(p, lib.loc = lib, quietly = TRUE))), logical(1))
print(data.frame(pkg = pk, loads = res, version = vapply(pk, function(p) tryCatch(as.character(packageVersion(p, lib.loc = lib)), error = function(e) NA_character_), character(1)), row.names = NULL))
cat("Caricano:", sum(res), "/", length(pk), "\n")
if (requireNamespace("sf", lib.loc = lib, quietly = TRUE)) print(sf::sf_extSoftVersion())
