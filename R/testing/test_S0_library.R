#!/usr/bin/env Rscript
# R/testing/test_S0_library.R — S0: verifica indipendente della libreria R di progetto (check C-S0.2)
# Uso (dalla root del repo): Rscript --vanilla R/testing/test_S0_library.R
# Stampa una riga PASS/FAIL per asserzione; exit status 1 se una FAIL.
lib <- file.path(getwd(), "renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(c(lib, .libPaths()))
pk_step1 <- c("sf", "deldir", "polylabelr", "imager", "testthat", "quarto")
pk_pipeline <- c("Matrix", "ClusterR", "dplyr", "sp", "gstat", "png", "pbapply",
                 "future", "future.apply", "scales", "RColorBrewer", "ggplot2", "tidyr",
                 "viridis", "reshape2", "hexbin", "gridExtra", "pheatmap",
                 "spatstat", "spdep", "fields", "MASS", "cluster", "tictoc",
                 "tidyverse", "rmarkdown", "knitr", "withr", "renv")
pk <- unique(c(pk_step1, pk_pipeline))
ip <- installed.packages()
rows <- lapply(pk, function(p) {
  ok <- suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))
  path <- tryCatch(find.package(p), error = function(e) NA_character_)
  in_proj <- !is.na(path) && startsWith(normalizePath(path), normalizePath(lib))
  ver <- tryCatch(as.character(packageVersion(p)), error = function(e) NA_character_)
  built <- if (p %in% rownames(ip)) ip[p, "Built"] else NA_character_
  data.frame(pkg = p, set = ifelse(p %in% pk_step1, "step1", "pipeline"), loads = ok,
             version = ver, in_project_lib = in_proj, built_R = built, stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)
status <- function(ok) ifelse(ok, "PASS", "FAIL")
for (i in seq_len(nrow(tab))) cat(sprintf("[%s] C-S0.2 load %-14s v%-10s lib=%s\n", status(tab$loads[i]), tab$pkg[i], tab$version[i], ifelse(tab$in_project_lib[i], "project", "system")))
ext <- sf::sf_extSoftVersion()
exp_ext <- c(GEOS = "3.12.1", GDAL = "3.8.4", PROJ = "9.4.0")
for (k in names(exp_ext)) cat(sprintf("[%s] C-S0.2 sf %s = %s (atteso %s)\n", status(identical(unname(ext[k]), exp_ext[[k]])), k, ext[k], exp_ext[[k]]))
cat(sprintf("[%s] C-S0.2 tutti i %d pacchetti caricano (%d/%d)\n", status(all(tab$loads)), length(pk), sum(tab$loads), length(pk)))
cat(sprintf("[INFO] pacchetti dalla libreria di progetto: %d/%d; dalla site-library di sistema: %d\n", sum(tab$in_project_lib), nrow(tab), sum(!tab$in_project_lib)))
cat(sprintf("[INFO] R %s; .libPaths() = %s\n", getRversion(), paste(.libPaths(), collapse = " | ")))
site_ok <- all(normalizePath(.libPaths()) %in% normalizePath(c(lib, .Library)))
cat(sprintf("[%s] C-S0.2 .libPaths() contiene solo libreria di progetto + libreria base di R (isolamento; richiede R_LIBS_SITE e R_LIBS_USER inesistenti come in run_S0.sh)\n", status(site_ok)))
dir.create("results", showWarnings = FALSE)
write.csv(tab, "results/S0_library_check.csv", row.names = FALSE)
if (!all(tab$loads)) quit(status = 1)
