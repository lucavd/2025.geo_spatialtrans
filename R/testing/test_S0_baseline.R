#!/usr/bin/env Rscript
# R/testing/test_S0_baseline.R — S0: confronto dei run di full_test.R (check C-S0.3, C-S0.4, C-S0.5, controprova 1)
# Prerequisiti (prodotti da R/testing/run_S0.sh):
#   R/testing/full_test_result_S0_run1.rds, _S0_run2.rds (seed 42), _S0_seed43.rds (seed 43)
#   results/full_simulation_data_S0_run1.rds, _S0_run2.rds, _S0_seed43.rds
# Uso: Rscript --vanilla R/testing/test_S0_baseline.R
lib <- file.path(getwd(), "renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(Matrix))
status <- function(ok) ifelse(isTRUE(ok), "PASS", "FAIL")
rd <- function(tag) readRDS(sprintf("R/testing/full_test_result_%s.rds", tag))
dp <- function(tag) sprintf("results/full_simulation_data_%s.rds", tag)
r1 <- rd("S0_run1"); r2 <- rd("S0_run2"); r3 <- rd("S0_seed43")
out <- list()
add <- function(id, ok, detail) { out[[length(out) + 1]] <<- data.frame(check = id, esito = status(ok), dettaglio = detail); cat(sprintf("[%s] %s — %s\n", status(ok), id, detail)) }

# C-S0.3: entrambi i run seed 42 hanno status PASS (7/7 asserzioni interne)
add("C-S0.3 run1 status", identical(r1$status, "PASS"), sprintf("status=%s, seed=%s", r1$status, r1$config$random_seed))
add("C-S0.3 run2 status", identical(r2$status, "PASS"), sprintf("status=%s, seed=%s", r2$status, r2$config$random_seed))

# C-S0.4: uguaglianza esatta dei riassunti fra run1 e run2 (seed 42)
keys <- c("n_cells_generated", "n_clusters_found", "matrix_dims", "sparsity_pct", "umi_mean", "umi_median", "umi_cv")
for (k in keys) add(paste("C-S0.4", k), identical(r1$results[[k]], r2$results[[k]]),
                    sprintf("run1=%s | run2=%s", paste(r1$results[[k]], collapse = "x"), paste(r2$results[[k]], collapse = "x")))
md5 <- tools::md5sum(c(dp("S0_run1"), dp("S0_run2"), dp("S0_seed43")))
add("C-S0.4 md5 full_simulation_data run1==run2", unname(md5[1]) == unname(md5[2]), sprintf("%s | %s", md5[1], md5[2]))
d1 <- readRDS(dp("S0_run1")); d2 <- readRDS(dp("S0_run2"))
add("C-S0.4 identical(expression) run1==run2", identical(d1$expression, d2$expression), sprintf("dim=%s, nnz=%d", paste(dim(d1$expression), collapse = "x"), Matrix::nnzero(d1$expression)))
add("C-S0.4 identical(coordinates) run1==run2", identical(d1$coordinates, d2$coordinates), sprintf("n=%d", nrow(d1$coordinates)))

# Controprova 1: seed 43 deve DIFFERIRE da seed 42
d3 <- readRDS(dp("S0_seed43"))
add("Controprova-1 seed43: status PASS", identical(r3$status, "PASS"), sprintf("status=%s, seed=%s", r3$status, r3$config$random_seed))
add("Controprova-1 seed43: md5 diverso da run1", unname(md5[3]) != unname(md5[1]), sprintf("%s vs %s", md5[3], md5[1]))
add("Controprova-1 seed43: expression diversa da run1", !identical(d3$expression, d1$expression),
    sprintf("n_cells run1=%d seed43=%d; sparsita run1=%.1f seed43=%.1f; umi_mean run1=%d seed43=%d",
            r1$results$n_cells_generated, r3$results$n_cells_generated, r1$results$sparsity_pct, r3$results$sparsity_pct, r1$results$umi_mean, r3$results$umi_mean))

# C-S0.5: pre-registrato NON riproducibile — il legacy full_test_result.rds (da full_test_visHD.R) differisce per config e n_cells
leg <- readRDS("R/testing/legacy/full_test_result.rds")
add("C-S0.5 legacy n_cells != run1 n_cells (atteso diverso)", leg$results$n_cells_generated != r1$results$n_cells_generated,
    sprintf("legacy=%d (cfg n_cells=%s, px=%s, grid=%s mm) vs run1=%d (cfg n_cells=%s, px=%s, grid=%s mm)",
            leg$results$n_cells_generated, leg$config$n_cells, leg$config$pixel_size_um, leg$config$fixed_grid_width_mm,
            r1$results$n_cells_generated, r1$config$n_cells, r1$config$pixel_size_um, r1$config$fixed_grid_width_mm))
legd <- readRDS("R/testing/legacy/full_simulation_data.rds")
add("C-S0.5 legacy full_simulation_data seed != 42 (atteso 1200)", identical(as.integer(legd$config$random_seed), 1200L), sprintf("seed legacy=%s", legd$config$random_seed))

tab <- do.call(rbind, out)
write.csv(tab, "results/S0_baseline_check.csv", row.names = FALSE)
summ <- data.frame(run = c("S0_run1", "S0_run2", "S0_seed43", "legacy_visHD"),
                   seed = c(r1$config$random_seed, r2$config$random_seed, r3$config$random_seed, leg$config$random_seed),
                   n_cells = c(r1$results$n_cells_generated, r2$results$n_cells_generated, r3$results$n_cells_generated, leg$results$n_cells_generated),
                   sparsity_pct = c(r1$results$sparsity_pct, r2$results$sparsity_pct, r3$results$sparsity_pct, leg$results$sparsity_pct),
                   umi_mean = c(r1$results$umi_mean, r2$results$umi_mean, r3$results$umi_mean, leg$results$umi_mean),
                   umi_median = c(r1$results$umi_median, r2$results$umi_median, r3$results$umi_median, leg$results$umi_median),
                   umi_cv = c(r1$results$umi_cv, r2$results$umi_cv, r3$results$umi_cv, leg$results$umi_cv),
                   elapsed_min = c(r1$results$elapsed_minutes, r2$results$elapsed_minutes, r3$results$elapsed_minutes, leg$results$elapsed_minutes),
                   md5_data = c(unname(md5), NA))
write.csv(summ, "results/S0_runs_summary.csv", row.names = FALSE)
print(summ)
cat(sprintf("\nTotale: %d PASS, %d FAIL\n", sum(tab$esito == "PASS"), sum(tab$esito == "FAIL")))
if (any(tab$esito == "FAIL")) quit(status = 1)
