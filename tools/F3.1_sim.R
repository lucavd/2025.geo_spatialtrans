# tools/F3.1_sim.R - F3.1: simulate_visium_hd_pieces() sulle aree dello Step 1 (livelli U6 di S1.6), un file per (area, seme).
# Uso: Rscript --vanilla tools/F3.1_sim.R <out_dir> <seed_base> <tags comma-separated> <workers> [k list]
# Seme dell'esecuzione (area, k): tabella Mersenne-Twister estratta da seed_base sulle 30 aree ordinate (k = 0..9).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos); library(parallel) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
# Profili: quelli di S2.2 (predefiniti); per il braccio di CP-3 la variabile d'ambiente F31_PROFILES_DIR indica i profili dal solo nucleo.
if (nzchar(Sys.getenv("F31_PROFILES_DIR"))) options(geo.spatialtrans.profiles_dir = Sys.getenv("F31_PROFILES_DIR"))
a <- commandArgs(TRUE); out_dir <- a[1]; seed_base <- as.integer(a[2]); tags <- strsplit(a[3], ",")[[1]]; W <- as.integer(a[4])
ks <- if (length(a) >= 5) as.integer(strsplit(a[5], ",")[[1]]) else 0L
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
lay_dir <- "/mnt/micron/geo_spatialtrans/S1.6/layers"
all_tags <- sub("_U6\\.rds$", "", sort(list.files(lay_dir, pattern = "_U6\\.rds$")))
set.seed(seed_base, kind = "Mersenne-Twister")
seed_tab <- matrix(sample.int(.Machine$integer.max, length(all_tags) * 10), ncol = 10, dimnames = list(all_tags, NULL))
jobs <- expand.grid(tag = tags, k = ks, stringsAsFactors = FALSE)
run1 <- function(i) {
  tag <- jobs$tag[i]; k <- jobs$k[i]; s <- seed_tab[tag, k + 1L]
  fo <- file.path(out_dir, sprintf("%s_k%d.rds", tag, k)); if (file.exists(fo)) return(paste("skip", tag, k))
  L <- readRDS(file.path(lay_dir, paste0(tag, "_U6.rds")))
  t0 <- proc.time()[["elapsed"]]
  out <- simulate_visium_hd_pieces(L, seed = s, frame_open = TRUE, workers = 1L, assemble = "full")
  out$f31 <- list(tag = tag, k = k, seed = s, seed_base = seed_base, elapsed = proc.time()[["elapsed"]] - t0,
                  archetype = visium_hd_archetype(L), profiles_dir = getOption("geo.spatialtrans.profiles_dir", "default"))
  saveRDS(out, paste0(fo, ".tmp")); file.rename(paste0(fo, ".tmp"), fo)
  sprintf("%s k%d seed %d: %.0f s", tag, k, s, out$f31$elapsed)
}
Sys.setenv(OMP_NUM_THREADS = "1"); if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
res <- mclapply(seq_len(nrow(jobs)), function(i) tryCatch(run1(i), error = function(e) paste("ERROR", jobs$tag[i], conditionMessage(e))),
                mc.cores = W, mc.preschedule = FALSE)
cat(unlist(res), sep = "\n")
o <- readRDS(file.path(out_dir, sprintf("%s_k%d.rds", jobs$tag[1], jobs$k[1])))
cat("names:", names(o), "\ntruth:", names(o$truth), "\n")
for (nm in c("overlap", "origin", "cells")) { cat(nm, ":\n"); print(head(o$truth[[nm]], 3)) }
cat("bins:\n"); print(head(o$bins, 3)); cat("dim counts_2um", dim(o$counts_2um), "\n"); print(table(o$truth$cells$cell_type))
