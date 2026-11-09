# tools/S2.5_PH2_timing.R - S2.5, analisi post hoc PH-2 (dichiarata dopo i risultati): tempo e memoria delle 5 aree del tumore del colon con i
# predefiniti finali dell'orchestratore (fattore di profondita' di R-2, 2e6 trascritti per blocco), un processo per area, a macchina libera.
# Uso: Rscript --vanilla tools/S2.5_PH2_timing.R <tag> <out.json>
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
data.table::setDTthreads(1L)
a <- commandArgs(TRUE); tag <- a[1]; fo <- a[2]; A <- strsplit(tag, "_")[[1]][1]
lay <- sort(list.files("/mnt/micron/geo_spatialtrans/S1.6/layers", pattern = "_U6\\.rds$", full.names = TRUE)); tags <- sub("_U6\\.rds$", "", basename(lay)); ti <- match(tag, tags)
RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(20261109L); SEEDS <- matrix(sample.int(1e8, 2L * length(tags)), length(tags), 2L)
L <- readRDS(lay[ti]); hd <- readRDS(file.path("/mnt/micron/geo_spatialtrans/S2.3/panel/hd", paste0(tag, ".rds")))
hwm <- function() { s <- readLines("/proc/self/status"); as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", s[grepl("^VmHWM", s)])) / 1024 }
l0 <- readLines("/proc/loadavg"); t0 <- proc.time()[["elapsed"]]
out <- simulate_visium_hd(L, seed = SEEDS[ti, 1L], frame_open = TRUE, hd = hd)
res <- list(tag = tag, chunk_size = out$metadata$chunk_size, chunk_transcripts = out$metadata$chunk_transcripts, depth_factor = out$truth$layers$depth$depth_factor,
            n_generated = out$metadata$n_generated, peak_mb = hwm(), elapsed_s = proc.time()[["elapsed"]] - t0, stages = as.list(out$metadata$elapsed_s),
            load_before = l0, load_after = readLines("/proc/loadavg"))
writeLines(jsonlite::toJSON(res, auto_unbox = TRUE, digits = NA), fo)
