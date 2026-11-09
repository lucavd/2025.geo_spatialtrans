# tools/S2.5_PH1_memory.R - S2.5, analisi post hoc PH-1 (dichiarata dopo i risultati): memoria e tempo dell'area piu' carica (A2_r4, main_r2 k = 0)
# con blocchi di 4, 2 e 1 milione di trascritti; un processo per esecuzione, in sequenza. Uso: Rscript --vanilla tools/S2.5_PH1_memory.R <ct> <out.json>
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
data.table::setDTthreads(1L)
a <- commandArgs(TRUE); ct <- as.numeric(a[1]); fo <- a[2]
tag <- "A2_r4"; A <- "A2"
lay <- sort(list.files("/mnt/micron/geo_spatialtrans/S1.6/layers", pattern = "_U6\\.rds$", full.names = TRUE)); tags <- sub("_U6\\.rds$", "", basename(lay)); ti <- match(tag, tags)
RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(20261109L); SEEDS <- matrix(sample.int(1e8, 2L * length(tags)), length(tags), 2L)
R2 <- read.csv("results/S2.5/S2.5_R2_depth.csv")
L <- readRDS(lay[ti]); hd <- readRDS(file.path("/mnt/micron/geo_spatialtrans/S2.3/panel/hd", paste0(tag, ".rds")))
hwm <- function() { s <- readLines("/proc/self/status"); as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", s[grepl("^VmHWM", s)])) / 1024 }
m0 <- hwm(); t0 <- proc.time()[["elapsed"]]
out <- simulate_visium_hd(L, seed = SEEDS[ti, 1L], archetype = A, frame_open = TRUE, hd = hd, depth_factor = R2$f_new[R2$archetype == A], chunk_transcripts = ct)
res <- list(tag = tag, chunk_transcripts = ct, chunk_size = out$metadata$chunk_size, n_generated = out$metadata$n_generated, peak_mb = hwm(), peak_mb_before = m0,
            elapsed_s = proc.time()[["elapsed"]] - t0, stages = as.list(out$metadata$elapsed_s),
            size_out_mb = as.numeric(object.size(out)) / 2^20, size_counts2_mb = as.numeric(object.size(out$counts_2um)) / 2^20,
            load = readLines("/proc/loadavg"))
writeLines(jsonlite::toJSON(res, auto_unbox = TRUE, digits = NA), fo); print(unlist(res[c("chunk_transcripts", "chunk_size", "peak_mb", "elapsed_s", "size_out_mb")]))
