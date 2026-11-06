# tools/S2.4_mem_stages.R - S2.4, analisi del dubbio sulla memoria (P-7): picco di memoria residente dopo ogni fase della catena,
# area A2_r1 seme k = 0, stessi semi del pannello. Uso: Rscript --vanilla tools/S2.4_mem_stages.R <tag> -> stampa una riga per fase.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
tag <- commandArgs(TRUE)[1]; A <- strsplit(tag, "_")[[1]][1]
hwm <- function() { s <- readLines("/proc/self/status"); as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", s[grepl("^VmHWM", s)])) / 1024 }
rss <- function() { s <- readLines("/proc/self/status"); as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", s[grepl("^VmRSS", s)])) / 1024 }
lay <- sort(list.files("/mnt/micron/geo_spatialtrans/S1.6/layers", pattern = "_U6\\.rds$", full.names = TRUE)); tags <- sub("_U6\\.rds$", "", basename(lay)); ti <- match(tag, tags)
J <- jsonlite::fromJSON(sprintf("/mnt/micron/geo_spatialtrans/S2.4/panel/main/%s_k0.json", tag)); sd_ <- J$seeds
rep_ <- function(s) cat(sprintf("%-28s picco %7.0f MB  residente %7.0f MB\n", s, hwm(), rss()))
L <- readRDS(lay[ti]); hd <- readRDS(sprintf("/mnt/micron/geo_spatialtrans/S2.3/panel/hd/%s.rds", tag)); rep_("livello e griglia")
cp <- capture_params(A); p <- transcript_params(A); p$w2 <- cp$w2; p$ell2 <- cp$ell2
St <- simulate_stripes(hd, cp, seed = sd_[1]); rep_("strisce")
pd <- load_cell_profiles(A); prof <- simulate_cell_profiles(L$cell_df[, c("cell_id", "cell_type")], pd, seed = sd_[2], depth_factor = pd$depth_factor * St$M); rep_("profili (x M)")
cat(sprintf("trascritti generati %.0f, M %.2f, dimensione profili %.0f MB\n", sum(prof$counts), St$M, as.numeric(object.size(prof)) / 2^20))
tx <- place_and_move_transcripts(L, hd, prof, p, seed = sd_[3], field = TRUE, frame_open = TRUE); rep_("posizione e spostamento")
cat(sprintf("terne %d, dimensione uscita %.0f MB\n", nrow(tx$transcripts), as.numeric(object.size(tx)) / 2^20))
rm(prof); invisible(gc())
cap <- capture_and_bin(tx, hd, St, list(b = cp$b), seed = sd_[4], overgeneration = St$M); rep_("cattura")
cat(sprintf("dimensione uscita della cattura %.0f MB; identico al pannello: %s\n", as.numeric(object.size(cap)) / 2^20,
            identical(as.integer(round(Matrix::colSums(cap$counts_2um))), readBin(sprintf("/mnt/micron/geo_spatialtrans/S2.4/panel/main/%s_k0.bin", tag), "integer", n = hd$grid$n_row * hd$grid$n_col, size = 4L))))
