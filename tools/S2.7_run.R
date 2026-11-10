# tools/S2.7_run.R - S2.7: esecuzioni del calcolo (docs/steps/S2.7.md). Uso (root del repo):
#   Rscript --vanilla tools/S2.7_run.R new <layer.rds> <seed> <piece_um|Inf> <workers> <outdir> <tag> <frame_open 0|1> <full|light> [k] [variant]
#   Rscript --vanilla tools/S2.7_run.R old <area_tag> <seed> <k> <outdir>          motore di S2.5 come la variante main_r2 del pannello di S2.5
# Uscite in <outdir>/<tag>[_k<k>]: .bin (UMI per quadratino di 2 um, int32), .g2.bin (geni per quadratino), .u8.bin, .g8.bin (blocchi di 8 um),
# .json (tag, archetype, k, griglia, tempi per fase, memoria, controlli, metadati), .md5 (impronte delle uscite, JSON), _c8.rds (light).
# Con "new" i file dei pezzi stanno in <outdir>/<tag>_pieces/ (rilanciabile); si cancellano a esecuzione riuscita salvo S27_KEEP=1.
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
peak_mb <- function() { s <- readLines("/proc/self/status"); as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", s[grepl("^VmHWM", s)])) / 1024 }
a <- commandArgs(TRUE); MODE <- a[1]
md5_of <- function(x) { f <- tempfile(); saveRDS(x, f, compress = FALSE); on.exit(unlink(f)); unname(tools::md5sum(f)) }
write_bins <- function(fo, X2, X8) {
  writeBin(as.integer(round(Matrix::colSums(X2))), paste0(fo, ".bin"), size = 4L); writeBin(as.integer(diff(X2@p)), paste0(fo, ".g2.bin"), size = 4L)
  X8 <- as(X8, "CsparseMatrix")
  writeBin(as.integer(round(Matrix::colSums(X8))), paste0(fo, ".u8.bin"), size = 4L); writeBin(as.integer(diff(X8@p)), paste0(fo, ".g8.bin"), size = 4L)
}
wjson <- function(x, f) { writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, digits = NA, null = "null"), paste0(f, ".tmp")); file.rename(paste0(f, ".tmp"), f) }
if (MODE == "new") {
  lay <- a[2]; seed <- as.numeric(a[3]); pu <- as.numeric(a[4]); W <- as.integer(a[5]); OUT <- a[6]; tag <- a[7]; fo_ <- a[8] == "1"; asm <- a[9]
  k <- if (length(a) >= 10) as.integer(a[10]) else NA_integer_; variant <- if (length(a) >= 11) a[11] else "new"
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  base <- file.path(OUT, if (is.na(k)) tag else sprintf("%s_k%d", tag, k))
  if (file.exists(paste0(base, ".json"))) quit(status = 0L)
  pd <- paste0(base, "_pieces")
  t0 <- proc.time()[["elapsed"]]
  out <- simulate_visium_hd_pieces(lay, seed = seed, piece_um = pu, frame_open = fo_, out_dir = pd, workers = W, assemble = asm, verbose = TRUE)
  te <- proc.time()[["elapsed"]] - t0
  m <- out$metadata; g <- m$grid
  ps <- m$piece_stats
  j <- list(tag = tag, archetype = m$archetype, k = k, variant = variant, seed = seed, layer = lay, piece_um = m$piece_um, n_pieces = m$n_pieces, workers = W,
            frame_open = fo_, grid = list(x0 = g$x0, y0 = g$y0, n_row = g$n_row, n_col = g$n_col, bin_um = g$bin_um, block = g$block),
            n_cells = m$n_cells, overgeneration = m$overgeneration, n_bins_over_ceiling = m$n_bins_over_ceiling, truncated_share = m$truncated_share,
            n_generated = m$n_generated, n_captured = m$n_captured, n_background = m$n_background, n_lost = m$n_lost, n_sent = m$n_sent_between_pieces,
            n_fallback = m$n_fallback, fallback_cell_ids = as.list(m$fallback_cell_ids), checks = m$checks,
            n_in_tissue = sum(out$bins$in_tissue), tissue_area_um2 = sum(out$bins$tissue_area_um2),
            elapsed = as.list(m$elapsed_s), elapsed_total_s = te, main_peak_mb = peak_mb(), main_rss_at_fork_mb = m$main_rss_at_fork_mb,
            piece_peak_max_mb = max(ps$peak_mb, na.rm = TRUE), piece_time_sum_s = sum(ps$elapsed_s, na.rm = TRUE),
            stage_time_sum_s = as.list(tapply(ps$elapsed_s, ps$stage, sum, na.rm = TRUE)), stage_peak_max_mb = as.list(tapply(ps$peak_mb, ps$stage, max, na.rm = TRUE)),
            n_reused = sum(ps$reused), load = paste(readLines("/proc/loadavg"), collapse = ""))
  utils::write.csv(ps, paste0(base, "_pieces.csv"), row.names = FALSE)
  if (asm == "full") {
    write_bins(base, out$counts_2um, out$counts_8um)
    j$n_tissue_no_territory <- sum(out$bins$in_tissue & !(out$bins$bin_id %in% out$truth$overlap$bin_id))
    tr <- out$truth
    h <- list(counts_2um = md5_of(out$counts_2um), counts_8um = md5_of(out$counts_8um), bins = md5_of(out$bins), bins_8um = md5_of(out$bins_8um),
              overlap = md5_of(tr$overlap), overlap_8um = md5_of(tr$overlap_8um), origin = md5_of(tr$origin), origin_8um = md5_of(tr$origin_8um),
              background = md5_of(tr$background), background_8um = md5_of(tr$background_8um), arrived_2um = md5_of(tr$arrived_2um),
              cells = md5_of(tr$cells), genes = md5_of(tr$genes), layers = md5_of(tr$layers))
    wjson(h, paste0(base, ".md5"))
  } else {
    writeBin(as.integer(out$umi_2um), paste0(base, ".bin"), size = 4L)
    saveRDS(out$counts_8um, paste0(base, "_c8.rds"))
    utils::write.csv(out$truth$cells, gzfile(paste0(base, "_cells.csv.gz")), row.names = FALSE)
    wjson(list(u2 = md5_of(out$umi_2um), counts_8um = md5_of(out$counts_8um), cells = md5_of(out$truth$cells), genes = md5_of(out$truth$genes)), paste0(base, ".md5"))
  }
  wjson(j, paste0(base, ".json"))
  if (!nzchar(Sys.getenv("S27_KEEP"))) unlink(pd, recursive = TRUE)
  cat(sprintf("[%s] %s, %d pezzi, %.0f s, picco principale %.0f MB, picco pezzi %.0f MB, controlli %s\n", basename(base), m$archetype, m$n_pieces, te, j$main_peak_mb,
              j$piece_peak_max_mb, all(unlist(m$checks))))
} else if (MODE == "old") {
  tag <- a[2]; seed <- as.numeric(a[3]); k <- as.integer(a[4]); OUT <- a[5]; A <- strsplit(tag, "_")[[1]][1]
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  base <- file.path(OUT, sprintf("%s_k%d", tag, k)); if (file.exists(paste0(base, ".json"))) quit(status = 0L)
  t0 <- proc.time()[["elapsed"]]
  L <- readRDS(sprintf("/mnt/micron/geo_spatialtrans/S1.6/layers/%s_U6.rds", tag)); hd <- readRDS(sprintf("/mnt/micron/geo_spatialtrans/S2.3/panel/hd/%s.rds", tag))
  out <- simulate_visium_hd(L, seed = seed, archetype = A, frame_open = TRUE, hd = hd, workers = 1L, chunk_transcripts = 4e6)
  te <- proc.time()[["elapsed"]] - t0
  write_bins(base, out$counts_2um, out$counts_8um); g <- out$metadata$grid
  wjson(list(tag = tag, archetype = A, k = k, variant = "old", seed = seed, grid = list(x0 = g$x0, y0 = g$y0, n_row = g$n_row, n_col = g$n_col, bin_um = g$bin_um, block = g$block),
             M = out$metadata$overgeneration, n_in_tissue = sum(hd$bins$in_tissue), n_tissue_no_territory = sum(hd$bins$in_tissue & !(hd$bins$bin_id %in% hd$truth$overlap$bin_id)),
             n_captured = out$metadata$n_captured, elapsed_total_s = te, peak_mb = peak_mb()), paste0(base, ".json"))
  cat(sprintf("[%s old] M %.3f, %.0f s\n", basename(base), out$metadata$overgeneration, te))
} else stop("modo sconosciuto")
