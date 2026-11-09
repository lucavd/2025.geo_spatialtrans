# tools/S2.5_panel.R - S2.5: simulatore completo simulate_visium_hd() sulle 30 aree dello Step 1 e controlli del pannello (G-2).
# Definizioni: docs/steps/S2.5.md («Il pannello», cancelli, previsioni). Uso (root del repo):
#   Rscript --vanilla tools/S2.5_panel.R <variante> [outdir] [workers]   |   ... check <outdir>
#   varianti: main (2 semi), depth_off (1 seme, stessi semi di main k = 0, fattore di profondita' 1; controprova),
#             explore (un'area per tessuto, seme dell'esplorazione), s24_identity (2 aree con i semi e i blocchi di S2.4: continuita'),
#             main_r2 (2 semi, stessi di main, fattore di profondita' della regola R-2 da results/S2.5/S2.5_R2_depth.csv o S25_R2)
# Scrive <outdir>/<variante>/<A>_<roi>_k<k>.bin (UMI per quadratino di 2 um, int32, ordine di bin_id; stesso formato di S2.3 e S2.4),
# .g2.bin (geni rilevati per quadratino), .u8.bin e .g8.bin (blocchi di 8 um, ordine di bin8_id), _out.rds (uscita senza le matrici),
# _counts2.rds (matrice a 2 um, solo main), .json (controlli, tempi, memoria). Griglie di S2.1 lette da S2.3 (sola lettura, si controllano).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
a <- commandArgs(TRUE); VAR <- a[1]; OUT <- if (length(a) >= 2) a[2] else "/mnt/micron/geo_spatialtrans/S2.5/panel"; W <- if (length(a) >= 3) as.integer(a[3]) else 4L
HD23 <- "/mnt/micron/geo_spatialtrans/S2.3/panel/hd"
lay <- sort(list.files("/mnt/micron/geo_spatialtrans/S1.6/layers", pattern = "_U6\\.rds$", full.names = TRUE))
tags <- sub("_U6\\.rds$", "", basename(lay))
# semi dell'utente per (area, k), estratti con Mersenne-Twister dal seme di base (DEC-057); l'orchestratore ne estrae i semi dei passi
RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(20261109L); SEEDS <- matrix(sample.int(1e8, 2L * length(tags)), length(tags), 2L)
set.seed(20261110L); SEEDS_EX <- sample.int(1e8, length(tags))
EXPLORE <- c("A1_r1", "A2_r1", "A3_r1", "A4_f1", "A5_r1", "A6_r1")
# continuita' con S2.4: semi dei passi del pannello di S2.4 (tools/S2.4_panel.R) e blocchi di 2000 cellule
RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(20261114L); S24 <- array(sample.int(1e8, 2L * length(tags) * 4L), c(4L, length(tags), 2L))
IDENT <- c("A2_r1", "A6_r1")
peak_mb <- function() { s <- readLines("/proc/self/status"); as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", s[grepl("^VmHWM", s)])) / 1024 }

if (VAR == "check") {
  J <- list()
  for (v in c("main", "depth_off", "main_r2")) for (f in list.files(file.path(OUT, v), pattern = "_k[01]\\.json$", full.names = TRUE)) J[[length(J) + 1]] <- jsonlite::fromJSON(f)
  nz <- function(x) if (is.null(x) || length(x) == 0) NA_real_ else x
  G2 <- do.call(rbind, lapply(J, function(j) data.frame(tag = j$tag, archetype = j$archetype, k = j$k, variant = j$variant, M = j$M, chunk_size = j$chunk_size,
    conservation = j$conservation_cells && j$conservation_genes && j$conservation_bins && j$conservation_8um && j$truth8_ok && j$profiles_match,
    n_na = j$n_na, z_captured = j$z_captured, z_background = nz(j$z_background), elapsed_total_s = j$elapsed_total_s, peak_mb = j$peak_mb)))
  write.csv(G2, file.path(OUT, "G2.csv"), row.names = FALSE)
  zt <- stats::qnorm(1 - 0.0005 / (nrow(G2) + sum(!is.na(G2$z_background))))
  cat(sprintf("G-2: %d esecuzioni; conservazione %d; NA %d; max |z| catturati %.2f, fondo %.2f; soglia %.2f\n",
              nrow(G2), sum(G2$conservation), sum(G2$n_na), max(abs(G2$z_captured)), max(abs(G2$z_background), na.rm = TRUE), zt))
  quit(status = if (all(G2$conservation) && sum(G2$n_na) == 0 && max(abs(G2$z_captured)) < zt && max(abs(G2$z_background), na.rm = TRUE) < zt) 0L else 1L)
}
stopifnot(VAR %in% c("main", "depth_off", "explore", "s24_identity", "main_r2"))
R2TAB <- if (VAR == "main_r2") read.csv(Sys.getenv("S25_R2", "results/S2.5/S2.5_R2_depth.csv")) else NULL
dir.create(file.path(OUT, VAR), recursive = TRUE, showWarnings = FALSE)
jobs <- switch(VAR, main_r2 = expand.grid(t = seq_along(tags), k = 0:1), main = expand.grid(t = seq_along(tags), k = if (nzchar(Sys.getenv("S25_KS"))) as.integer(strsplit(Sys.getenv("S25_KS"), ",")[[1]]) else 0:1),
               depth_off = expand.grid(t = seq_along(tags), k = 0L), explore = expand.grid(t = match(EXPLORE, tags), k = 0L),
               s24_identity = expand.grid(t = match(IDENT, tags), k = 0L))
one <- function(q) {
  ti <- jobs$t[q]; k <- jobs$k[q]; tag <- tags[ti]; A <- strsplit(tag, "_")[[1]][1]
  fo <- file.path(OUT, VAR, sprintf("%s_k%d", tag, k))
  if (file.exists(paste0(fo, ".json"))) return(fo)
  t0 <- proc.time()[["elapsed"]]
  L <- readRDS(lay[ti]); hd <- readRDS(file.path(HD23, paste0(tag, ".rds")))
  args <- list(cell_layer = L, archetype = A, frame_open = TRUE, hd = hd, workers = 1L, chunk_transcripts = 4e6)   # valore della regola R-1 con cui e' girato il pannello
  if (VAR == "s24_identity") { args$component_seeds <- setNames(S24[, ti, 1L], c("stripes", "profiles", "place", "capture")); args$chunk_size <- 2000L }
  else args$seed <- if (VAR == "explore") SEEDS_EX[ti] else SEEDS[ti, k + 1L]
  # profondita' esplicita per ogni variante (il predefinito dell'orchestratore e' cambiato dopo R-2): main, explore, s24_identity con il fattore
  # di S2.2 dei profili, con cui sono girate
  if (VAR %in% c("main", "explore", "s24_identity")) args$depth_factor <- load_cell_profiles(A)$depth_factor
  if (VAR == "depth_off") args$depth_factor <- 1
  if (VAR == "main_r2") { args$depth_factor <- R2TAB$f_new[R2TAB$archetype == A]; stopifnot(length(args$depth_factor) == 1L, is.finite(args$depth_factor)) }
  out <- do.call(simulate_visium_hd, args)
  te <- proc.time()[["elapsed"]]
  g <- out$metadata$grid; nb <- g$n_row * g$n_col
  X2 <- out$counts_2um; X8 <- out$counts_8um
  u2 <- as.integer(round(Matrix::colSums(X2))); g2 <- as.integer(diff(X2@p))
  u8 <- as.integer(round(Matrix::colSums(X8))); g8 <- as.integer(diff(as(X8, "CsparseMatrix")@p))
  tr <- out$truth; cl <- tr$cells; gl <- tr$genes
  ob <- numeric(nb); tb <- rowsum(tr$origin$n, tr$origin$bin_id); ob[as.integer(rownames(tb))] <- tb[, 1]
  bgv <- numeric(nb); bgv[tr$background$bin_id] <- tr$background$n
  cvec <- out$bins$c; C <- max(cvec); pk <- if (out$metadata$layers$stripes) cvec / C else rep(1, nb)
  Ecap <- sum(tr$arrived_2um * pk); Vcap <- sum(tr$arrived_2um * pk * (1 - pk))
  Ebg <- if (out$metadata$layers$background) tr$layers$background$b * sum(cvec) else 0
  b8_of <- function(id) { i <- (id - 1L) %/% g$n_col; j <- (id - 1L) %% g$n_col; as.integer((i %/% g$block) * (g$n_col %/% g$block) + j %/% g$block + 1L) }
  n8 <- (g$n_row %/% g$block) * (g$n_col %/% g$block)
  tot8 <- function(id, n) { v <- numeric(n8); s <- rowsum(as.numeric(n), b8_of(id)); v[as.integer(rownames(s))] <- s[, 1]; v }
  t8 <- function(id8, n) { v <- numeric(n8); s <- rowsum(as.numeric(n), id8); v[as.integer(rownames(s))] <- s[, 1]; v }
  truth8_ok <- isTRUE(all(tot8(tr$origin$bin_id, tr$origin$n) == t8(tr$origin_8um$bin8_id, tr$origin_8um$n))) &&
               isTRUE(all(tot8(tr$origin$bin_id, tr$origin$n_born_nucleus) == t8(tr$origin_8um$bin8_id, tr$origin_8um$n_born_nucleus))) &&
               isTRUE(all(abs(tot8(tr$overlap$bin_id, tr$overlap$area_territory_um2) - t8(tr$overlap_8um$bin8_id, tr$overlap_8um$area_territory_um2)) < 1e-9)) &&
               isTRUE(all(tot8(tr$background$bin_id, tr$background$n) == t8(tr$background_8um$bin8_id, tr$background_8um$n)))
  chk <- list(tag = tag, archetype = A, k = k, variant = VAR, seed = if (is.null(args$seed)) NA else args$seed, component_seeds = as.list(out$metadata$component_seeds),
    b = tr$layers$background$b, depth_factor = tr$layers$depth$depth_factor,
    M = out$metadata$overgeneration, chunk_size = out$metadata$chunk_size, n_cells = nrow(L$cell_df),
    n_generated = sum(cl$n_generated), n_lost = sum(cl$n_lost), n_thinned = sum(cl$n_thinned), n_captured = sum(cl$n_captured), n_background = sum(bgv),
    n_target_R4 = sum(cl$n_expected_captured),
    profiles_match = isTRUE(all(cl$n_generated == as.numeric(Matrix::colSums(tr$cell_counts)))),
    conservation_cells = all(cl$n_generated == cl$n_captured + cl$n_thinned + cl$n_lost),
    conservation_genes = all(gl$n_generated == gl$n_captured + gl$n_thinned + gl$n_lost) && all(as.numeric(Matrix::rowSums(X2)) == gl$n_captured + gl$n_background),
    conservation_bins = all(u2 == ob + bgv) && sum(tr$arrived_2um) == sum(cl$n_generated) - sum(cl$n_lost),
    conservation_8um = { D8 <- X8 - sum_to_8um(X2, g); D8 <- as(D8, "CsparseMatrix"); length(D8@x) == 0L || all(D8@x == 0) } && all(u8 == as.numeric(sum_to_8um(matrix(u2, 1), g))),
    truth8_ok = truth8_ok,
    n_in_tissue = sum(hd$bins$in_tissue), n_tissue_no_territory = sum(hd$bins$in_tissue & !(hd$bins$bin_id %in% hd$truth$overlap$bin_id)),
    n_na = sum(is.na(u2)) + sum(is.na(tr$origin$n)) + sum(is.na(cl$n_captured)),
    E_captured = Ecap, z_captured = (sum(cl$n_captured) - Ecap) / sqrt(max(Vcap, 1e-12)),
    E_background = Ebg, z_background = if (Ebg > 0) (sum(bgv) - Ebg) / sqrt(Ebg) else NA,
    grid = list(x0 = g$x0, y0 = g$y0, n_row = g$n_row, n_col = g$n_col, bin_um = g$bin_um, block = g$block),
    elapsed = as.list(out$metadata$elapsed_s), elapsed_total_s = proc.time()[["elapsed"]] - t0, elapsed_sim_s = te - t0, peak_mb = peak_mb())
  writeBin(u2, paste0(fo, ".bin"), size = 4L); writeBin(g2, paste0(fo, ".g2.bin"), size = 4L)
  writeBin(u8, paste0(fo, ".u8.bin"), size = 4L); writeBin(g8, paste0(fo, ".g8.bin"), size = 4L)
  if (VAR %in% c("main", "main_r2")) saveRDS(X2, paste0(fo, "_counts2.rds"))
  out$counts_2um <- NULL; out$counts_8um <- NULL; out$truth$overlap <- NULL; saveRDS(out, paste0(fo, "_out.rds"))
  writeLines(jsonlite::toJSON(chk, auto_unbox = TRUE, digits = NA), paste0(fo, ".json.tmp")); file.rename(paste0(fo, ".json.tmp"), paste0(fo, ".json"))
  cat(sprintf("[%s k%d %s] M %.2f blocchi %d generati %.0f catturati %.0f fondo %.0f, %.0f s, %.0f MB\n", tag, k, VAR, chk$M, chk$chunk_size, chk$n_generated, chk$n_captured, chk$n_background, chk$elapsed_total_s, chk$peak_mb))
  fo
}
todo <- which(!file.exists(file.path(OUT, VAR, sprintf("%s_k%d.json", tags[jobs$t], jobs$k))))
r <- parallel::mclapply(todo, function(q) tryCatch(one(q), error = function(e) paste("ERRORE", tags[jobs$t[q]], conditionMessage(e))), mc.cores = W, mc.preschedule = FALSE)
print(unlist(r)[grepl("^ERRORE", unlist(r))])
