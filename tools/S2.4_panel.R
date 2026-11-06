# tools/S2.4_panel.R - S2.4: catena completa (strisce, profili sovra-generati, posizione e spostamento con la seconda coda, cattura) sulle
# 30 aree dello Step 1 e cancello G-2. Definizioni: docs/steps/S2.4.md («Il pannello», G-2, P-7, correzione C-1 sui semi).
# Uso (root del repo): Rscript --vanilla tools/S2.4_panel.R <variante> [outdir] [workers]   |   ... check [outdir]
#   varianti: main (2 semi), stripes_off (1 seme, stessi semi di main k = 0; CP-1), tail2_off (1 seme, w2 = 0; CP-2)
# Scrive <outdir>/<variante>/<A>_<roi>_k<k>.bin (UMI per quadratino di 2 um, int32, ordine di bin_id), .json (griglia, controlli, tempi),
# _cap.rds (uscita di capture_and_bin senza la matrice a 2 um). Griglie di S2.1 lette da S2.3 (<S23>/hd/<A>_<roi>.rds, sola lettura).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
a <- commandArgs(TRUE); VAR <- a[1]; OUT <- if (length(a) >= 2) a[2] else "/mnt/micron/geo_spatialtrans/S2.4/panel"; W <- if (length(a) >= 3) as.integer(a[3]) else 4L
HD23 <- "/mnt/micron/geo_spatialtrans/S2.3/panel/hd"
lay <- sort(list.files("/mnt/micron/geo_spatialtrans/S1.6/layers", pattern = "_U6\\.rds$", full.names = TRUE))
tags <- sub("_U6\\.rds$", "", basename(lay))
# correzione C-1: semi di ogni chiamata estratti con Mersenne-Twister dal seme di base 20261114 (niente progressioni aritmetiche),
# in ordine fisso: k (0, 1) x area (1..30) x componente (1 strisce, 2 profili, 3 posizione e spostamento, 4 cattura)
RNGkind("Mersenne-Twister", "Inversion", "Rejection"); set.seed(20261114L)
SEEDS <- array(sample.int(1e8, 2L * length(tags) * 4L), c(4L, length(tags), 2L))
seed_of <- function(k, ti, comp) SEEDS[comp, ti, k + 1L]
peak_mb <- function() { s <- readLines("/proc/self/status"); as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", s[grepl("^VmHWM", s)])) / 1024 }

if (VAR == "check") {
  J <- list()
  for (v in c("main", "stripes_off", "tail2_off")) for (f in list.files(file.path(OUT, v), pattern = "_k[01]\\.json$", full.names = TRUE)) J[[length(J) + 1]] <- jsonlite::fromJSON(f)
  nz <- function(x) if (is.null(x) || length(x) == 0) NA_real_ else x
  G2 <- do.call(rbind, lapply(J, function(j) data.frame(tag = j$tag, archetype = j$archetype, k = j$k, variant = j$variant, M = j$M,
    conservation = j$conservation_cells && j$conservation_genes && j$conservation_bins && j$profiles_match, n_na = j$n_na,
    z_captured = j$z_captured, z_background = nz(j$z_background), elapsed_total_s = j$elapsed_total_s, peak_mb = j$peak_mb)))
  write.csv(G2, file.path(OUT, "G2.csv"), row.names = FALSE)
  zt <- stats::qnorm(1 - 0.0005 / nrow(G2))
  cat(sprintf("G-2: %d esecuzioni; conservazione %d; NA %d; max |z| catturati %.2f, fondo %.2f (soglia di Bonferroni sulle esecuzioni, 3 errori Monte Carlo: %.2f)\n",
              nrow(G2), sum(G2$conservation), sum(G2$n_na), max(abs(G2$z_captured)), max(abs(G2$z_background), na.rm = TRUE), 3))
  quit(status = if (all(G2$conservation) && sum(G2$n_na) == 0 && max(abs(G2$z_captured)) < 3 && max(abs(G2$z_background), na.rm = TRUE) < 3) 0L else 1L)
}
stopifnot(VAR %in% c("main", "stripes_off", "tail2_off"))
dir.create(file.path(OUT, VAR), recursive = TRUE, showWarnings = FALSE)
ks <- if (VAR == "main") { if (nzchar(Sys.getenv("S24_KS"))) as.integer(strsplit(Sys.getenv("S24_KS"), ",")[[1]]) else 0:1 } else 0L
jobs <- expand.grid(t = seq_along(tags), k = ks)
one <- function(q) {
  ti <- jobs$t[q]; k <- jobs$k[q]; tag <- tags[ti]; A <- strsplit(tag, "_")[[1]][1]
  fo <- file.path(OUT, VAR, sprintf("%s_k%d", tag, k))
  if (file.exists(paste0(fo, ".json"))) return(fo)
  t0 <- proc.time()[["elapsed"]]
  L <- readRDS(lay[ti]); hd <- readRDS(file.path(HD23, paste0(tag, ".rds")))
  cp <- capture_params(A); p <- transcript_params(A)
  p$w2 <- if (VAR == "tail2_off") 0 else cp$w2; p$ell2 <- cp$ell2
  ts <- proc.time()[["elapsed"]]
  St <- if (VAR == "stripes_off") NULL else simulate_stripes(hd, cp, seed = seed_of(k, ti, 1L))
  M <- if (is.null(St)) 1 else St$M
  tp <- proc.time()[["elapsed"]]
  pd <- load_cell_profiles(A)
  prof <- simulate_cell_profiles(L$cell_df[, c("cell_id", "cell_type")], pd, seed = seed_of(k, ti, 2L), depth_factor = pd$depth_factor * M, workers = 1L)
  tm <- proc.time()[["elapsed"]]
  tx <- place_and_move_transcripts(L, hd, prof, p, seed = seed_of(k, ti, 3L), field = TRUE, frame_open = TRUE, workers = 1L)
  tc <- proc.time()[["elapsed"]]
  cap <- capture_and_bin(tx, hd, St, list(b = cp$b), seed = seed_of(k, ti, 4L), overgeneration = M, background = TRUE, workers = 1L)
  te <- proc.time()[["elapsed"]]
  nb <- hd$grid$n_row * hd$grid$n_col
  cvec <- if (is.null(St)) rep(1, nb) else St$c; pk <- cvec / (if (is.null(St)) 1 else St$C)
  tr <- tx$transcripts
  Ecap <- sum(tr$n * pk[tr$bin_id]); Vcap <- sum(tr$n * pk[tr$bin_id] * (1 - pk[tr$bin_id]))
  cl <- cap$truth$cells; gl <- cap$truth$genes
  Ebg <- cp$b * sum(cvec)
  tot <- as.integer(round(Matrix::colSums(cap$counts_2um)))
  ob <- numeric(nb); tb <- tapply(cap$truth$origin$n, cap$truth$origin$bin_id, sum); ob[as.integer(names(tb))] <- tb
  bgv <- numeric(nb); bgv[cap$truth$background$bin_id] <- cap$truth$background$n
  modes <- NULL
  chk <- list(tag = tag, archetype = A, k = k, variant = VAR, seeds = SEEDS[, ti, k + 1L], M = M, C = if (is.null(St)) 1 else St$C,
    n_cells = nrow(L$cell_df), n_generated = sum(cl$n_generated), n_lost = sum(cl$n_lost), n_thinned = sum(cl$n_thinned), n_captured = sum(cl$n_captured),
    n_background = sum(bgv), n_target_R4 = sum(Matrix::colSums(prof$counts)) / M,
    profiles_match = isTRUE(all(cl$n_generated == as.numeric(Matrix::colSums(prof$counts)))),
    conservation_cells = all(cl$n_generated == cl$n_captured + cl$n_thinned + cl$n_lost),
    conservation_genes = all(gl$n_generated == gl$n_captured + gl$n_thinned + gl$n_lost) && all(as.numeric(Matrix::rowSums(cap$counts_2um)) == gl$n_captured + gl$n_background),
    conservation_bins = all(tot == ob + bgv) && all(as.numeric(Matrix::colSums(cap$counts_8um)) == as.numeric(sum_to_8um(matrix(tot, 1), hd$grid))),
    n_na = sum(is.na(tr$bin_id)) + sum(is.na(cap$truth$origin$n)) + sum(is.na(tot)),
    E_captured = Ecap, z_captured = (sum(cl$n_captured) - Ecap) / sqrt(max(Vcap, 1e-12)),
    E_background = Ebg, z_background = if (Ebg > 0) (sum(bgv) - Ebg) / sqrt(Ebg) else NA,
    w2 = p$w2, ell2 = p$ell2, b = cp$b,
    grid = list(x0 = hd$grid$x0, y0 = hd$grid$y0, n_row = hd$grid$n_row, n_col = hd$grid$n_col, bin_um = hd$grid$bin_um),
    elapsed_stripes_s = tp - ts, elapsed_profiles_s = tm - tp, elapsed_place_s = tc - tm, elapsed_capture_s = te - tc,
    elapsed_total_s = proc.time()[["elapsed"]] - t0, peak_mb = peak_mb())
  writeBin(tot, paste0(fo, ".bin"), size = 4L)
  cap$counts_2um <- NULL; saveRDS(cap, paste0(fo, "_cap.rds"))
  writeLines(jsonlite::toJSON(chk, auto_unbox = TRUE, digits = NA), paste0(fo, ".json.tmp")); file.rename(paste0(fo, ".json.tmp"), paste0(fo, ".json"))
  cat(sprintf("[%s k%d %s] M %.2f generati %d catturati %d fondo %d persi %d, %.0f s, %.0f MB\n", tag, k, VAR, M, chk$n_generated, chk$n_captured, chk$n_background, chk$n_lost, chk$elapsed_total_s, chk$peak_mb))
  fo
}
todo <- which(!file.exists(file.path(OUT, VAR, sprintf("%s_k%d.json", tags[jobs$t], jobs$k))))
r <- parallel::mclapply(todo, function(q) tryCatch(one(q), error = function(e) paste("ERRORE", tags[jobs$t[q]], conditionMessage(e))), mc.cores = W, mc.preschedule = FALSE)
print(unlist(r)[grepl("^ERRORE", unlist(r))])
