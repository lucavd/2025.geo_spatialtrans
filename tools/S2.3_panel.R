# tools/S2.3_panel.R - S2.3: place_and_move_transcripts() sulle 30 aree dello Step 1 (pannello) e cancello G-2.
# Definizioni: docs/steps/S2.3.md («La controprova sulle aree simulate», G-2). Uso (root del repo):
#   Rscript --vanilla tools/S2.3_panel.R <variante> [outdir] [workers]
#   varianti: main (2 semi, campo di R-1), field_off (1 seme; G-2 quota nel nucleo; CP-1), narrow (1 seme, w = 0, D -> 0; CP-2), r1 (1 seme, r = 1; CP-3)
# Scrive <outdir>/<variante>/<A>_<roi>_k<k>.bin (totali per quadratino, int32, ordine di bin_id, fondo di Poisson incluso), .json (griglia, controlli),
# e <outdir>/<variante>/G2.csv. Griglie di S2.1 in <outdir>/hd/<A>_<roi>.rds (deterministiche).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
a <- commandArgs(TRUE); VAR <- a[1]; OUT <- if (length(a) >= 2) a[2] else "/mnt/micron/geo_spatialtrans/S2.3/panel"; W <- if (length(a) >= 3) as.integer(a[3]) else 4L
stopifnot(VAR %in% c("main", "field_off", "narrow", "r1"))
dir.create(file.path(OUT, VAR), recursive = TRUE, showWarnings = FALSE); dir.create(file.path(OUT, "hd"), showWarnings = FALSE)
lay <- sort(list.files("/mnt/micron/geo_spatialtrans/S1.6/layers", pattern = "_U6\\.rds$", full.names = TRUE))
tags <- sub("_U6\\.rds$", "", basename(lay))
R1 <- if (file.exists("results/S2.3/S2.3_R1_choice.csv")) read.csv("results/S2.3/S2.3_R1_choice.csv") else NULL
SEEDS_PROF <- c(20261108L, 20261109L)
ks <- if (VAR == "main") 0:1 else 0L
jobs <- expand.grid(t = seq_along(tags), k = ks)
one <- function(q) {
  ti <- jobs$t[q]; k <- jobs$k[q]; tag <- tags[ti]; A <- strsplit(tag, "_")[[1]][1]
  fo <- file.path(OUT, VAR, sprintf("%s_k%d", tag, k))
  if (file.exists(paste0(fo, ".json"))) return(fo)
  t0 <- proc.time()[["elapsed"]]
  L <- readRDS(lay[ti])
  fh <- file.path(OUT, "hd", paste0(tag, ".rds"))
  hd <- if (file.exists(fh)) readRDS(fh) else { h <- build_hd_grid(L); saveRDS(h, paste0(fh, ".tmp")); file.rename(paste0(fh, ".tmp"), fh); h }
  win <- load_cell_profiles(A)$noise
  prof <- readRDS(sprintf("/mnt/micron/geo_spatialtrans/S2.2/panel/%s_s%d_%s.rds", tag, SEEDS_PROF[k + 1L], win))
  p <- transcript_params(A)
  field <- VAR %in% c("main", "narrow", "r1")
  if (field) { if (is.null(R1)) stop("manca results/S2.3/S2.3_R1_choice.csv (regola R-1)"); rr <- R1[R1$archetype == A, ]; p$field_s2 <- rr$s2; p$field_L <- rr$L }
  if (VAR == "narrow") { p$w <- 0; p$D <- 1e-9 }
  if (VAR == "r1") p$r <- 1
  seed <- 20261112L + 1000L * k + ti
  tx <- place_and_move_transcripts(L, hd, prof, p, seed = seed, field = field, frame_open = TRUE, workers = 1L)
  nb <- hd$grid$n_row * hd$grid$n_col
  tot <- tabulate(rep(tx$origin$bin_id, tx$origin$n), nb)
  RNGkind("L'Ecuyer-CMRG"); set.seed(seed + 500000L)
  bg <- stats::rpois(nb, p$b)
  Nc <- tx$cells$n_transcripts
  per_cell <- tapply(tx$origin$n, factor(tx$origin$cell_id, levels = L$cell_df$cell_id), sum); per_cell[is.na(per_cell)] <- 0
  nuc_born <- tapply(tx$origin$n_born_nucleus, factor(tx$origin$cell_id, levels = L$cell_df$cell_id), sum); nuc_born[is.na(nuc_born)] <- 0
  nuc_born <- as.numeric(nuc_born) + tx$lost$n_lost_born_nucleus
  chk <- list(tag = tag, archetype = A, k = k, variant = VAR, seed = seed, n_cells = nrow(L$cell_df), n_transcripts = sum(Nc),
    n_lost = sum(tx$lost$n_lost), conservation_cells = all(as.numeric(per_cell) + tx$lost$n_lost == Nc),
    conservation_total = sum(tx$origin$n) + sum(tx$lost$n_lost) == sum(Nc), conservation_genes = sum(tx$transcripts$n) == sum(tx$origin$n) &&
      all(as.numeric(tapply(tx$transcripts$n, factor(tx$transcripts$gene, levels = seq_len(nrow(prof$counts))), sum, default = 0)) + tx$lost_by_gene == as.numeric(Matrix::rowSums(prof$counts))),
    n_na = sum(is.na(tx$transcripts$bin_id)) + sum(is.na(tx$origin$n)) + sum(is.na(tx$lost$n_lost)),
    bins_in_grid = all(tx$origin$bin_id >= 1 & tx$origin$bin_id <= nb),
    nuc_born = sum(nuc_born), nuc_expected_no_field = sum(Nc * tx$cells$p_nucleus_no_field),
    nuc_var_no_field = sum(Nc * tx$cells$p_nucleus_no_field * (1 - tx$cells$p_nucleus_no_field)),
    grid = list(x0 = hd$grid$x0, y0 = hd$grid$y0, n_row = hd$grid$n_row, n_col = hd$grid$n_col, bin_um = hd$grid$bin_um),
    b = p$b, field_s2 = if (field) p$field_s2 else 0, field_L = if (field) p$field_L else NA, r = p$r, w = p$w, D = p$D,
    elapsed_function_s = tx$metadata$elapsed_s, elapsed_total_s = proc.time()[["elapsed"]] - t0)
  writeBin(as.integer(tot + bg), paste0(fo, ".bin"), size = 4L)
  writeLines(jsonlite::toJSON(chk, auto_unbox = TRUE, digits = NA), paste0(fo, ".json.tmp")); file.rename(paste0(fo, ".json.tmp"), paste0(fo, ".json"))
  cat(sprintf("[%s k%d %s] %d trascritti, persi %d, %.0f s\n", tag, k, VAR, sum(Nc), sum(tx$lost$n_lost), chk$elapsed_total_s))
  fo
}
todo <- which(!file.exists(file.path(OUT, VAR, sprintf("%s_k%d.json", tags[jobs$t], jobs$k))))
r <- parallel::mclapply(todo, function(q) tryCatch(one(q), error = function(e) paste("ERRORE", tags[jobs$t[q]], conditionMessage(e))), mc.cores = W, mc.preschedule = FALSE)
print(unlist(r)[grepl("^ERRORE", unlist(r))])
J <- lapply(list.files(file.path(OUT, VAR), pattern = "\\.json$", full.names = TRUE), function(f) jsonlite::fromJSON(f))
G2 <- do.call(rbind, lapply(J, function(j) data.frame(tag = j$tag, archetype = j$archetype, k = j$k, n_transcripts = j$n_transcripts, n_lost = j$n_lost,
  conservation = j$conservation_cells && j$conservation_total && j$conservation_genes, n_na = j$n_na, bins_in_grid = j$bins_in_grid,
  nuc_born = j$nuc_born, nuc_expected_no_field = j$nuc_expected_no_field, nuc_var_no_field = j$nuc_var_no_field,
  elapsed_function_s = j$elapsed_function_s)))
write.csv(G2, file.path(OUT, VAR, "G2.csv"), row.names = FALSE)
cat(sprintf("G-2 (%s): %d esecuzioni; conservazione %d; NA %d; quadratini nella griglia %d\n", VAR, nrow(G2), sum(G2$conservation), sum(G2$n_na), sum(G2$bins_in_grid)))
if (VAR == "field_off") {
  z <- tapply(seq_len(nrow(G2)), G2$archetype, function(ix) (sum(G2$nuc_born[ix]) - sum(G2$nuc_expected_no_field[ix])) / sqrt(sum(G2$nuc_var_no_field[ix])))
  print(z); cat(sprintf("G-2 quota nel nucleo: max |z| %.2f, soglia di Bonferroni %.2f\n", max(abs(z)), qnorm(1 - 0.0005 / length(z))))
}
