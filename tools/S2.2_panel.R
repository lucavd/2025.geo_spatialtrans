# S2.2, cancello G-2: simulate_cell_profiles() sulle 30 aree dello Step 1 (S1.6, predefinito dei preset), 2 semi, modello del vincitore
# della gara e modello M (controllo dei profili ricostruiti). Uso (root del repo):
#   Rscript --vanilla tools/S2.2_panel.R check [outdir=/mnt/micron/geo_spatialtrans/S2.2/panel] [workers=4]
# Scrive <outdir>/<A>_<roi>_s<seme>_<modello>.rds e results/S2.2/S2.2_G2_panel.csv (o <outdir>/S2.2_G2_panel.csv se outdir non e' quello predefinito).
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
a <- commandArgs(TRUE); OUT <- if (length(a) >= 2) a[2] else "/mnt/micron/geo_spatialtrans/S2.2/panel"; W <- if (length(a) >= 3) as.integer(a[3]) else 4L
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SEEDS <- c(20261108L, 20261109L)
lay <- sort(list.files("/mnt/micron/geo_spatialtrans/S1.6/layers", pattern = "_U6\\.rds$", full.names = TRUE))
rows <- list(); pooled <- list()
for (p in lay) {
  tag <- sub("_U6\\.rds$", "", basename(p)); A <- strsplit(tag, "_")[[1]][1]
  prof <- load_cell_profiles(A); cd <- readRDS(p)$cell_df
  for (sd in SEEDS) for (nz in unique(c(prof$noise, "M"))) {
    t0 <- Sys.time()
    s <- simulate_cell_profiles(cd[, c("cell_id", "cell_type")], prof, seed = sd + match(tag, sub("_U6\\.rds$", "", basename(lay))), noise = nz, workers = W)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    saveRDS(s, file.path(OUT, sprintf("%s_s%d_%s.rds", tag, sd, nz)))
    ok_sum <- all(Matrix::colSums(s$counts) == s$cells$n_transcripts)
    ok_type <- identical(s$cells$cell_type, as.character(cd$cell_type)) && all(startsWith(s$cells$subtype, paste0(s$cells$cell_type, "|")))
    rows[[length(rows) + 1]] <- data.frame(archetype = A, roi_id = sub("^[^_]*_", "", tag), seed = sd, noise = nz, winner_noise = prof$noise,
      n_cells = nrow(cd), n_na = sum(is.na(s$counts@x)) + sum(is.na(s$cells$n_transcripts)), min_total = min(s$cells$n_transcripts),
      sum_equals_total = ok_sum, type_and_subtype_ok = ok_type, seconds = el)
    key <- paste(A, sd, nz)
    agg <- Matrix::sparseMatrix(i = seq_len(nrow(s$cells)), j = match(s$cells$subtype, prof$subtypes$subtype), x = 1, dims = c(nrow(s$cells), nrow(prof$subtypes)))
    pooled[[key]] <- list(A = A, seed = sd, noise = nz, prof = prof,
      sums = if (is.null(pooled[[key]])) as.matrix(s$counts %*% agg) else pooled[[key]]$sums + as.matrix(s$counts %*% agg),
      n_sub = (if (is.null(pooled[[key]])) 0 else pooled[[key]]$n_sub) + Matrix::colSums(agg),
      n_type = c(if (!is.null(pooled[[key]])) pooled[[key]]$n_type, table(s$cells$cell_type)))
    cat(sprintf("[%s s%d %s] %d cellule, %.1f s\n", tag, sd, nz, nrow(cd), el))
  }
}
tab <- do.call(rbind, rows)
# quote dei sottotipi (tutti i modelli) e profili ricostruiti (solo M), per tessuto e seme, sottotipi con >= 200 cellule
rec <- list()
for (key in names(pooled)) { x <- pooled[[key]]; pr <- x$prof
  nt <- tapply(x$n_type, names(x$n_type), sum)
  zs <- vapply(seq_len(nrow(pr$subtypes)), function(i) { t <- pr$subtypes$cell_type[i]; p <- pr$subtypes$share[i]; n <- nt[[t]]
    if (p <= 0 || p >= 1) return(if (p >= 1) (x$n_sub[i] - n) else x$n_sub[i]); (x$n_sub[i] - n * p) / sqrt(n * p * (1 - p)) }, numeric(1))
  zp <- NA_real_; n_tested <- 0L
  if (x$noise == "M") { zz <- c()
    for (i in which(x$n_sub >= 200)) { Tt <- sum(x$sums[, i]); ex <- Tt * pr$P[i, ]; k <- ex >= 20
      zz <- c(zz, (x$sums[k, i] - ex[k]) / sqrt(ex[k] * (1 - pr$P[i, k]))); n_tested <- n_tested + 1L }
    zp <- if (length(zz)) max(abs(zz)) else NA_real_ }
  rec[[key]] <- data.frame(archetype = x$A, seed = x$seed, noise = x$noise, max_abs_z_share = max(abs(zs)), max_abs_z_profile = zp,
                           n_subtypes_tested = n_tested, n_genes_tested = if (x$noise == "M") NA else NA)
}
rec <- do.call(rbind, rec)
fo <- if (OUT == "/mnt/micron/geo_spatialtrans/S2.2/panel") "results/S2.2" else OUT
write.csv(tab, file.path(fo, "S2.2_G2_panel.csv"), row.names = FALSE)
write.csv(rec, file.path(fo, "S2.2_G2_recovery.csv"), row.names = FALSE)
print(rec)
cat(sprintf("G-2: %d esecuzioni; somma = totale %d/%d; tipo e sottotipo %d/%d; NA %d; quote |z| max %.2f; profili (M) |z| max %.2f\n",
            nrow(tab), sum(tab$sum_equals_total), nrow(tab), sum(tab$type_and_subtype_ok), nrow(tab), sum(tab$n_na),
            max(rec$max_abs_z_share), max(rec$max_abs_z_profile, na.rm = TRUE)))
