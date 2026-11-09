# tools/S2.5_truth8_check.R - S2.5, G-2 rafforzato dopo la revisione G-5: verita' a 8 um per coppia (blocco, cellula) su tutte le esecuzioni del
# pannello, con l'indice del blocco della tabella dei quadratini di S2.1 (hd$bins$bin8_id), indipendente dall'aritmetica dell'orchestratore.
# Uso (root del repo): Rscript --vanilla tools/S2.5_truth8_check.R <panel_dir> <out.csv> [workers]
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix) })
a <- commandArgs(TRUE); PD <- a[1]; OUT <- a[2]; W <- if (length(a) >= 3) as.integer(a[3]) else 4L
HD23 <- "/mnt/micron/geo_spatialtrans/S2.3/panel/hd"
fs <- unlist(lapply(c("main", "depth_off", "main_r2"), function(v) list.files(file.path(PD, v), pattern = "_k[01]_out\\.rds$", full.names = TRUE)))
pairs <- function(b8, cell, M) { u <- sort(unique(cell)); key <- (as.numeric(b8) - 1) * (length(u) + 1) + match(cell, u)
  s <- rowsum(M, key); list(key = as.numeric(rownames(s)), s = s, u = u) }
one <- function(f) {
  o <- readRDS(f); tag <- sub("_k[01]_out\\.rds$", "", basename(f)); hd <- readRDS(file.path(HD23, paste0(tag, ".rds")))
  tr <- o$truth; b8of <- function(id) hd$bins$bin8_id[match(id, hd$bins$bin_id)]
  cmp <- function(df, df8, cols, tol) {
    p <- pairs(b8of(df$bin_id), df$cell_id, as.matrix(df[cols]))
    k8 <- (as.numeric(df8$bin8_id) - 1) * (length(p$u) + 1) + match(df8$cell_id, p$u)
    if (anyNA(k8) || length(k8) != length(p$key) || anyDuplicated(k8)) return(FALSE)
    m <- match(p$key, k8); !anyNA(m) && all(abs(as.matrix(df8[m, cols]) - p$s) <= tol)
  }
  ok_or <- cmp(tr$origin, tr$origin_8um, c("n", "n_born_nucleus"), 0)
  ok_ov <- cmp(hd$truth$overlap, tr$overlap_8um, c("area_territory_um2", "area_nucleus_um2"), 1e-9)
  ok_fr <- all(abs(tr$overlap_8um$frac_territory - tr$overlap_8um$area_territory_um2 / 64) < 1e-12) && all(abs(tr$overlap_8um$frac_nucleus - tr$overlap_8um$area_nucleus_um2 / 64) < 1e-12)
  bg <- rowsum(as.numeric(tr$background$n), b8of(tr$background$bin_id)); bg8 <- setNames(as.numeric(tr$background_8um$n), tr$background_8um$bin8_id)
  ok_bg <- length(bg8) == nrow(bg) && all(bg8[rownames(bg)] == bg[, 1])
  data.frame(file = basename(dirname(f)), run = basename(f), origin_pairs = nrow(tr$origin_8um), overlap_pairs = nrow(tr$overlap_8um),
             origin_ok = ok_or, overlap_ok = ok_ov, fractions_ok = ok_fr, background_ok = ok_bg)
}
R <- do.call(rbind, parallel::mclapply(fs, function(f) tryCatch(one(f), error = function(e) data.frame(file = f, run = conditionMessage(e), origin_pairs = NA, overlap_pairs = NA,
  origin_ok = FALSE, overlap_ok = FALSE, fractions_ok = FALSE, background_ok = FALSE)), mc.cores = W, mc.preschedule = FALSE))
write.csv(R, OUT, row.names = FALSE)
cat(sprintf("verita' a 8 um per coppia: %d esecuzioni; origine %d, sovrapposizione %d, frazioni %d, fondo %d esatte\n", nrow(R), sum(R$origin_ok), sum(R$overlap_ok), sum(R$fractions_ok), sum(R$background_ok)))
quit(status = if (all(R$origin_ok & R$overlap_ok & R$fractions_ok & R$background_ok)) 0L else 1L)
