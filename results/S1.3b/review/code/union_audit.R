
# review/code/union_audit.R — 15 nulli con unione V0 sbagliata: basta fondere i vertici non identici (senza giunzioni a T)?
REAL_ARGS <- commandArgs(TRUE)
RV <- "/mnt/micron/geo_spatialtrans/S1.3b/review_code"; dir.create(RV, showWarnings = FALSE)
args <- c("none", file.path(RV, "s13out")); commandArgs <- function(trailingOnly = TRUE) args
suppressMessages(source("tools/S1.3_roi.R")); rm(commandArgs); source("tools/S1.3b_overlay.R")
src <- readLines("tools/S1.3b_run.R"); eval(parse(text = src[grep("^regen <- function", src):(grep("^reg_features <- function", src) - 1)]))
suppressPackageStartupMessages(library(parallel))
roi_i <- function(A, roi) which(rois$archetype == A & rois$roi_id == roi)
T <- read.csv("results/S1.3b/review/code/null_long.csv"); t0 <- T[T$v == "V0", ]
F <- t0[abs(t0$ov)/t0$A > 1e-9 | t0$symd/t0$A > 1e-9, c("arch", "roi", "model", "rep")]
uerr <- function(terr, reg, A) { u <- st_union(terr); a <- sum(as.numeric(st_area(terr)))
  c(ov_rel = (a - as.numeric(st_area(u))) / A, symd_rel = sum(as.numeric(st_area(st_sym_difference(u, reg)))) / A) }
one <- function(q) { cs <- F[q, ]; G <- regen(roi_i(cs$arch, cs$roi), cs$model, cs$rep); A <- G$w$area_um2
  e0 <- uerr(G$terr, G$reg, A)
  geoms <- c(list(G$reg[[1]]), s3b_sfg_list(G$terr))
  nf <- s3b_node(geoms, 1e-6, tjunction = FALSE); e_fuse <- uerr(st_sfc(nf$geoms[-1]), st_sfc(nf$geoms[1]), A)
  nt <- s3b_node(geoms, 1e-7, tjunction = TRUE); e_v3 <- uerr(st_sfc(nt$geoms[-1]), st_sfc(nt$geoms[1]), A)
  # posizione dei vertici non identici: sul bordo della regione?
  v <- s3b_geoms_df(geoms); fu <- s3b_fuse(v$x, v$y, 1e-6); cr <- v$g[fu$I] != v$g[fu$J] & fu$d > 0
  onreg <- if (any(cr)) as.numeric(st_distance(st_sfc(lapply(fu$I[cr], function(k) st_point(c(v$x[k], v$y[k])))), st_boundary(G$reg))) <= 1e-6 else logical(0)
  involves_reg <- if (any(cr)) sum(v$g[fu$I[cr]] == 1 | v$g[fu$J[cr]] == 1) else 0
  data.frame(cs, A = A, ov0 = e0[[1]], sd0 = e0[[2]], ov_fuse = e_fuse[[1]], sd_fuse = e_fuse[[2]], ov_v3 = e_v3[[1]], sd_v3 = e_v3[[2]],
             n_diff = sum(cr), n_diff_onreg = sum(onreg), n_diff_with_region_vertex = involves_reg, n_tj_v3 = nt$info[["n_tjunction"]], d_max = if (any(cr)) max(fu$d[cr]) else 0)
}
res <- mclapply(seq_len(nrow(F)), function(q) tryCatch(one(q), error = function(e) data.frame(F[q, ], err = conditionMessage(e))), mc.cores = 15, mc.preschedule = FALSE)
R <- do.call(rbind, lapply(res, function(r) { if (!"err" %in% names(r)) r$err <- NA; r }))
write.csv(R, "results/S1.3b/review/code/union_audit.csv", row.names = FALSE)
R$fixed_fuse <- abs(R$ov_fuse) <= 1e-9 & R$sd_fuse <= 1e-9; R$fixed_v3 <- abs(R$ov_v3) <= 1e-9 & R$sd_v3 <= 1e-9
print(R, row.names = FALSE)
cat(sprintf("fusion-only fixes %d/15; V3(1e-7) fixes %d/15\n", sum(R$fixed_fuse, na.rm = TRUE), sum(R$fixed_v3, na.rm = TRUE)))
