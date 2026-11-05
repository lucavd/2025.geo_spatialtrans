# tools/S2.1_cp.R - S2.1, controprova CP-1: la nostra griglia e hd_overlap() contro l'assegnazione dei quadratini ai nuclei di Space Ranger 4.
# Definizioni: docs/steps/S2.1.md (Controprova). Uso (root del repo):
#   Rscript --vanilla tools/S2.1_cp.R explore      E-1 su tumore del colon r1: corrispondenza degli id, orientamento, affine (non guarda CP-1)
#   Rscript --vanilla tools/S2.1_cp.R run <KEY>     KEY = A23 | A4 | A5 | A6 | A1 (A1: solo i controlli della griglia, niente nuclei)
#   Rscript --vanilla tools/S2.1_cp.R collect       tabelle in results/S2.1/
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos); library(arrow); library(Matrix) })
for (f in sort(list.files("R", pattern = "^(04b5|05b).*\\.R$", full.names = TRUE))) source(f)
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
args <- commandArgs(TRUE); STAGE <- args[1]
M <- "/mnt/micron/geo_spatialtrans"; D <- file.path(M, "data_real"); OUT <- file.path(M, Sys.getenv("S21_RUN", "S2.1/run1"), "cp")
RES <- Sys.getenv("S21_RES", "results/S2.1"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE); dir.create(RES, recursive = TRUE, showWarnings = FALSE)
DS <- list(A1 = list(ds = "Visium_HD_Mouse_Small_Intestine", A = "A1"), A23 = list(ds = "Visium_HD_Human_Colon_Cancer", A = c("A2", "A3")),
           A4 = list(ds = "Visium_HD_Human_Lymph_Node_FFPE", A = "A4"), A5 = list(ds = "Visium_HD_6p5mm_Mouse_Brain", A = "A5"),
           A6 = list(ds = "Visium_HD_6p5mm_Human_Heart", A = "A6"))
rois <- read.csv("results/R2/R2_rois_checked.csv")
sq_dir <- function(ds, s) file.path(D, ds, "binned_outputs_x/binned_outputs", s, "spatial")
parse_bc <- function(x) { p <- regmatches(x, regexec("_(\\d+)_(\\d+)-1$", x)); m <- do.call(rbind, p); cbind(as.integer(m[, 2]), as.integer(m[, 3])) }
# trasformazione affine pixel (colonna, riga) -> convenzione della griglia (2 j + 1, 2 i + 1) um, minimi quadrati su tutti i quadratini
fit_affine <- function(tp, flip = "none") {
  i <- tp$array_row; j <- tp$array_col
  tgt <- switch(flip, none = cbind(u = 2 * j + 1, v = 2 * i + 1), transpose = cbind(u = 2 * i + 1, v = 2 * j + 1),
                half = cbind(u = 2 * j + 2, v = 2 * i + 2), reverse = cbind(u = 2 * j + 1, v = 2 * (max(i) - i) + 1))
  X <- cbind(1, tp$pxl_col_in_fullres, tp$pxl_row_in_fullres)
  B <- solve(crossprod(X), crossprod(X, tgt)); res <- tgt - X %*% B
  list(B = B, max_res_um = max(abs(res)), flip = flip)
}
apply_affine <- function(xy, af) cbind(1, xy) %*% af$B
# poligoni del geojson nelle coordinate della griglia
to_grid <- function(geom, af) sf::st_sfc(lapply(geom, function(p) sf::st_polygon(lapply(p, function(r) apply_affine(r, af)))))
load_section <- function(ds) {
  t0 <- proc.time()[[3]]
  tp <- as.data.frame(arrow::read_parquet(file.path(sq_dir(ds, "square_002um"), "tissue_positions.parquet")))
  sf_json <- jsonlite::fromJSON(file.path(sq_dir(ds, "square_002um"), "scalefactors_json.json"))
  list(tp = tp, mpp = sf_json$microns_per_pixel, t_load = proc.time()[[3]] - t0)
}
load_sr <- function(ds) {
  nuc <- sf::st_read(file.path(D, ds, "segmented_outputs_x/segmented_outputs/nucleus_segmentations.geojson"), quiet = TRUE)
  sf::st_crs(nuc) <- NA                                   # coordinate in pixel dell'immagine, non geografiche
  bm <- as.data.frame(arrow::read_parquet(file.path(D, ds, paste0(ds, "_barcode_mappings.parquet")), col_select = c("square_002um", "square_008um", "cell_id", "in_nucleus")))
  list(nuc = nuc, bm = bm)
}
grid_for <- function(u, v) list(x0 = 2 * floor(min(u) / 2) - 8, y0 = 2 * floor(min(v) / 2) - 8, bin_um = 2,
                                n_row = as.integer(4 * ceiling((max(v) - 2 * floor(min(v) / 2) + 16) / 8)),
                                n_col = as.integer(4 * ceiling((max(u) - 2 * floor(min(u) / 2) + 16) / 8)), block = 4L, margin_um = 0)
# per un'area: coppie (quadratino SR, nucleo) nostre e di Space Ranger
roi_cp <- function(R, sec, sr, af) {
  bb <- sf::st_bbox(sr$nuc)
  box <- sf::st_as_sfc(sf::st_bbox(c(xmin = R$c0, ymin = R$r0, xmax = R$c1, ymax = R$r1)))
  k <- unlist(sf::st_intersects(box, sr$nuc)); nuc <- sr$nuc[k, ]
  G <- to_grid(sf::st_geometry(nuc), af); xy <- sf::st_coordinates(G)
  grid <- grid_for(xy[, 1], xy[, 2])
  ov <- hd_overlap(grid, G, NULL, tile_um = 500)
  ov$cell <- as.integer(nuc$cell_id[ov$cell_id])
  bi <- (ov$bin_id - 1L) %/% grid$n_col; bj <- (ov$bin_id - 1L) %% grid$n_col
  # indice del quadratino nella convenzione: i = (y - y0) / 2 -> riga, j -> colonna (riga e colonna di Space Ranger se la griglia e' giusta)
  ov$i <- as.integer(bi + round(grid$y0 / 2)); ov$j <- as.integer(bj + round(grid$x0 / 2)); ov$frac <- ov$area_territory_um2 / 4
  cid <- sprintf("cellid_%09d-1", nuc$cell_id)
  b <- sr$bm[sr$bm$in_nucleus & sr$bm$cell_id %in% cid, ]
  ij <- parse_bc(b$square_002um)
  srp <- data.frame(i = ij[, 1], j = ij[, 2], cell = as.integer(sub("cellid_0*(\\d+)-1", "\\1", b$cell_id)))
  list(ov = ov, sr = srp, n_nuc = nrow(nuc))
}
cp_stats <- function(x) {
  k <- x$ov$frac >= 0.5; ours <- unique(sprintf("%d_%d_%d", x$ov$i[k], x$ov$j[k], x$ov$cell[k]))
  srs <- unique(sprintf("%d_%d_%d", x$sr$i, x$sr$j, x$sr$cell))
  c(n_nuclei = x$n_nuc, n_pairs = nrow(x$ov), frac_mid = mean(x$ov$frac > 0.1 & x$ov$frac < 0.9),
    n_ours = length(ours), n_sr = length(srs), jaccard = length(intersect(ours, srs)) / length(union(ours, srs)))
}
if (STAGE == "explore") {
  sec <- load_section(DS$A23$ds); sr <- load_sr(DS$A23$ds)
  cat("E-1 colonne tissue_positions:", paste(names(sec$tp), collapse = ", "), "\n")
  af <- fit_affine(sec$tp); cat(sprintf("E-1 affine: max scarto %.4f um; coefficienti\n", af$max_res_um)); print(af$B)
  A <- af$B[2:3, ]; cat(sprintf("E-1 passo della griglia in pixel per quadratino: %.4f; con microns_per_pixel %.6f -> %.5f um; determinante %.3e\n",
                                2 / sqrt(abs(det(A))), sec$mpp, 2 / sqrt(abs(det(A))) * sec$mpp, det(A)))
  cid <- sprintf("cellid_%09d-1", sr$nuc$cell_id[1:2000]); cat(sprintf("E-1 id del geojson trovati in barcode_mappings (primi 2000): %d\n", sum(cid %in% sr$bm$cell_id)))
  cat("E-1 esempi cell_id di barcode_mappings:", head(unique(sr$bm$cell_id[!is.na(sr$bm$cell_id)]), 3), "\n")
  cat(sprintf("E-1 nuclei nel geojson %d; righe in_nucleus %d\n", nrow(sr$nuc), sum(sr$bm$in_nucleus)))
  R <- rois[rois$archetype == "A2" & rois$roi_id == "r1", ]
  cat(sprintf("E-1 tumore r1: riquadro pixel colonne %d-%d righe %d-%d; um_per_px %.5f\n", R$c0, R$c1, R$r0, R$r1, R$um_per_px))
  quit(status = 0)
}
if (STAGE %in% c("run", "gridonly")) {
  KEY <- args[2]; info <- DS[[KEY]]; t0 <- proc.time()[[3]]
  sec <- load_section(info$ds); tp <- sec$tp
  af <- fit_affine(tp); A <- af$B[2:3, ]; pitch_um <- 2 / sqrt(abs(det(A))) * sec$mpp
  # blocchi di 8 um: posizioni dei centri di 8 um = media dei 16 figli
  tp8 <- as.data.frame(arrow::read_parquet(file.path(sq_dir(info$ds, "square_008um"), "tissue_positions.parquet")))
  key8 <- (tp$array_row %/% 4L) * 1e5 + tp$array_col %/% 4L
  mc <- tapply(tp$pxl_col_in_fullres, key8, mean); mr <- tapply(tp$pxl_row_in_fullres, key8, mean); nch <- tapply(key8, key8, length)
  k8 <- match(tp8$array_row * 1e5 + tp8$array_col, as.numeric(names(mc)))
  grid_rows <- data.frame(key = KEY, dataset = info$ds, n_bins = nrow(tp), affine_max_res_um = af$max_res_um, mpp = sec$mpp, pitch_um = pitch_um,
                          n_8um = nrow(tp8), n_8um_matched = sum(!is.na(k8)), n_8um_16children = sum(nch[k8] == 16, na.rm = TRUE),
                          max_dev_8um_px = max(abs(c(mc[k8] - tp8$pxl_col_in_fullres, mr[k8] - tp8$pxl_row_in_fullres)), na.rm = TRUE))
  full <- !is.na(k8) & nch[k8] == 16; edge <- tp8$array_row == max(tp8$array_row) | tp8$array_col == max(tp8$array_col)
  grid_rows$max_dev_8um_px_16 <- max(abs(c(mc[k8][full] - tp8$pxl_col_in_fullres[full], mr[k8][full] - tp8$pxl_row_in_fullres[full])))
  grid_rows$n_2um_side <- max(tp$array_row) + 1L; grid_rows$n_partial_8um <- sum(!full); grid_rows$n_partial_on_last_row_col <- sum(!full & edge)
  if (KEY != "A1") {
    sr <- if (STAGE == "run") load_sr(info$ds) else list(bm = as.data.frame(arrow::read_parquet(file.path(D, info$ds, paste0(info$ds, "_barcode_mappings.parquet")), col_select = c("square_002um", "square_008um"))))
    ij <- parse_bc(sr$bm$square_002um); ij8 <- parse_bc(sr$bm$square_008um)
    grid_rows$n_bm <- nrow(sr$bm); grid_rows$n_bm_block_mismatch <- sum(ij8[, 1] != ij[, 1] %/% 4L | ij8[, 2] != ij[, 2] %/% 4L)
    out <- list()
    if (STAGE == "run") for (A_ in info$A) for (r in which(rois$archetype == A_)) {
      R <- rois[r, ]; pad <- 200
      tpr <- tp[tp$pxl_col_in_fullres >= R$c0 - pad & tp$pxl_col_in_fullres < R$c1 + pad & tp$pxl_row_in_fullres >= R$r0 - pad & tp$pxl_row_in_fullres < R$r1 + pad, ]
      for (fl in c("none", "transpose", "half", "reverse")) {
        afr <- fit_affine(tpr, fl)                        # trasformazione per area (correzione 2 prima del calcolo)
        x <- roi_cp(R, sec, sr, afr); s <- cp_stats(x)
        out[[length(out) + 1]] <- data.frame(archetype = A_, roi_id = R$roi_id, grid = fl, affine_max_res_um = afr$max_res_um,
                                             affine_max_res_px = afr$max_res_um / sec$mpp, n_bins_fit = nrow(tpr), t(s))
        if (fl == "none") saveRDS(x, file.path(OUT, sprintf("%s_%s_cp.rds", A_, R$roi_id)))
        cat(sprintf("[%s] %s %s %s: scarto %.3f px, jaccard %.5f, a meta' %.5f\n", KEY, A_, R$roi_id, fl, afr$max_res_um / sec$mpp, s["jaccard"], s["frac_mid"]), flush = TRUE)
      }
    }
    if (STAGE == "run") write.csv(do.call(rbind, out), file.path(OUT, sprintf("cp_%s.csv", KEY)), row.names = FALSE)
  }
  grid_rows$elapsed_s <- proc.time()[[3]] - t0
  write.csv(grid_rows, file.path(OUT, sprintf("grid_%s.csv", KEY)), row.names = FALSE)
  print(grid_rows)
}
if (STAGE == "collect") {
  cpt <- do.call(rbind, lapply(list.files(OUT, "^cp_.*csv$", full.names = TRUE), read.csv))
  gl <- lapply(list.files(OUT, "^grid_.*csv$", full.names = TRUE), function(f) { x <- read.csv(f); x[setdiff(names(x), "elapsed_s")] })
  allc <- unique(unlist(lapply(gl, names)))
  grt <- do.call(rbind, lapply(gl, function(x) { for (m in setdiff(allc, names(x))) x[[m]] <- NA; x[allc] }))
  cpt <- cpt[order(cpt$grid, cpt$archetype, cpt$roi_id), ]
  write.csv(cpt, file.path(RES, "S2.1_CP1_roi.csv"), row.names = FALSE); write.csv(grt, file.path(RES, "S2.1_CP1_grid.csv"), row.names = FALSE)
  ok <- cpt$grid == "none"
  a <- all(cpt$frac_mid[ok] <= 0.001); b <- all(cpt$jaccard[ok] >= 0.999)
  cc <- all(cpt$jaccard[!ok] < 0.9)                                   # correzione prima del calcolo (docs/steps/S2.1.md): CP-1a deve fallire solo
  wrong_fail_a <- tapply(cpt$frac_mid[!ok] > 0.001, cpt$grid[!ok], all) # con la griglia spostata di mezzo quadratino
  blk <- all(grt$n_bm_block_mismatch[!is.na(grt$n_bm_block_mismatch)] == 0) && all(grt$n_8um_matched == grt$n_8um) && all(grt$n_partial_8um == grt$n_partial_on_last_row_col) &&
    all(grt$max_dev_8um_px_16 < 0.5) && all(abs(grt$pitch_um / 2 - 1) < 0.005) && all(cpt$affine_max_res_px[ok] < 0.5)
  cat(sprintf("CP-1a %s: quota a meta' max %.6f (25 aree)\n", if (a) "PASS" else "FAIL", max(cpt$frac_mid[ok])))
  cat(sprintf("CP-1b %s: jaccard min %.6f (25 aree)\n", if (b) "PASS" else "FAIL", min(cpt$jaccard[ok])))
  cat(sprintf("CP-1c %s: jaccard max delle griglie sbagliate %.4f; CP-1a fallisce per griglia: %s\n", if (cc && isTRUE(wrong_fail_a[["half"]])) "PASS" else "FAIL",
              max(cpt$jaccard[!ok]), paste(names(wrong_fail_a), wrong_fail_a, collapse = " ")))
  cat(sprintf("CP-1 blocchi, passo e scarto per area %s; scarto dell'affine su tutta la sezione (giudicato come scritto, < 0.5 px): %s\n",
              if (blk) "PASS" else "FAIL", paste(sprintf("%s %.2f px", grt$key, grt$affine_max_res_um / grt$mpp), collapse = ", "))); print(grt)
}
