# S1.1 — costruzione degli input del pannello per extract_regions()
# Uso (dalla root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla tools/S1.1_inputs.R [gruppo ...]
# gruppi: adv syn real perf (default: tutti). Rilanciabile: salta gli input gia' presenti.
# Output: /mnt/micron/geo_spatialtrans/S1.1/inputs/<id>.rds (fuori git)
#         results/S1.1/S1.1_inputs.csv (tracciato: id, gruppo, n_px, md5 dell'RDS, sorgenti)
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(dplyr) })
source("R/03b_generate_synthetic_tissue.R")
source("R/04_clustering.R")

OUT  <- "/mnt/micron/geo_spatialtrans/S1.1/inputs"
DATA <- "data/real/datasets"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
dir.create("results/S1.1", recursive = TRUE, showWarnings = FALSE)
args <- commandArgs(trailingOnly = TRUE)
groups <- if (length(args)) args else c("adv", "syn", "real", "perf")

save_input <- function(id, group, clust, pixel_size_um, cluster_col = "intensity_cluster",
                       meta = list(), truth = NULL) {
  f <- file.path(OUT, paste0(id, ".rds"))
  obj <- list(id = id, group = group, clust = clust, pixel_size_um = pixel_size_um,
              cluster_col = cluster_col, meta = meta, truth = truth)
  saveRDS(obj, f)
  cat(sprintf("  %-34s %9d px  px=%g\n", id, nrow(clust), pixel_size_um))
}
done <- function(id) file.exists(file.path(OUT, paste0(id, ".rds")))
from_matrix <- function(m) {          # m[riga = y, colonna = x], 0 = fondo
  idx <- which(m != 0, arr.ind = TRUE)
  data.frame(x = idx[, 2], y = idx[, 1], intensity_cluster = as.integer(m[idx]))
}
threshold_df <- function(syn, thr = 0.7) {   # convenzione full_test.R righe 189-191
  d <- syn$img_df[syn$img_df$intensity / 255 < thr, ]
  d$value <- d$intensity / 255
  d[, c("x", "y", "value")]
}

# ---- avversari ----------------------------------------------------------------
if ("adv" %in% groups) {
  cat("[adv]\n")
  m <- outer(1:50, 1:50, function(i, j) ((i + j) %% 2) + 1)
  save_input("adv_checker", "adv", from_matrix(m), 1, meta = list(expect_components = 2500))
  m <- matrix(0, 12, 12); m[2:11, 2:11] <- 1; m[5:8, 5:8] <- 0
  save_input("adv_donut", "adv", from_matrix(m), 1, meta = list(expect_area = 100 - 16, expect_holes = 1))
  m[6:7, 6:7] <- 2
  save_input("adv_donut_island", "adv", from_matrix(m), 1, meta = list(expect_regions = 2))
  m <- matrix(0, 6, 6); m[1:2, 1:2] <- 1; m[3:4, 3:4] <- 1
  save_input("adv_corner", "adv", from_matrix(m), 1, meta = list(expect_components = 2))
  m <- matrix(0, 6, 6); m[1, 1:5] <- 1; m[1:5, 1] <- 1; m[1:5, 5] <- 1; m[5, 1:3] <- 1; m[4, 4] <- 1
  save_input("adv_pinch_hole", "adv", from_matrix(m), 1, meta = list(expect_components = 1, expect_area = 16))
  m <- matrix(1, 80, 100)
  save_input("adv_single_cluster", "adv", from_matrix(m), 1, meta = list(expect_regions = 1, expect_area = 8000))
  m <- matrix(0, 20, 200); m[10, 1:200] <- 1; m[1:20, 100] <- 2
  save_input("adv_line", "adv", from_matrix(m), 1, meta = list(expect_components = 3))
  m <- matrix(1, 60, 90); m[, 31:60] <- 2; m[, 61:90] <- 3; m[21:40, 1:90] <- 4
  save_input("adv_border", "adv", from_matrix(m), 1, meta = list(expect_components = 7))
  m <- matrix(0, 40, 40); m[5:13, 5:13] <- 1; m[20:29, 20:29] <- 2
  save_input("adv_blobs_9_10", "adv", from_matrix(m), 1,
             meta = list(expect_regions_at_100 = 1, expect_excluded_at_100 = 1))
  save_input("adv_pixel35", "adv", data.frame(x = 3L, y = 5L, intensity_cluster = 1L), 1,
             meta = list(expect_bbox = c(2, 4, 3, 5)))
  save_input("adv_L", "adv", data.frame(x = c(1:5, 1L, 1L), y = c(rep(1L, 5), 2L, 3L), intensity_cluster = 1L), 1,
             meta = list(expect_bbox = c(0, 0, 5, 3)))
}

# ---- sintetici ----------------------------------------------------------------
if ("syn" %in% groups) {
  cat("[syn]\n")
  for (cx in 1:4) {
    id <- sprintf("I%d_syn600_c%d", cx, cx); if (done(id)) next
    syn <- generate_synthetic_tissue(600, 600, complexity = cx, seed = 42)
    cl  <- cluster_image(threshold_df(syn), k_cell_types = 4, random_seed = 42, spatial_weight = 0.6)
    save_input(id, "syn", cl, 1, meta = list(complexity = cx, w = 600, h = 600, k = 4))
  }
  id <- "I5_syn600_c2_labels"
  if (!done(id)) {
    syn <- generate_synthetic_tissue(600, 600, complexity = 2, seed = 42, return_labels = TRUE)
    lm  <- syn$label_matrix
    save_input(id, "syn", from_matrix(lm), 1, meta = list(complexity = 2, w = 600, h = 600),
               truth = as.data.frame(table(patch = as.vector(lm)), stringsAsFactors = FALSE))
  }
  id <- "I6_syn6800x6500_c2"
  if (!done(id)) {
    syn <- generate_synthetic_tissue(6800, 6500, complexity = 2, seed = 42)
    cl  <- cluster_image(threshold_df(syn), k_cell_types = 4, random_seed = 42, spatial_weight = 0.6)
    save_input(id, "syn", cl, 1, meta = list(complexity = 2, w = 6800, h = 6500, k = 4, expect_stride = 6))
    rm(syn, cl); gc()
  }
}

# ---- reali (graphclust 8 µm) + nulli ----------------------------------------
if ("real" %in% groups) {
  cat("[real]\n")
  ds_all <- c(A1 = "Visium_HD_Mouse_Small_Intestine", A2A3 = "Visium_HD_Human_Colon_Cancer",
              A4 = "Visium_HD_Human_Lymph_Node_FFPE", A5 = "Visium_HD_6p5mm_Mouse_Brain",
              A6 = "Visium_HD_6p5mm_Human_Heart")
  rois <- read.csv("results/R2/R2_rois_checked.csv", stringsAsFactors = FALSE)
  src_rows <- list()
  for (tag in names(ds_all)) {
    ds <- ds_all[[tag]]
    sq <- file.path(DATA, ds, "binned_outputs_x/binned_outputs/square_008um")
    f_cl  <- file.path(sq, "analysis/clustering/gene_expression_graphclust/clusters.csv")
    f_pos <- file.path(sq, "spatial/tissue_positions.parquet")
    src_rows[[ds]] <- data.frame(dataset = ds, file = c(f_cl, f_pos),
                                 md5 = unname(tools::md5sum(c(f_cl, f_pos))))
    cl  <- read.csv(f_cl, stringsAsFactors = FALSE)
    pos <- as.data.frame(arrow::read_parquet(f_pos))
    d <- merge(pos, cl, by.x = "barcode", by.y = "Barcode")
    d <- d[order(d$array_row, d$array_col), ]
    clust <- data.frame(x = d$array_col + 1L, y = d$array_row + 1L, cluster = d$Cluster)
    meta <- list(dataset = ds, archetype = tag, n_in_tissue = sum(pos$in_tissue == 1),
                 n_clustered = nrow(clust), n_clusters = length(unique(clust$cluster)))
    id <- paste0("real_full_", tag)
    if (!done(id)) save_input(id, "real_full", clust, 8, cluster_col = "cluster", meta = meta)
    id <- paste0("null_full_", tag)
    if (!done(id)) {
      set.seed(20261007)
      nul <- clust; nul$cluster <- sample(nul$cluster)
      save_input(id, "null_full", nul, 8, cluster_col = "cluster", meta = c(meta, list(null = "permutazione delle etichette, seed 20261007")))
    }
    for (k in which(rois$dataset == ds)) {
      r <- rois[k, ]
      id <- sprintf("roi_%s_%s", r$archetype, r$roi_id); if (done(id)) next
      sel <- d$pxl_col_in_fullres >= r$c0 & d$pxl_col_in_fullres <= r$c1 &
             d$pxl_row_in_fullres >= r$r0 & d$pxl_row_in_fullres <= r$r1
      save_input(id, "real_roi", clust[sel, ], 8, cluster_col = "cluster",
                 meta = list(dataset = ds, archetype = r$archetype, roi_id = r$roi_id, label = r$label))
    }
  }
  write.csv(do.call(rbind, src_rows), "results/S1.1/S1.1_real_sources_md5.csv", row.names = FALSE)
}

# ---- prestazioni (griglie piene 4000x4000) ------------------------------------
if ("perf" %in% groups) {
  cat("[perf]\n")
  id <- "P1_4000_c2_labels"
  if (!done(id)) {
    syn <- generate_synthetic_tissue(4000, 4000, complexity = 2, seed = 42, return_labels = TRUE)
    save_input(id, "perf", from_matrix(syn$label_matrix), 1, meta = list(kind = "patch Voronoi Manhattan"))
    rm(syn); gc()
  }
  id <- "P2_4000_c4_quant4"
  if (!done(id)) {
    syn <- generate_synthetic_tissue(4000, 4000, complexity = 4, seed = 42)
    im <- syn$img_matrix
    q <- quantile(im, c(0.25, 0.5, 0.75), names = FALSE)
    lab <- matrix(findInterval(im, q) + 1L, nrow = nrow(im))
    save_input(id, "perf", from_matrix(lab), 1, meta = list(kind = "complessita' 4 quantizzata in 4 classi (quartili)"))
    rm(syn, im, lab); gc()
  }
}

# ---- indice -------------------------------------------------------------------
fs <- list.files(OUT, pattern = "\\.rds$", full.names = TRUE)
idx <- do.call(rbind, lapply(fs, function(f) {
  o <- readRDS(f)
  data.frame(id = o$id, group = o$group, n_px = nrow(o$clust), pixel_size_um = o$pixel_size_um,
             cluster_col = o$cluster_col, md5_rds = unname(tools::md5sum(f)))
}))
write.csv(idx, "results/S1.1/S1.1_inputs.csv", row.names = FALSE)
cat(sprintf("[indice] %d input\n", nrow(idx)))
