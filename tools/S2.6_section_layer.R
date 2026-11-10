# tools/S2.6_section_layer.R - S2.6: livello cellulare di una sezione intera dalla maschera del tessuto di Space Ranger (quadratini di 8 um).
# Definizioni: docs/steps/S2.6.md. Pixel = quadratino di 8 um (x = array_col + 1, y = array_row + 1, pixel_size_um = 8); un solo gruppo di
# regioni; preset del tessuto; argomenti di place_cells() come i livelli U6 di S1.6 (semina dei preset). Al livello si aggiungono quattro
# ancore (quadratini di 1e-4 um) agli angoli dell'area di cattura [0, 6704] um, cosi' la griglia di simulate_visium_hd() copre l'area intera.
# Uso: Rscript --vanilla tools/S2.6_section_layer.R <A> <seed> <outdir>
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
peak_mb <- function() { s <- readLines("/proc/self/status"); as.numeric(sub("[^0-9]*([0-9]+).*", "\\1", s[grepl("^VmHWM", s)])) / 1024 }
PRESET <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
SECTION <- c(A1 = "Visium_HD_Mouse_Small_Intestine", A2 = "Visium_HD_Human_Colon_Cancer", A4 = "Visium_HD_Human_Lymph_Node_FFPE",
             A5 = "Visium_HD_6p5mm_Mouse_Brain", A6 = "Visium_HD_6p5mm_Human_Heart")
a <- commandArgs(TRUE); A <- a[1]; seed <- as.integer(a[2]); OUT <- a[3]
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
tag <- sprintf("sec_%s_s%d", A, seed); fo <- file.path(OUT, paste0(tag, ".rds"))
if (file.exists(fo)) quit(status = 0)
f <- file.path("/mnt/micron/geo_spatialtrans/data_real", SECTION[[A]], "binned_outputs_x/binned_outputs/square_008um/spatial/tissue_positions.parquet")
tp <- as.data.frame(arrow::read_parquet(f, col_select = c("in_tissue", "array_row", "array_col")))
W <- 8 * (max(tp$array_col) + 1); Hh <- 8 * (max(tp$array_row) + 1)
tp <- tp[tp$in_tissue == 1, ]
clust <- data.frame(x = tp$array_col + 1L, y = tp$array_row + 1L, intensity_cluster = 1L)
gc(); t0 <- proc.time()[["elapsed"]]; pf <- tempfile(); Rprof(pf, interval = 0.5)
L <- suppressWarnings(place_cells(clust, PRESET[[A]], pixel_size_um = 8, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
                  small_components = "exclude", random_seed = seed, types = "multi", size_arm = "recentred",
                  spatial = "uniform", size_resid = "gaussian", strict = FALSE, verbose = FALSE))
te <- proc.time()[["elapsed"]] - t0; Rprof(NULL)
s <- tryCatch(summaryRprof(pf)$by.total, error = function(e) NULL); rn <- gsub('"', "", rownames(s))
keys <- c("extract_regions", "seed_centroids", "tessellate_voronoi", "derive_cell_geometry", "validate_cell_layer")
prof <- setNames(lapply(keys, function(k) { w <- which(rn == k); if (length(w)) s$total.time[w[1]] else 0 }), keys)
an <- 1e-4
anc <- sf::st_sfc(lapply(list(c(0, 0), c(W - an, 0), c(0, Hh - an), c(W - an, Hh - an)),
                         function(p) sf::st_polygon(list(rbind(p, p + c(an, 0), p + c(an, an), p + c(0, an), p)))))
L$region_polygons <- c(sf::st_set_crs(L$region_polygons, NA), anc)
L$metadata$s26_capture_box <- c(0, 0, W, Hh)
saveRDS(L, fo)
chk <- list(tag = tag, archetype = A, section = SECTION[[A]], seed = seed, n_cells = nrow(L$cell_df), n_regions = nrow(L$region_df),
            tissue_mm2_layer = sum(L$region_df$area_um2) / 1e6, tissue_mm2_spaceranger = nrow(tp) * 64 / 1e6, capture_box = c(0, 0, W, Hh),
            elapsed_s = te, peak_mb = peak_mb(), rprof_total_s = prof, load = readLines("/proc/loadavg"))
writeLines(jsonlite::toJSON(chk, auto_unbox = TRUE, digits = NA), file.path(OUT, paste0(tag, ".json")))
cat(sprintf("[%s] %d cellule, %d regioni, %.1f mm2, %.0f s, %.0f MB\n", tag, chk$n_cells, chk$n_regions, chk$tissue_mm2_layer, te, chk$peak_mb))
