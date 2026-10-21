# tools/S1.5_posthoc.R - S1.5, analisi post hoc (dopo i verdetti; non cambiano previsioni, intervalli ne' verdetti).
# PH-1 contatti e PH-2 rapporto nucleo/territorio: tre bracci sugli stessi ROI e semi, stesso preset (conversione "preserve_area"):
#   sim_manual = pannello (semina RSA d*, densita' manuale)          -> gia' in S1.5/gen/preserve_area e raster/preserve_area
#   sim_segdens = semina RSA d* alla densita' dei nuclei segmentati del ROI (stessa taglia assoluta del nucleo)
#   realgen = generatori e territori reali del ROI (S1.4/real, R3b/S1.4/real), nuclei del preset con la densita' del ROI
# Uso: Rscript --vanilla tools/S1.5_posthoc.R gen  (poi raster + tools/S1.5_windows.py con i bracci posthoc)
source("tools/S1.5_bench.R")
rr <- read.csv(file.path(RES, "S1.5_real_refs_roi.csv"))
ct_at <- function(A, dens) { sp <- tissue_preset_spec(PRESET[[A]]); s <- setNames(sp$sources$value, sp$sources$parameter)
  ct <- sp$cell_types; ct$density <- dens; ct$nucleus_to_eq_ratio <- s[["rho"]] * (dens / s[["density_model"]])^((1 - s[["beta"]]) / 2); ct }
if (STAGE == "none" && length(args) == 0) STAGE <- "gen"
ph_dir <- file.path(OUT, "posthoc"); for (d in c("sim_segdens", "realgen")) { dir.create(file.path(ph_dir, "units", d), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(OUT, "nuc", d), showWarnings = FALSE) }
jobs <- expand.grid(i = rois$roi_idx, s = SEEDS, arm = c("sim_segdens", "realgen"), stringsAsFactors = FALSE)
par_run(jobs, function(j) {
  i <- j$i; s <- j$s; A <- rois$archetype[i]; roi <- rois$roi_id[i]
  f <- file.path(ph_dir, "units", j$arm, sprintf("%s_%s_s%02d.rds", A, roi, s)); if (file.exists(f)) return(invisible())
  r <- rr[rr$archetype == A & rr$roi_id == roi, ]; dens <- r$density_seg_mm2; ct <- ct_at(A, dens)
  w <- r3_roi_window(A, roi, rois)
  if (j$arm == "sim_segdens") {
    L <- suppressWarnings(place_cells(clust_of(w), cell_types = ct, pixel_size_um = w$px, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
                                      small_components = "exclude", random_seed = seed_of(i, s), keep_tiles = TRUE, strict = FALSE, verbose = FALSE))
    cd <- L$cell_df; td <- L$territory_df
    in_frame <- geos::geos_within(geos::as_geos_geometry(L$tiles), geos::as_geos_geometry(w$frame))
    it <- !td$clipped | (in_frame & td$territory_area / td$tile_area > 0.999); nuc <- L$cell_nuclei
  } else {
    base <- if (A %in% c("A4", "A6")) file.path(M, "R3b/S1.4/real") else file.path(M, "S1.4/real")
    z <- readRDS(file.path(base, sprintf("%s_%s_%s.rds", A, roi, r$method)))
    cen <- data.frame(cell_id = seq_len(z$n), region_id = 1L, x = z$cells$x, y = z$cells$y, cell_type = ct$cell_type)
    tess <- list(cell_territories = z$terr, territory_df = data.frame(cell_id = seq_len(z$n), territory_area = z$cells$a_t))
    dg <- derive_cell_geometry(cen, tess, ct, random_seed = stage_seeds(seed_of(i, s))[["size"]], strict = FALSE, verbose = FALSE)
    cd <- dg$cell_df; it <- z$cells$interior; nuc <- dg$cell_nuclei
  }
  nc <- cd$nucleus_area / cd$territory_area
  out <- data.frame(archetype = A, roi_id = roi, seed = s, arm = j$arm, density = dens, n = nrow(cd), median_area_nuc = median(cd$nucleus_area[it]),
                    median_nc = median(nc[it]), frac_nc095 = mean(nc[it] >= 0.95), median_area_terr = median(cd$territory_area[it]),
                    cv_terr = sd(cd$territory_area[it]) / mean(cd$territory_area[it]))
  if (paste(A, roi) %in% WIN_ROI) { g <- geos::as_geos_geometry(nuc)
    write_parquet(data.frame(idx = seq_along(g), wkb = I(geos::geos_write_wkb(g))), file.path(OUT, "nuc", j$arm, sprintf("%s_%s_s%02d.parquet", A, roi, s))) }
  saveRDS(out, f)
}, "posthoc")
tab <- do.call(rbind, lapply(list.files(file.path(ph_dir, "units"), recursive = TRUE, full.names = TRUE), readRDS))
write.csv(tab, file.path(RES, "S1.5_posthoc_units.csv"), row.names = FALSE); cat("unita' post hoc:", nrow(tab), "\n")
