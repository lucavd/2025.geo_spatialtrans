# tools/riesame_B_sim.R - riesame, Parte B: braccio B (banale), definito in docs/riesame/REGOLA_B.md; scelte di attuazione in
# docs/riesame/DECISIONI_B.md (righe B-1...B-12). Usa le funzioni del pacchetto (R/04b*, R/05b*) senza modificarle.
#   Strato cellulare: centri da un processo di Poisson omogeneo con la densita' del preset (.TP_TABLE$density, nuclei/mm^2) in ogni
#   regione dell'area di S (regioni dei livelli U6 di S1.6, le stesse 5 geometrie per tessuto); territori di Voronoi
#   (tessellate_voronoi(), corner_smoothing = 0 come S); nucleo = poligono regolare di 64 lati con area uguale all'area nucleare
#   media del preset (media dei nuclei dei 5 livelli U6 del tessuto), centrato nel centro e ritagliato nel territorio (senza
#   riportarlo all'area); tipi estratti iid dalla composizione del preset (tissue_preset_spec()$composition).
#   Step 2: build_hd_grid(); simulate_stripes() con capture_params() del tessuto; simulate_cell_profiles() con un profilo per tipo
#   (use_subtypes = FALSE, media dei sottotipi), conteggi multinomiali (noise = "M"), profondita' del preset visium_hd_depth_factor()
#   per la sovra-generazione M delle strisce; place_and_move_transcripts() con r = 1 (nascita uniforme nel territorio), campo spento,
#   spostamento spento (w = 0, w2 = 0, sigma_in = sigma_e = 1e-9 um: la funzione chiede sigma > 0); capture_and_bin() con il fondo b
#   del tessuto. Uscita con i campi letti da tools/F3.1_export_sim.R.
# Semi: tabella Mersenne-Twister estratta dal seme base sulle 10 aree ordinate (A2_r1...A3_r5) x 2 semi (k = 0, 1); per ogni
#   esecuzione 6 semi dei passi (poisson, types, stripes, profiles, place, capture) con sample.int(1e8, 6) dal seme dell'esecuzione.
# Uso (dalla radice del repository): Rscript --vanilla tools/riesame_B_sim.R <out_dir> <seed_base> <tags,> <workers> <k list>
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos); library(parallel) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
a <- commandArgs(TRUE); out_dir <- a[1]; seed_base <- as.integer(a[2]); tags <- strsplit(a[3], ",")[[1]]; W <- as.integer(a[4])
ks <- as.integer(strsplit(a[5], ",")[[1]])
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
lay_dir <- "/mnt/micron/geo_spatialtrans/S1.6/layers"; NV <- 64L
ALL <- c(paste0("A2_r", 1:5), paste0("A3_r", 1:5))
set.seed(seed_base, kind = "Mersenne-Twister")
seed_tab <- matrix(sample.int(.Machine$integer.max, length(ALL) * 2L), ncol = 2L, dimnames = list(ALL, NULL))
write.csv(data.frame(tag = rep(ALL, 2), k = rep(0:1, each = length(ALL)), seed = c(seed_tab)), file.path(out_dir, "seed_table.csv"), row.names = FALSE)
preset_of <- function(A) .TP_TABLE$preset[.TP_TABLE$archetype == A]
# area nucleare media del preset: media di cell_df$nucleus_area sui 5 livelli U6 del tessuto (B-4)
nuc_mean <- sapply(c("A2", "A3"), function(A) mean(unlist(lapply(1:5, function(k)
  readRDS(file.path(lay_dir, sprintf("%s_r%d_U6.rds", A, k)))$cell_df$nucleus_area))))
poisson_in <- function(poly, lam) {   # processo di Poisson omogeneo di intensita' lam (per um^2) nel poligono (sfc di lunghezza 1)
  Ar <- as.numeric(sf::st_area(poly)); n <- stats::rpois(1L, lam * Ar); bb <- sf::st_bbox(poly)
  Ab <- (bb[["xmax"]] - bb[["xmin"]]) * (bb[["ymax"]] - bb[["ymin"]]); x <- numeric(0); y <- numeric(0)
  while (length(x) < n) {
    m <- as.integer(ceiling((n - length(x)) * 1.2 * Ab / Ar)) + 10L
    px <- stats::runif(m, bb[["xmin"]], bb[["xmax"]]); py <- stats::runif(m, bb[["ymin"]], bb[["ymax"]])
    pts <- sf::st_as_sfc(sprintf("POINT (%.17g %.17g)", px, py)); ins <- sf::st_intersects(poly, pts)[[1]]
    ok <- seq_len(m) %in% ins; x <- c(x, px[ok]); y <- c(y, py[ok])
  }
  list(x = x[seq_len(n)], y = y[seq_len(n)])
}
b_layer <- function(tag, sd) {
  A <- substr(tag, 1, 2); pre <- preset_of(A); L <- readRDS(file.path(lay_dir, paste0(tag, "_U6.rds")))
  dens <- .TP_TABLE$density[.TP_TABLE$archetype == A]; comp <- tissue_preset_spec(pre)$composition
  set.seed(sd[["poisson"]], kind = "Mersenne-Twister")
  P <- lapply(seq_len(nrow(L$region_df)), function(r) poisson_in(L$region_polygons[r], dens / 1e6))
  rid <- rep(L$region_df$region_id, vapply(P, function(p) length(p$x), 1L))
  cen <- data.frame(cell_id = seq_along(rid), region_id = rid, x = unlist(lapply(P, `[[`, "x")), y = unlist(lapply(P, `[[`, "y")))
  set.seed(sd[["types"]], kind = "Mersenne-Twister")
  ctype <- sample(comp$cell_type, nrow(cen), replace = TRUE, prob = comp$fraction)
  tv <- tessellate_voronoi(cen, list(region_df = L$region_df, region_polygons = L$region_polygons), corner_smoothing = 0, verbose = FALSE)
  th <- 2 * pi * (0:(NV - 1L)) / NV; R <- sqrt(nuc_mean[[A]] / pi) * .dcg_poly_factor(NV)   # area del 64-gono = area media
  wkt <- vapply(seq_len(nrow(cen)), function(i) { xx <- cen$x[i] + R * cos(th); yy <- cen$y[i] + R * sin(th)
    paste0("POLYGON ((", paste(sprintf("%.17g %.17g", c(xx, xx[1]), c(yy, yy[1])), collapse = ", "), "))") }, "")
  gn <- geos::geos_intersection(geos::geos_read_wkt(wkt), geos::as_geos_geometry(sf::st_set_crs(tv$cell_territories, NA)))
  nuc <- sf::st_as_sfc(geos::geos_write_wkb(gn)); gt <- as.character(sf::st_geometry_type(nuc))
  if (any(!gt %in% c("POLYGON", "MULTIPOLYGON"))) nuc <- sf::st_collection_extract(nuc, "POLYGON")
  if (any(as.character(sf::st_geometry_type(nuc)) == "MULTIPOLYGON")) nuc <- sf::st_cast(nuc, "MULTIPOLYGON")
  a_n <- as.numeric(sf::st_area(nuc)); a_t <- as.numeric(sf::st_area(tv$cell_territories))
  cd <- data.frame(cell_id = cen$cell_id, region_id = cen$region_id, cell_type = ctype, x = cen$x, y = cen$y, territory_area = a_t,
                   nucleus_cx = cen$x, nucleus_cy = cen$y, nucleus_area = a_n, nucleus_area_target = nuc_mean[[A]],
                   nucleus_clipped = a_n < nuc_mean[[A]] - 1e-6, stringsAsFactors = FALSE)
  list(cell_df = cd, region_df = L$region_df, region_polygons = L$region_polygons, cell_territories = tv$cell_territories, cell_nuclei = nuc,
       territory_df = tv$territory_df, region_check = tv$region_check,
       metadata = list(pixel_size_um = L$metadata$pixel_size_um, tissue_preset_used = pre, random_seed = sd[["poisson"]], n_cells = nrow(cd),
                       n_regions = nrow(L$region_df), package_version = .geo_version(), arm = "riesame B",
                       B = list(density_mm2 = dens, nucleus_area_mean = nuc_mean[[A]], n_vertices = NV, composition = comp,
                                tessellation_info = tv$info, n_nuclei_multipart = sum(sf::st_geometry_type(nuc) == "MULTIPOLYGON" & lengths(nuc) > 1))))
}
run1 <- function(tag, k) {
  s <- seed_tab[tag, k + 1L]; fo <- file.path(out_dir, sprintf("%s_k%d.rds", tag, k)); if (file.exists(fo)) return(paste("skip", tag, k))
  t0 <- proc.time()[["elapsed"]]; A <- substr(tag, 1, 2)
  set.seed(s, kind = "Mersenne-Twister"); sd <- setNames(sample.int(1e8, 6L), c("poisson", "types", "stripes", "profiles", "place", "capture"))
  BL <- b_layer(tag, sd); t1 <- proc.time()[["elapsed"]]
  hd <- build_hd_grid(BL)
  cp <- capture_params(A); p <- transcript_params(A)
  p$r <- 1; p$w <- 0; p$w2 <- 0; p$ell2 <- cp$ell2; p$sigma_in <- 1e-9; p$sigma_e <- 1e-9
  St <- simulate_stripes(hd, cp, seed = sd[["stripes"]]); M <- St$M
  profiles <- load_cell_profiles(A); df0 <- visium_hd_depth_factor(A)
  prof <- simulate_cell_profiles(BL$cell_df[, c("cell_id", "cell_type")], profiles, seed = sd[["profiles"]], depth_factor = df0 * M,
                                 noise = "M", use_subtypes = FALSE)
  cs <- visium_hd_chunk_size(Matrix::colSums(prof$counts))
  tx <- place_and_move_transcripts(BL, hd, prof, p, seed = sd[["place"]], field = FALSE, frame_open = TRUE, chunk_size = cs)
  cap <- capture_and_bin(tx, hd, St, list(b = cp$b), seed = sd[["capture"]], overgeneration = M, background = TRUE)
  # controllo: senza spostamento ogni trascritto catturato cade in un quadratino che tocca il territorio della sua cellula
  og <- cap$truth$origin; ov <- hd$truth$overlap
  own <- paste(og$bin_id, og$cell_id) %in% paste(ov$bin_id[ov$area_territory_um2 > 0], ov$cell_id[ov$area_territory_um2 > 0])
  chk <- list(n_cells = nrow(BL$cell_df), density_realised_mm2 = nrow(BL$cell_df) / sum(as.numeric(sf::st_area(BL$region_polygons))) * 1e6,
              type_share = as.list(table(BL$cell_df$cell_type) / nrow(BL$cell_df)), nucleus_area_mean_realised = mean(BL$cell_df$nucleus_area),
              share_clipped = mean(BL$cell_df$nucleus_clipped), share_captured_in_own_territory_bin = sum(og$n[own]) / sum(og$n),
              n_lost = sum(cap$truth$cells$n_lost), M = M, depth_factor = df0, chunk_size = cs)
  out <- list(counts_2um = cap$counts_2um, bins = cap$bins,
              truth = list(overlap = ov, origin = og, background = cap$truth$background, cell_counts = prof$counts, cells = prof$cells),
              cell_layer = BL[c("cell_df", "region_df", "territory_df", "region_check", "metadata")],
              f31 = list(tag = tag, k = k, seed = s, seed_base = seed_base, step_seeds = sd, archetype = A, arm = "B",
                         elapsed = proc.time()[["elapsed"]] - t0, elapsed_layer = t1 - t0, check = chk,
                         params = list(transcripts = p, capture = cp, field = FALSE, frame_open = TRUE, noise = "M", use_subtypes = FALSE)))
  saveRDS(out, paste0(fo, ".tmp")); file.rename(paste0(fo, ".tmp"), fo)
  jsonlite::write_json(c(list(tag = tag, k = k, seed = s, elapsed = out$f31$elapsed), chk), sub("\\.rds$", "_check.json", fo), auto_unbox = TRUE, digits = NA)
  sprintf("%s k%d seed %d: %.0f s, %d cells, own-territory share %.6f", tag, k, s, out$f31$elapsed, chk$n_cells, chk$share_captured_in_own_territory_bin)
}
Sys.setenv(OMP_NUM_THREADS = "1"); if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1L)
jobs <- expand.grid(tag = tags, k = ks, stringsAsFactors = FALSE)
res <- mclapply(seq_len(nrow(jobs)), function(i) tryCatch(run1(jobs$tag[i], jobs$k[i]), error = function(e) paste("ERROR", jobs$tag[i], jobs$k[i], conditionMessage(e))),
                mc.cores = W, mc.preschedule = FALSE)
cat(unlist(res), sep = "\n")
