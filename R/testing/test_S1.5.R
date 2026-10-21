# R/testing/test_S1.5.R - S1.5: test di orchestratore, preset, controlli, scala, figura, casi limite (cancello G-1; mutanti G-2).
# Pre-registrazione: docs/steps/S1.5.md (commit 6042a89, 710f5df). Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S15_MUTANT=none|M1..M9] Rscript --vanilla R/testing/test_S1.5.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
source("tools/S1.5_faults.R")
MUT <- Sys.getenv("S15_MUTANT", "none")
VIS <- "results/visualizations/S1.5"; dir.create(VIS, recursive = TRUE, showWarnings = FALSE)   # BL-013: fuori dai file tracciati
res <- list()
ck <- function(name, ok, val = "") { ok <- isTRUE(ok); res[[length(res) + 1]] <<- ok
  cat(sprintf("%s %s %s\n", if (ok) "PASS" else "FAIL", name, paste(val, collapse = " "))) }
safe <- function(expr) tryCatch(expr, error = function(e) structure(list(err = conditionMessage(e)), class = "t_err"))
isok <- function(o) !inherits(o, "t_err"); msg <- function(o) if (isok(o)) "" else o$err
st_of <- function(v, id) v$table$status[v$table$check == id]
wkb_md5 <- function(L) { f <- tempfile(); saveRDS(list(L$cell_df, sf::st_as_binary(L$cell_territories), sf::st_as_binary(L$cell_nuclei)), f, compress = FALSE)
  unname(tools::md5sum(f)) }
toy <- function(nx = 200, ny = 150) { cl <- expand.grid(x = 1:nx, y = 1:ny); cl$intensity_cluster <- ifelse(cl$x < nx * 0.45, 1L, 2L); cl }
PC <- function(...) place_cells(..., verbose = FALSE)

# ---- O: orchestratore ----
cl0 <- toy()
L <- safe(PC(cl0, "epithelial", pixel_size_um = 1, random_seed = 11))
ck("O1 place_cells su tessuto di prova: nessun errore, controlli senza ERROR", isok(L) && L$metadata$validation$ok, msg(L))
if (isok(L)) {
  need <- c("cell_id", "region_id", "cluster_id", "cell_type", "x", "y", "territory_area", "eq_radius", "nucleus_shape", "nucleus_radius",
            "nucleus_major_r", "nucleus_minor_r", "nucleus_orient", "nucleus_cx", "nucleus_cy", "nucleus_area", "cytoplasm_area",
            "is_minority_type", "nucleus_pale")
  ck("O2 cell_df con le colonne del design §4.2 + is_minority_type + nucleus_pale", all(need %in% names(L$cell_df)), setdiff(need, names(L$cell_df)))
  ck("O3 cell_layer con gli elementi del design §4.1", all(c("cell_df", "region_df", "region_polygons", "cell_territories", "cell_nuclei", "metadata") %in% names(L)))
  mk <- c("pixel_size_um", "pixel_size_source", "n_cells", "n_regions", "tissue_preset_used", "corner_smoothing", "cell_types_resolved",
          "region_composition_resolved", "random_seed", "package_version", "timestamp", "validation")
  ck("O4 metadata con le chiavi del design §4.1", all(mk %in% names(L$metadata)), setdiff(mk, names(L$metadata)))
  n <- nrow(L$cell_df)
  ck("O5 cell_id 1..N nello stesso ordine in cell_df, territori e nuclei", identical(L$cell_df$cell_id, seq_len(n)) && length(L$cell_territories) == n &&
     length(L$cell_nuclei) == n && all(abs(as.numeric(st_area(L$cell_territories)) - L$cell_df$territory_area) <= 1e-9 * L$cell_df$territory_area))
  ck("O6 region_df con dominant_cell_type, achieved_density, n_cells", all(c("dominant_cell_type", "achieved_density", "n_cells", "target_density_weighted") %in% names(L$region_df)))
  ss <- stage_seeds(11)
  ck("O7 semi degli stadi distinti (semina, taglia, pallidi)", length(unique(ss)) == 3, paste(ss, collapse = " "))
  ck("O8 is_minority_type tutto FALSE con un tipo per tessuto", !any(L$cell_df$is_minority_type))
  pf <- L$cell_types_resolved$pale_fraction; fr <- mean(L$cell_df$nucleus_pale)
  ck("O9 nuclei pallidi: frazione osservata entro 4 errori binomiali dalla frazione del preset", abs(fr - pf) <= 4 * sqrt(pf * (1 - pf) / n), sprintf("%.3f vs %.3f", fr, pf))
  # equivalenza con la catena manuale
  sp <- tissue_preset_spec("epithelial"); s <- stage_seeds(11)
  cm <- merge_small_components(cl0, 1, 100, "intensity_cluster", 1L)
  rg <- extract_regions(cm, 1, 100, cluster_col = "intensity_cluster", simplify_tol_um = 0.5, stride = 1L, verbose = FALSE)
  comp <- data.frame(cluster_id = rg$region_df$cluster_id, cell_type = sp$cell_types$cell_type, fraction = 1); comp <- unique(comp)
  sc <- seed_centroids(rg, sp$cell_types, comp, random_seed = s[["seed"]], verbose = FALSE)
  tv <- tessellate_voronoi(sc$centroids, rg, verbose = FALSE)
  dg <- derive_cell_geometry(sc$centroids, tv, sp$cell_types, random_seed = s[["size"]], verbose = FALSE)
  same <- identical(dg$cell_df, L$cell_df[, names(dg$cell_df)]) && identical(st_as_binary(dg$cell_nuclei), st_as_binary(L$cell_nuclei)) &&
    identical(st_as_binary(tv$cell_territories), st_as_binary(L$cell_territories))
  ck("O10 place_cells identico alla catena manuale con i semi documentati", same)
  L2 <- PC(cl0, "epithelial", pixel_size_um = 1, random_seed = 11)
  ck("O11 stesso seme -> uscita identica (md5 di cell_df e WKB)", wkb_md5(L) == wkb_md5(L2))
  L3 <- PC(cl0, "epithelial", pixel_size_um = 1, random_seed = 12)
  ck("O12 seme diverso -> n_cells entro ±5 % (design §6.2)", abs(nrow(L3$cell_df) / n - 1) <= 0.05, sprintf("%d vs %d", nrow(L3$cell_df), n))
  ck("O13 seme diverso -> posizioni diverse", !identical(L3$cell_df$x[1:10], L$cell_df$x[1:10]))
}

# ---- P: preset contro BIO_REFERENCES (letti dal file) ----
bio <- readLines("docs/BIO_REFERENCES.md"); bio <- bio[grepl("^\\| B-[0-9]{3} \\|", bio)]
bval <- function(id) { r <- bio[startsWith(bio, paste0("| ", id, " |"))]; trimws(strsplit(r, "\\|")[[1]][5]) }
bnum <- function(id) as.numeric(gsub("[^0-9.]", "", bval(id)))
binterval <- function(id) { r <- bio[startsWith(bio, paste0("| ", id, " |"))]; f <- trimws(strsplit(r, "\\|")[[1]]); x <- regmatches(f[6], regexpr("\\[[0-9., ]+\\]", f[6]))
  as.numeric(strsplit(gsub("[][ ]", "", x), ",")[[1]]) }
# tetto phi: B-121...B-126 lo riportano nelle note ("q99 N/C x", con tetto 1: DEC-005), B-143 e B-154 nel valore
bphi <- function(id) { v <- suppressWarnings(as.numeric(bval(id))); if (is.finite(v)) return(v)
  r <- bio[startsWith(bio, paste0("| ", id, " |"))]; m <- regmatches(r, regexpr("q99 N/C [0-9.]+", r)); min(1, as.numeric(sub("q99 N/C ", "", m))) }
cal <- read.csv("results/S1.4/S1.4_calibration.csv"); r3b <- read.csv("results/R3b/R3b_refs.csv")
dm_ref <- c(cal$density_mm2[cal$archetype == "A1" & cal$method == "cellpose_rgb"], cal$density_mm2[cal$archetype == "A2" & cal$method == "cellpose_rgb"],
            cal$density_mm2[cal$archetype == "A3" & cal$method == "cellpose_rgb"], r3b$consensus[r3b$archetype == "A4" & r3b$metric == "density_mm2"],
            cal$density_mm2[cal$archetype == "A5" & cal$method == "cellpose_rgb"], r3b$consensus[r3b$archetype == "A6" & r3b$metric == "density_mm2"])
bad <- character(0)
for (k in seq_along(tissue_presets())) {
  p <- tissue_preset_spec(tissue_presets()[k]); s <- p$sources; v <- setNames(s$value, s$parameter); src <- setNames(s$source, s$parameter)
  chk <- function(par, ref) if (!isTRUE(abs(v[[par]] - ref) <= 1e-9 * max(1, abs(ref)))) bad <<- c(bad, sprintf("%s:%s %g vs %g", p$preset, par, v[[par]], ref))
  chk("density", bnum(src[["density"]])); chk("density_evident", binterval(src[["density"]])[1])
  for (par in c("min_dist_um", "rho", "beta", "sigma")) chk(par, bnum(src[[par]]))
  chk("phi", bphi(src[["phi"]]))
  iv <- if (grepl("^\\[", bval(src[["c8_lo"]]))) as.numeric(strsplit(gsub("[][ ]", "", bval(src[["c8_lo"]])), ",")[[1]]) else NA
  chk("c8_lo", iv[1]); chk("c8_hi", iv[2]); chk("density_model", round(dm_ref[k], 1))
  ct <- p$cell_types
  rho_exp <- v[["rho"]] * (v[["density"]] / v[["density_model"]])^((1 - v[["beta"]]) / 2)
  if (abs(ct$nucleus_to_eq_ratio - rho_exp) > 1e-12) bad <- c(bad, paste(p$preset, "rho'"))
  if (abs(ct$pale_fraction - (1 - v[["density_evident"]] / v[["density"]])) > 1e-12) bad <- c(bad, paste(p$preset, "pale"))
  if (abs(tissue_preset_spec(tissue_presets()[k], "none")$cell_types$nucleus_to_eq_ratio - v[["rho"]]) > 1e-12) bad <- c(bad, paste(p$preset, "rho none"))
}
ck("P1 ogni valore dei 6 preset uguale alla riga di BIO_REFERENCES o al file di calibrazione che cita; rho' e frazione pallida dalle formule", !length(bad), bad)
ck("P2 preset sconosciuto -> errore", !isok(safe(tissue_preset_spec("liver"))))

# ---- V: controlli su un caso buono e su guasti costruiti ----
if (isok(L)) {
  V0 <- validate_cell_layer(L, c8_interval = c(2.81, 7.71))
  ck("V0 caso buono: tutti PASS", all(V0$table$status == "PASS"), paste(V0$table$check[V0$table$status != "PASS"]))
  npr <- nrow(.v5_pairs(L$cell_territories))
  ck("V0b C2 valuta le coppie nei due ordini (BL-079)", grepl(sprintf("intersezioni valutate %d", 2 * npr), V0$table$note[V0$table$check == "C2"]), V0$table$note[V0$table$check == "C2"])
  for (ty in FAULT_TYPES) for (d in c(0.01, 1)) {
    fz <- safe(make_fault(L, ty, d))
    if (!isok(fz)) { ck(sprintf("V-%s-%g guasto costruito", ty, d), FALSE, msg(fz)); next }
    vv <- safe(validate_cell_layer(fz$cl, c8_interval = c(2.81, 7.71)))
    ck(sprintf("V-%s-%g rilevato da %s (ERROR)", ty, d, fz$expect), isok(vv) && identical(st_of(vv, fz$expect), "ERROR"), if (isok(vv)) paste(vv$table$check, vv$table$status, collapse = " ") else msg(vv))
  }
  Lg <- L; Lg$region_check$gap_area[1] <- NaN
  ck("V-nan-gap NaN nella somma di C2 -> ERROR (BL-095)", identical(st_of(validate_cell_layer(Lg), "C2"), "ERROR"))
  L0 <- L; L0$region_df <- rbind(L0$region_df, transform(L0$region_df[1, ], region_id = 99L)); L0$region_polygons <- c(L0$region_polygons, L0$region_polygons[1] + c(1000, 0)); L0$region_df$n_cells[nrow(L0$region_df)] <- 0L
  v0 <- validate_cell_layer(L0, min_region_area_um2 = 0)
  ck("V-empty regione senza cellule: C6 WARN, C2 non la conta", identical(st_of(v0, "C6"), "WARN") && identical(st_of(v0, "C2"), "PASS"))
  Lo <- L; Lo$cell_nuclei <- Lo$cell_nuclei[c(2:length(Lo$cell_nuclei), 1)]
  ck("V-order nuclei fuori ordine -> C3 ERROR", identical(st_of(validate_cell_layer(Lo), "C3"), "ERROR"))
  Lc <- L; Lc$cell_df$cell_id <- rev(Lc$cell_df$cell_id)
  ck("V-cellid cell_id non 1..N -> C9 ERROR", identical(st_of(validate_cell_layer(Lc), "C9"), "ERROR"))
  L1 <- L; L1$region_df$achieved_density[1] <- 2 * L1$region_df$target_density_weighted[1]
  ck("V-C1 densita' doppia -> WARN", identical(st_of(validate_cell_layer(L1), "C1"), "WARN"))
  L4 <- L; L4$region_composition <- rbind(L4$region_composition, data.frame(cluster_id = L4$region_df$cluster_id[1], cell_type = "altro", fraction = 1))
  ck("V-C4 frazione di tipo mancata -> WARN", identical(st_of(validate_cell_layer(L4), "C4"), "WARN"))
  L6 <- L; L6$region_df$n_cells[1] <- 0L
  ck("V-C6 regione senza cellule -> WARN", identical(st_of(validate_cell_layer(L6, min_region_area_um2 = 0), "C6"), "WARN"))
  L7 <- L; L7$cell_df$eq_radius <- 3 * L7$cell_df$eq_radius
  ck("V-C7 eq_r triplo -> WARN", identical(st_of(validate_cell_layer(L7), "C7"), "WARN"))
}
# overlay_prec
a <- st_sfc(st_polygon(list(rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0))))) + c(5000, 5000)
b <- st_sfc(st_polygon(list(rbind(c(5, 5), c(20, 5), c(20, 20), c(5, 20), c(5, 5))))) + c(5000, 5000)
ok_ov <- abs(overlay_prec("intersection", a, b) - 25) < 1e-9 && abs(overlay_prec("difference", a, b) - 75) < 1e-9 &&
  abs(overlay_prec("union", a, b) - 300) < 1e-9 && abs(as.numeric(st_area(overlay_prec("intersection", a, b, area = FALSE))) - 25) < 1e-9 &&
  all(abs(st_bbox(overlay_prec("intersection", a, b, area = FALSE))[c("xmin", "ymin")] - 5005) < 1e-9)
ck("V-ov overlay_prec: aree esatte, geometria via WKB in coordinate originali (BL-116)", ok_ov)

# ---- S: scala ----
sj <- list.files("/mnt/micron/geo_spatialtrans/data_real", pattern = "^scalefactors_json.json$", recursive = TRUE, full.names = TRUE)
sj <- sj[grepl("square_002um", sj)]
okm <- length(sj) == 5 && all(vapply(sj, function(f) { e <- estimate_pixel_size(NULL, NULL, method = "metadata", scalefactors_json = f)
  abs(e$pixel_size_um - jsonlite::fromJSON(f)$microns_per_pixel) <= 1e-9 && e$pixel_size_um > 0.2 && e$pixel_size_um < 0.6 }, TRUE))
ck("S1 metodo dai metadati sui 5 dataset reali = microns_per_pixel (entro 1e-9), controllato sulla griglia dei bin", okm, length(sj))
# immagine sintetica con nuclei a scala nota: dischi scuri di raggio 3 µm a 4 000 / mm², 0.5 µm/px
disk_img <- function(s_um, n_side_um = 400, dens = 4000, r_um = 3, seed = 5) {
  set.seed(seed); npx <- round(n_side_um / s_um); nn <- round(dens * (n_side_um / 1000)^2); P <- matrix(0, 0, 2); tries <- 0
  while (nrow(P) < nn && tries < 1e6) { tries <- tries + 1; q <- runif(2, 0, n_side_um); if (!nrow(P) || min(colSums((t(P) - q)^2)) > (2.6 * r_um)^2) P <- rbind(P, q) }
  g <- expand.grid(x = 1:npx, y = 1:npx); cx <- (g$x - 0.5) * s_um; cy <- (g$y - 0.5) * s_um; v <- rep(0.82, nrow(g))
  for (k in seq_len(nrow(P))) { w <- (cx - P[k, 1])^2 + (cy - P[k, 2])^2 <= r_um^2; v[w] <- 0.35 }
  g$value <- v; g$intensity_cluster <- 1L; g }
ctS <- data.frame(cell_type = "t", density = 4000)
e1 <- safe(estimate_pixel_size(disk_img(0.5), ctS))
ck("S2 dischi a 0.5 µm/px: scala stimata entro ±20 %, IC presente", isok(e1) && !e1$refused && abs(e1$pixel_size_um / 0.5 - 1) <= 0.2 && all(is.finite(e1$ci)),
   if (isok(e1)) sprintf("%.4g (N %d, sigma* %.2f)", e1$pixel_size_um, e1$n_features, e1$sigma_px) else msg(e1))
smooth <- expand.grid(x = 1:300, y = 1:300); smooth$value <- 0.4 + 0.3 * (smooth$x / 300); smooth$intensity_cluster <- 1L
e2 <- safe(estimate_pixel_size(smooth, ctS))
ck("S3 tessuto liscio senza nuclei: rifiuto, nessun numero", isok(e2) && isTRUE(e2$refused) && is.na(e2$pixel_size_um), if (isok(e2)) e2$reason else msg(e2))
e3 <- safe(estimate_pixel_size(disk_img(10, n_side_um = 2000, seed = 6), ctS))
ck("S4 dischi a 10 µm/px (nuclei sotto il pixel): rifiuto", isok(e3) && isTRUE(e3$refused), if (isok(e3)) e3$reason else msg(e3))
de <- disk_img(0.5); Le <- safe(withCallingHandlers(PC(de, cell_types = data.frame(cell_type = "t", density = 4000, nucleus_to_eq_ratio = 0.4), pixel_size_um = "estimate", random_seed = 3),
  warning = function(w) invokeRestart("muffleWarning")))
ck("S5 place_cells(pixel_size_um = \"estimate\"): fonte \"estimated\" in metadata, scala entro ±20 %", isok(Le) && Le$metadata$pixel_size_source == "estimated" &&
   abs(Le$metadata$pixel_size_um / 0.5 - 1) <= 0.2, if (isok(Le)) Le$metadata$pixel_size_um else msg(Le))
w5 <- tryCatch({ PC(de, cell_types = data.frame(cell_type = "t", density = 4000, nucleus_to_eq_ratio = 0.4), pixel_size_um = "estimate", random_seed = 3); "" },
  warning = function(w) conditionMessage(w), error = function(e) "")
ck("S6 stima della scala mai silenziosa: avviso", grepl("scala stimata", w5), w5)
Ls <- safe(PC(smooth, cell_types = ctS, pixel_size_um = "estimate"))
ck("S7 stima rifiutata -> place_cells si ferma con errore", !isok(Ls) && grepl("rifiutata", msg(Ls)), msg(Ls))

# ---- F: figura ----
if (isok(L)) {
  okp <- all(vapply(c("voronoi", "types", "density", "regions"), function(m) { p <- plot_cell_layer(L, m)
    f <- file.path(VIS, paste0("plot_", m, ".png")); ggplot2::ggsave(f, p, width = 5, height = 4, dpi = 80, bg = "white"); inherits(p, "ggplot") && file.exists(f) }, TRUE))
  ck("F1 plot_cell_layer: 4 modi, PNG fuori dai file tracciati (BL-013)", okp)
}

# ---- E: casi limite del design §6.2 ----
ce <- toy(); ce$intensity_cluster[ce$x %in% 30:34 & ce$y %in% 30:34] <- 2L            # 25 µm² di cluster 2 dentro il cluster 1
Ex <- PC(ce, "epithelial", pixel_size_um = 1, small_components = "exclude", random_seed = 2)
ck("E1 regione sotto la soglia di area: esclusa (exclude)", nrow(Ex$region_df) == 2 && nrow(Ex$metadata$excluded_df) == 1)
Em <- PC(ce, "epithelial", pixel_size_um = 1, random_seed = 2)
ck("E1b BL-043: con merge la componente piccola passa al cluster adiacente", nrow(Em$region_df) == 2 && nrow(Em$metadata$excluded_df) == 0 && Em$metadata$merge_info$n_px_moved == 25)
ci <- toy(); ci <- ci[!(ci$x %in% 60:70 & ci$y %in% 60:70) | (ci$x %in% 64:66 & ci$y %in% 64:66), ]; ci$intensity_cluster[ci$x %in% 64:66 & ci$y %in% 64:66] <- 2L
Ei <- PC(ci, "epithelial", pixel_size_um = 1, random_seed = 2)
ck("E1c BL-043: isola che tocca solo il fondo resta esclusa", nrow(Ei$metadata$excluded_df) == 1)
ctu <- data.frame(cell_type = c("a", "b"), density = c(5000, 3000), nucleus_to_eq_ratio = c(0.4, 0.5))
w2 <- tryCatch({ PC(toy(), cell_types = ctu, region_composition = data.frame(cluster_id = c(1, 1, 2), cell_type = c("a", "b", "a"), fraction = c(0.5, 0.3, 1)), pixel_size_um = 1); "" },
  warning = function(w) conditionMessage(w))
ck("E2 frazioni che non sommano a 1: normalizzazione con avviso", nchar(w2) > 0, w2)
e3 <- safe(PC(toy(), cell_types = ctu, region_composition = data.frame(cluster_id = 1:2, cell_type = c("a", "zz"), fraction = 1), pixel_size_um = 1))
ck("E3 tipo nella composizione ma non nel catalogo: errore puntuale", !isok(e3) && grepl("zz", msg(e3)), msg(e3))
big <- expand.grid(x = 1:4000, y = 1:4000); big$intensity_cluster <- 1L
gc(reset = TRUE); tb <- system.time(Lb <- safe(PC(big, "columnar_muscle", pixel_size_um = 1, random_seed = 4)))[["elapsed"]]; mb <- sum(gc()[, 6])
ck("E4 tessuto enorme 4000 x 4000 px (cuore): nessun errore, < 5 min, memoria R < 8 GB (il tempo nel linfonodo e' la previsione P-8)",
   isok(Lb) && tb < 300 && mb < 8000, sprintf("%.1f s, %.0f MB, %s celle", tb, mb, if (isok(Lb)) nrow(Lb$cell_df) else msg(Lb)))
rm(big, Lb); invisible(gc())
c1 <- expand.grid(x = 1:120, y = 1:120); c1$intensity_cluster <- 3L
E5 <- PC(c1, "lymphoid", pixel_size_um = 1, random_seed = 2)
ck("E5 tessuto con un solo cluster: una regione", nrow(E5$region_df) == 1 && E5$metadata$validation$ok)
E6 <- safe(PC(toy(), "lymphoid", pixel_size_um = 1, corner_smoothing = 1, random_seed = 2))
ck("E6 corner_smoothing = 1 con territori piccoli: nessun poligono degenere, controlli senza ERROR", isok(E6) && all(E6$cell_df$territory_area > 0) && E6$metadata$validation$ok &&
   sum(E6$region_check$gap_area) > 0, msg(E6))
E7 <- PC(toy(), cell_types = ctu, region_composition = data.frame(cluster_id = c(1, 1, 2, 2), cell_type = c("a", "b", "a", "b"), fraction = c(1, 0, 0.7, 0.3)), pixel_size_um = 1, random_seed = 2)
ck("E7 fraction = 0: quel tipo non e' piazzato in quel cluster", !any(E7$cell_df$cell_type == "b" & E7$cell_df$cluster_id == 1) && any(E7$cell_df$cell_type == "b"))
ck("E7b is_minority_type vero solo per il tipo minoritario", all(E7$cell_df$is_minority_type == (E7$cell_df$cell_type == "b")))
e8 <- safe(place_cells(toy(), "epithelial", verbose = FALSE))
ck("E8 pixel_size_um omesso: errore puntuale (C-A)", !isok(e8) && grepl("obbligatorio", msg(e8)), msg(e8))
ck("E9 \"estimate\" su immagine a scala nota: vedi S2, S5, S6 (controprova CP-1 per le immagini reali)", isok(e1) && isok(Le))
E10 <- PC(toy(), "epithelial", pixel_size_um = 10, random_seed = 2)
ck("E10 scala sbagliata di 10 volte: C8 scatta (WARN con scala dell'utente)", identical(st_of(E10$metadata$validation, "C8"), "WARN"))
ck("E10b con scala stimata C8 e' ERROR", identical(st_of(validate_cell_layer(E10, c8_interval = c(2.81, 7.71), pixel_size_source = "estimated"), "C8"), "ERROR"))
# aggiunti dopo il primo giro (nota D-3 dello step): cosa C8 vede davvero
E10c <- PC(toy(2000, 1500), cell_types = data.frame(cell_type = "u", density = 12.419, nucleus_to_eq_ratio = 0.4), pixel_size_um = 1, random_seed = 2)
ck("E10c C8 scatta con la densita' del catalogo in unita' sbagliate (cellule/mm² scambiate per cellule/µm² x 1e-3)", identical(st_of(E10c$metadata$validation, "C8"), "WARN"),
   E10c$metadata$validation$table$value[E10c$metadata$validation$table$check == "C8"])
E10d <- PC(toy(), "epithelial", pixel_size_um = 1, random_seed = 2)
ck("E10d spiegazione: con la densita' in unita' assolute l'eq_r mediano non dipende dalla scala (px 1 vs 10 entro il 5 %)",
   abs(median(E10$cell_df$eq_radius) / median(E10d$cell_df$eq_radius) - 1) <= 0.05, sprintf("%.3f vs %.3f um", median(E10$cell_df$eq_radius), median(E10d$cell_df$eq_radius)))
E11 <- PC(toy(), "epithelial", pixel_size_um = 0.5, random_seed = 2)
ck("E11 la scala si applica: area delle regioni = n_px * 0.25 µm²", isTRUE(all.equal(sum(E11$region_df$area_um2), nrow(toy()) * 0.25)), sum(E11$region_df$area_um2))
# BL-042 passo dei pixel
cs <- expand.grid(x = seq(1, 399, 2), y = seq(1, 299, 2)); cs$intensity_cluster <- 1L
e12 <- safe(PC(cs, "epithelial", pixel_size_um = 1))
ck("E12 BL-042: passo 2 rilevato e non dichiarato -> errore", !isok(e12) && grepl("BL-042", msg(e12)), msg(e12))
attr(cs, "stride") <- 2L; e13 <- safe(PC(cs, "epithelial", pixel_size_um = 1, random_seed = 2))
ck("E13 BL-042: passo dai metadati (attributo stride) -> regione dell'area giusta", isok(e13) && isTRUE(all.equal(sum(e13$region_df$area_um2), nrow(cs) * 4)), msg(e13))

nf <- sum(!unlist(res)); cat(sprintf("RESULT %d/%d PASS (mutante %s)\n", length(res) - nf, length(res), MUT))
quit(status = if (nf) 1 else 0)
