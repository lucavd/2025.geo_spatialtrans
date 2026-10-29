# R/testing/test_S1.6.R - S1.6: test del consolidamento del codice dello Step 1 (cancello G-1; mutanti G-2).
# Pre-registrazione: docs/steps/S1.6.md. Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S16_RSRC=<cartella con R/>]
#   Rscript --vanilla R/testing/test_S1.6.R      I mutanti sono copie modificate del codice (R/testing/mutants_S1.6.py, run_S1.6.sh mutants).
.libPaths(Sys.getenv("S16_LIB", "renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"))   # S16_LIB: libreria del restore da zero (G-4)
suppressPackageStartupMessages({ library(sf); library(geos) })
SRC <- Sys.getenv("S16_RSRC", ".")
for (f in sort(list.files(file.path(SRC, "R"), pattern = "^04b.*\\.R$", full.names = TRUE))) source(f)
res <- list()
ck <- function(name, ok, val = "") { ok <- isTRUE(ok); res[[length(res) + 1]] <<- ok
  cat(sprintf("%s %s %s\n", if (ok) "PASS" else "FAIL", name, paste(val, collapse = " "))) }
safe <- function(expr) tryCatch(expr, error = function(e) structure(list(err = conditionMessage(e)), class = "t_err"))
PRE <- c(A1 = "epithelial", A2 = "tumor_microenv", A3 = "stromal_rich", A4 = "lymphoid", A5 = "brain_cortex", A6 = "columnar_muscle")
square <- function(side) list(region_df = data.frame(region_id = 1L, cluster_id = 1L, area_um2 = side^2),
  region_polygons = sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(side, 0), c(side, side), c(0, side), c(0, 0))))))
spec_c <- function(A) tissue_preset_spec(PRE[[A]], size_arm = "group", spatial = "compartments", min_dist_mode = "nucleus")
comp_of <- function(sp) data.frame(cluster_id = 1L, cell_type = sp$composition$cell_type, fraction = sp$composition$fraction)
H <- function(pv = 0, lv = NA, sd = 0, ld = NA, cap = NA) { h <- data.frame(void_frac = pv, void_ell_um = lv, dens_sigma = sd, dens_ell_um = ld)
  if (!is.na(cap)) h$dens_cap <- cap; h }
seed_c <- function(A, side, seed, het, ...) { sp <- spec_c(A)
  seed_centroids(square(side), sp$cell_types, comp_of(sp), random_seed = seed, spatial = "compartments", compartments = sp$compartments,
                 heterogeneity = het, verbose = FALSE, ...) }
clust_sq <- function(side) { g <- expand.grid(x = seq_len(side), y = seq_len(side)); g$intensity_cluster <- 1L; g }
side_of <- function(A) if (A %in% c("A5", "A6")) 600 else 300
# ---- T1: predefiniti senza quote (BL-170, DEC-025) ----
ck("T1a seed_centroids e place_cells: quote spente per predefinito (senza opzioni di sessione)",
   is.null(getOption("geo.spatialtrans.quotas")) && identical(eval(formals(seed_centroids)$quotas), FALSE) && identical(eval(formals(place_cells)$quotas), FALSE))
q0 <- lapply(names(PRE), function(A) seed_c(A, side_of(A), 10, NULL)); q1 <- seed_c("A4", 400, 10, NULL, quotas = TRUE)
rej0 <- vapply(q0, function(s) s$info$compartments$n_quota_reject, 0); qon <- vapply(q0, function(s) isTRUE(s$info$compartments$quotas), TRUE)
ck("T1b con il predefinito 0 rifiuti per quota nei 6 preset; con quotas = TRUE rifiuti > 0 (linfonodo, 400 um)",
   all(rej0 == 0) && !any(qon) && q1$info$compartments$n_quota_reject > 0,
   sprintf("rifiuti %s / %d", paste(rej0, collapse = ","), q1$info$compartments$n_quota_reject))
ck("T1c quotas deve essere TRUE o FALSE (\"TRUE\" o NA rifiutati)", inherits(safe(seed_c("A2", 200, 1, NULL, quotas = "TRUE")), "t_err") &&
   inherits(safe(seed_c("A2", 200, 1, NULL, quotas = NA)), "t_err"))
# ---- T2: limite esatto del diradamento (BL-161) ----
t2 <- do.call(rbind, lapply(names(PRE), function(A) do.call(rbind, lapply(1:3, function(s) {
  sc <- seed_c(A, side_of(A), 7000 + s, H(0.1, 16)); ci <- sc$info$compartments; r <- ci$regions[[1]]
  rho <- sc$type_df$density
  data.frame(A = A, seed = s, exceed = ci$n_exceed_wmax, bound = r$bound, wmax = r$wmax, rmax = max(rho), wgrid = r$w_grid_max, failed = sum(sc$info$n_failed)) }))))
ck("T2a limite esatto: 0 superamenti, limite = max rho_t (1 + 1e-12) e >= massimo sulla griglia (6 preset x 3 semi, vuoti 0.1 / 16 um)",
   all(t2$exceed == 0) && all(t2$bound == "exact") && all(abs(t2$wmax / (t2$rmax * (1 + 1e-12)) - 1) < 1e-14) && all(t2$wgrid <= t2$wmax),
   sprintf("superamenti %d, rapporto limite/griglia %.3f-%.3f", sum(t2$exceed), min(t2$wmax / t2$wgrid), max(t2$wmax / t2$wgrid)))
d1 <- seed_c("A2", 300, 71, H(0.1, 16, 0.3, 40)); d2 <- seed_c("A2", 300, 71, H(0.1, 16, 0.3, 40, cap = TRUE))
r1 <- d1$info$compartments$regions[[1]]; r2 <- d2$info$compartments$regions[[1]]; rmx <- max(d2$type_df$density)
ck("T2b campo della densita' senza tetto: ripiego sulla griglia registrato; con il tetto: limite esatto max rho_t x H",
   r1$bound == "grid" && abs(r1$wmax / (1.1 * r1$w_grid_max) - 1) < 1e-12 && r2$bound == "exact" &&
     abs(r2$wmax / (rmx * (if (is.na(r2$cap$s)) 1 else r2$cap$H) * (1 + 1e-12)) - 1) < 1e-12 && d2$info$compartments$n_exceed_wmax == 0,
   sprintf("%s / %s", r1$bound, r2$bound))
g0 <- seed_c("A2", 300, 72, H(0.1, 16), thinning_bound = "grid")$info$compartments$regions[[1]]
ck("T2c thinning_bound = \"grid\": 1.1 volte il massimo sulla griglia (comportamento di S1.5c-S1.5e)", g0$bound == "grid" && abs(g0$wmax / (1.1 * g0$w_grid_max) - 1) < 1e-12)
cnt <- function(b) t(vapply(1:100, function(s) { sc <- seed_c("A2", 250, 8000 + s, H(0.1, 16), thinning_bound = b)
  as.numeric(table(factor(sc$centroids$cell_type, levels = sc$type_df$cell_type))) }, numeric(2)))
ce <- cnt("exact"); cg <- cnt("grid"); z <- (colMeans(ce) - colMeans(cg)) / sqrt(apply(ce, 2, var) / 100 + apply(cg, 2, var) / 100)
ck("T2d stessa legge: conteggi medi per tipo con limite esatto e con la griglia entro 3 errori Monte Carlo (100 semi, tumore del colon)",
   all(abs(z) < 3), sprintf("z %.2f %.2f", z[1], z[2]))
# T2e (revisione RC-03): senza eterogeneita' limite esatto; massimo indipendente di w su 2e5 punti <= limite
t2e <- vapply(names(PRE), function(A) { sp <- spec_c(A); sc <- seed_c(A, side_of(A), 7100, NULL); r <- sc$info$compartments$regions[[1]]
  f <- sp$composition$fraction[match(sc$type_df$cell_type, sp$composition$cell_type)]; f <- f / sum(f)
  s <- sp$compartments$s[match(sc$type_df$cell_type, sp$compartments$cell_type)]; s[is.na(s)] <- 0
  set.seed(7); x <- runif(2e5, 0, side_of(A)); y <- runif(2e5, 0, side_of(A))
  w <- .sc_prob(.sc_eta(r$field_types, s, log(f), x, y), r$a, sc$type_df$density)$w
  c(exact = r$bound == "exact", ok = max(w) <= r$wmax, ratio = max(w) / r$wmax, exceed = sc$info$compartments$n_exceed_wmax) }, numeric(4))
ck("T2e senza eterogeneita': limite esatto e massimo di w su 2e5 punti indipendenti <= limite (6 preset)",
   all(t2e["exact", ] == 1) && all(t2e["ok", ] == 1) && all(t2e["exceed", ] == 0), sprintf("%s %.3f", colnames(t2e), t2e["ratio", ]))
# T2f (revisione RC-03): due regioni con tipi diversi: il limite usa i tipi della regione, non il catalogo
two <- list(region_df = data.frame(region_id = 1:2, cluster_id = 1:2, area_um2 = c(4e4, 4e4)),
  region_polygons = sf::st_sfc(sf::st_polygon(list(rbind(c(0, 0), c(200, 0), c(200, 200), c(0, 200), c(0, 0)))),
                               sf::st_polygon(list(rbind(c(200, 0), c(400, 0), c(400, 200), c(200, 200), c(200, 0))))))
ct2 <- data.frame(cell_type = c("a", "b", "c"), density = c(4000, 2000, 12000), min_dist_um = 3)
cm2 <- data.frame(cluster_id = c(1L, 1L, 2L, 2L), cell_type = c("a", "b", "a", "c"), fraction = c(0.5, 0.5, 0.5, 0.5))
s2f <- seed_centroids(two, ct2, cm2, random_seed = 5, spatial = "compartments", compartments = data.frame(cell_type = c("a", "b", "c"), ell_um = 20, s = 2), verbose = FALSE)
w2 <- vapply(s2f$info$compartments$regions, function(r) r$wmax, 0)
ck("T2f limite per regione = max rho_t dei tipi della regione (4000 e 12000, non 12000 e 12000)",
   all(abs(w2 / (c(4000, 12000) * (1 + 1e-12)) - 1) < 1e-14) && s2f$info$compartments$n_exceed_wmax == 0, sprintf("%.6g", w2))
# ---- T3: opzioni di sessione per riprodurre gli step passati ----
op <- options(geo.spatialtrans.quotas = TRUE, geo.spatialtrans.thinning_bound = "grid")
a3 <- seed_c("A3", 300, 31, H(0.1, 16)); options(op); b3 <- seed_c("A3", 300, 31, H(0.1, 16), quotas = TRUE, thinning_bound = "grid")
ck("T3a options(quotas = TRUE, thinning_bound = \"grid\") = argomenti espliciti (bit per bit)", identical(a3$centroids, b3$centroids) &&
   a3$info$compartments$quotas && a3$info$compartments$regions[[1]]$bound == "grid")
op <- options(geo.spatialtrans.preset_types = "single")
s3 <- tissue_preset_spec("epithelial"); options(op); s3b <- tissue_preset_spec("epithelial", types = "single")
ck("T3b options(preset_types = \"single\") = types = \"single\" (preset a un tipo di S1.5)", identical(s3, s3b) && identical(s3$types, "single") &&
   identical(tissue_preset_spec("epithelial")$types, "multi"))
# T3c, T3d (revisione RC-02): place_cells() inoltra quotas, thinning_bound e legge l'opzione dei preset a un tipo
pcc <- function(...) suppressWarnings(place_cells(clust_sq(300), "lymphoid", pixel_size_um = 1, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
  spatial = "compartments", min_dist_mode = "nucleus", size_arm = "group", random_seed = 61, verbose = FALSE, strict = FALSE, ...))
ci_d <- pcc()$metadata$seed_info$compartments; ci_l <- pcc(quotas = TRUE, thinning_bound = "grid")$metadata$seed_info$compartments
op <- options(geo.spatialtrans.quotas = TRUE, geo.spatialtrans.thinning_bound = "grid"); ci_o <- pcc()$metadata$seed_info$compartments; options(op)
bd <- function(ci) unique(vapply(ci$regions, function(z) z$bound, ""))
ck("T3c place_cells(): predefiniti senza quote e limite esatto; argomenti e opzioni di sessione arrivano alla semina (linfonodo)",
   !ci_d$quotas && identical(bd(ci_d), "exact") && ci_l$quotas && identical(bd(ci_l), "grid") && ci_l$n_quota_reject > 0 &&
   ci_o$quotas && identical(bd(ci_o), "grid"), sprintf("%s/%s, %s/%s, %s/%s", ci_d$quotas, bd(ci_d), ci_l$quotas, bd(ci_l), ci_o$quotas, bd(ci_o)))
pcs <- function(...) suppressWarnings(place_cells(clust_sq(200), "epithelial", pixel_size_um = 1, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
  random_seed = 62, verbose = FALSE, strict = FALSE, ...))
op <- options(geo.spatialtrans.preset_types = "single"); Lo <- pcs(); options(op); Le <- pcs(types = "single")
lmd5 <- function(L) { f <- tempfile(); saveRDS(list(L$cell_df, sf::st_as_binary(L$cell_territories), sf::st_as_binary(L$cell_nuclei)), f, compress = FALSE)
  on.exit(unlink(f)); unname(tools::md5sum(f)) }
ck("T3d place_cells() con options(preset_types = \"single\") = types = \"single\" (md5), catalogo epithelial_all",
   identical(lmd5(Lo), lmd5(Le)) && identical(Lo$cell_types_resolved$cell_type, "epithelial_all"))
# ---- T4: C7 con la densita' della miscela quando la semina e' cieca al tipo (BL-156) ----
L4 <- suppressWarnings(place_cells(clust_sq(600), "brain_cortex", pixel_size_um = 1, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
                                   random_seed = 44, verbose = FALSE, strict = FALSE))
c7 <- function(v) v$table[v$table$check == "C7", ]
v_b <- c7(L4$metadata$validation)
L4x <- L4; L4x$cell_types_resolved$density <- L4x$cell_types_resolved$density * 4
v_bx <- c7(validate_cell_layer(L4x, type_order = "blind")); v_sx <- c7(validate_cell_layer(L4x, type_order = "sequential"))
L4m <- L4; L4m$region_df$target_density_weighted <- L4m$region_df$target_density_weighted * 4
v_bm <- c7(validate_cell_layer(L4m, type_order = "blind"))
ck("T4a semina cieca (cervello, preset per tipo): C7 PASS con la densita' della miscela; tipo puro x 4 ignorato; miscela x 4 -> WARN",
   identical(L4$metadata$type_order, "blind") && v_b$status == "PASS" && grepl("miscela", v_b$threshold) && v_bx$status == "PASS" && v_bm$status == "WARN",
   sprintf("%s | %s | %s", v_b$status, v_bx$status, v_bm$status))
ck("T4b semina in sequenza: C7 con la densita' del tipo (tipo x 4 -> WARN, design §6.1)", v_sx$status == "WARN" && !grepl("miscela", v_sx$threshold))
L4s <- suppressWarnings(place_cells(clust_sq(300), "epithelial", pixel_size_um = 1, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
                                    types = "single", random_seed = 45, verbose = FALSE, strict = FALSE))
ck("T4c preset a un tipo: semina in sequenza e C7 con la densita' del tipo, come in S1.5", identical(L4s$metadata$type_order, "sequential") &&
   !grepl("miscela", c7(L4s$metadata$validation)$threshold))
# T4d (revisione RC-01, RC-16): compartimenti con la densita' del tipo; validate_cell_layer() senza type_order = validazione interna
Lc <- pcc(); vci <- c7(Lc$metadata$validation); vcn <- c7(validate_cell_layer(Lc)); vbn <- c7(validate_cell_layer(L4))
ck("T4d compartimenti: C7 con la densita' del tipo (anche senza type_order); semina cieca senza type_order = riga interna",
   !grepl("miscela", vci$threshold) && identical(vci$threshold, vcn$threshold) && identical(vci$value, vcn$value) &&
   identical(vbn$value, v_b$value) && identical(vbn$threshold, v_b$threshold) && identical(Lc$metadata$seed_info$type_order, "compartments"))
# T4e (revisione RC-18): semina cieca nei 6 preset: C7 PASS con la densita' della miscela
t4e <- vapply(names(PRE), function(A) { L <- suppressWarnings(place_cells(clust_sq(side_of(A)), PRE[[A]], pixel_size_um = 1, min_region_area_um2 = 0,
  simplify_tol_um = 0, stride = 1L, random_seed = 47, verbose = FALSE, strict = FALSE)); v <- c7(L$metadata$validation)
  paste(v$status, grepl("miscela", v$threshold)) }, "")
ck("T4e semina cieca, 6 preset: C7 PASS con la miscela", all(t4e == "PASS TRUE"), t4e)
w4f <- character(0)
invisible(withCallingHandlers(place_cells(clust_sq(300), "lymphoid", pixel_size_um = 1, min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L,
  spatial = "compartments", min_dist_mode = "nucleus", size_arm = "group", random_seed = 61, verbose = FALSE, strict = FALSE, type_order = "blind"),
  warning = function(w) { w4f <<- c(w4f, conditionMessage(w)); invokeRestart("muffleWarning") }))
ck("T4f place_cells(): type_order con i compartimenti e' ignorato con un avviso", any(grepl("type_order", w4f)))
# ---- T5: type_order per i cataloghi dell'utente (BL-155) ----
ctu <- data.frame(cell_type = c("a", "b"), density = c(9000, 3000), min_dist_um = 4, nucleus_shape = "circle", nucleus_radius = NA_real_,
                  nucleus_to_eq_ratio = 0.4, nucleus_area_slope = 0.2, nucleus_area_sdlog = 0.3, nucleus_aspect = 1, nucleus_orient = "random",
                  nucleus_nc_max = 0.8, stringsAsFactors = FALSE)
cmu <- data.frame(cluster_id = 1L, cell_type = c("a", "b"), fraction = c(0.7, 0.3))
pcu <- function(ct, cm, ord) suppressWarnings(place_cells(clust_sq(200), cell_types = ct, region_composition = cm, pixel_size_um = 1,
  min_region_area_um2 = 0, simplify_tol_um = 0, stride = 1L, random_seed = 46, type_order = ord, verbose = FALSE, strict = FALSE))
geo_md5 <- function(L) { f <- tempfile(); saveRDS(sf::st_as_binary(L$cell_territories), f, compress = FALSE); on.exit(unlink(f)); unname(tools::md5sum(f)) }
b1 <- pcu(ctu, cmu, "blind"); b2 <- pcu(ctu[2:1, ], cmu[2:1, ], "blind"); s1 <- pcu(ctu, cmu, "sequential"); s2 <- pcu(ctu[2:1, ], cmu[2:1, ], "sequential")
dflt <- pcu(ctu, cmu, NULL)
ck("T5a type_order = \"blind\": territori identici con il catalogo permutato; \"sequential\" no (controllo)",
   identical(geo_md5(b1), geo_md5(b2)) && !identical(geo_md5(s1), geo_md5(s2)) && identical(b1$metadata$type_order, "blind"))
ck("T5b predefinito per i cataloghi dell'utente invariato (in sequenza)", identical(dflt$metadata$type_order, "sequential") && identical(geo_md5(dflt), geo_md5(s1)))
# ---- T6: nessun interruttore dei mutanti nel codice di produzione (BL-155) ----
src <- unlist(lapply(list.files(file.path(SRC, "R"), pattern = "\\.R$", full.names = TRUE), readLines))
ck("T6 nessun Sys.getenv(\"S1...\") in R/*.R", !any(grepl("Sys.getenv\\(\"S1", src)), sum(grepl("Sys.getenv\\(\"S1", src)))
# ---- T7: i test di S1.1-S1.4 scrivono solo nella cartella dei log (BL-165; revisione RC-04, RC-05) ----
ff <- function(f) if (file.exists(file.path(SRC, f))) file.path(SRC, f) else f        # i mutanti possono toccare anche i test (N6, N16)
wr <- "write\\.csv|saveRDS|writeLines|write_parquet|png\\(|pdf\\(|ggsave|file\\.copy"
tw <- unlist(lapply(sprintf("R/testing/test_S1.%d.R", 1:4), function(f) { x <- readLines(ff(f)); x[grepl(wr, x) & grepl("\\bRES\\b|results/", x)] }))
ok7 <- all(vapply(sprintf("R/testing/test_S1.%d.R", 1:4), function(f) any(grepl("GEO_TEST_OUT", readLines(ff(f)))), TRUE))
o11 <- readLines(ff("R/testing/test_S1.1.R")); okd <- any(grepl('^OUTD <- "/mnt/micron/geo_spatialtrans/S1.1/outputs"', o11)) && !any(grepl("^OUTD <- Sys.getenv", o11))
ck("T7 test_S1.1-S1.4: nessuna scrittura in results/ (le tabelle vanno in GEO_TEST_OUT); le RDS di S1.1 restano in /mnt/micron", !length(tw) && ok7 && okd, tw)
nf <- sum(!unlist(res)); cat(sprintf("RESULT %d/%d PASS\n", length(res) - nf, length(res)))
