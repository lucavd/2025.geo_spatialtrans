# R/testing/test_S1.4d.R — S1.4d: test del modello "clip" di derive_cell_geometry() (cancello G-1; mutanti G-2).
# Pre-registrazione: docs/steps/S1.4d.md (commit 5118bce). Uso (root del repo):
#   R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent [S14D_MUTANT=none|M1..M7] Rscript --vanilla R/testing/test_S1.4d.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(sf); library(geos) })
source("tools/S1.4d_mutants.R"); source("tools/S1.4c_shape.R")          # S1.4c_shape: prototipo per U11
MUT <- Sys.getenv("S14D_MUTANT", "none"); E <- dcg_mutate_d(dcg_env_d(), MUT)
res <- list()
ck <- function(name, ok, val = "") { ok <- isTRUE(ok); res[[length(res) + 1]] <<- ok
  cat(sprintf("%s %s %s\n", if (ok) "PASS" else "FAIL", name, paste(val, collapse = " "))) }
safe <- function(expr) tryCatch(expr, error = function(e) structure(list(err = conditionMessage(e)), class = "dcg_err"))
isok <- function(o) !inherits(o, "dcg_err")
msg <- function(o) if (isok(o)) "" else o$err
DCG <- function(...) E$derive_cell_geometry(..., verbose = FALSE)
sqp <- function(x0, y0, s) st_polygon(list(rbind(c(x0, y0), c(x0 + s, y0), c(x0 + s, y0 + s), c(x0, y0 + s), c(x0, y0))))
mk <- function(polys, x, y, type = "c") { g <- st_sfc(polys); n <- length(g)
  list(cen = data.frame(cell_id = seq_len(n), region_id = 1L, x = x, y = y, cell_type = rep(type, length.out = n), stringsAsFactors = FALSE),
       tess = list(cell_territories = g, territory_df = data.frame(cell_id = seq_len(n), territory_area = as.numeric(st_area(g))))) }
CT <- function(A, ...) data.frame(cell_type = "c", density = 1000, nucleus_radius = sqrt(A / pi), ..., stringsAsFactors = FALSE)
rel <- function(a, b) abs(a - b) / b
inside <- function(o, x) all(diag(st_within(o$cell_nuclei, x$tess$cell_territories, sparse = FALSE)))

# U1 nucleo che non tocca il bordo: area = bersaglio, non ritagliato, 0 iterazioni
x1 <- mk(list(sqp(0, 0, 10)), 5, 5); o <- safe(DCG(x1$cen, x1$tess, CT(50)))
ck("U1 senza ritaglio: area = bersaglio, 0 iterazioni", isok(o) && rel(o$cell_df$nucleus_area, 50) <= 1e-6 && !o$cell_df$nucleus_clipped &&
   o$cell_df$nucleus_n_iter == 0, msg(o))
# U2 nucleo piu' grande del disco inscritto: ritagliato e ingrandito fino all'area bersaglio, dentro
o <- safe(DCG(x1$cen, x1$tess, CT(90)))
ck("U2 ritagliato: area = bersaglio (rel <= 1e-6), dentro", isok(o) && o$cell_df$nucleus_clipped && rel(o$cell_df$nucleus_area, 90) <= 1e-6 && inside(o, x1),
   if (isok(o)) signif(rel(o$cell_df$nucleus_area, 90), 3) else msg(o))
# U3 tetto phi = 0.95, bersaglio 120 > A_T = 100: area = phi * area del territorio rientrato
o <- safe(DCG(x1$cen, x1$tess, CT(120, nucleus_nc_max = 0.95)))
ck("U3 tetto: capped, area = phi * area rientrata", isok(o) && o$cell_df$nucleus_capped && abs(o$cell_df$nucleus_area_cap - 95) < 4e-4 &&
   rel(o$cell_df$nucleus_area, o$cell_df$nucleus_area_cap) <= 1e-6, if (isok(o)) sprintf("cap=%.6f area=%.6f", o$cell_df$nucleus_area_cap, o$cell_df$nucleus_area) else msg(o))
# U4 due quadrati adiacenti con nuclei grandi: si toccano senza sovrapporsi, ciascuno nel proprio territorio
x4 <- mk(list(sqp(0, 0, 10), sqp(10, 0, 10)), c(5, 15), c(5, 5)); o <- safe(DCG(x4$cen, x4$tess, CT(95)))
if (isok(o)) { g <- as_geos_geometry(o$cell_nuclei); ov <- geos_area(geos_intersection(g[1], g[2])); dd <- geos_distance(g[1], g[2]) } else { ov <- NA; dd <- NA }
ck("U4 nessuna sovrapposizione (<= 1e-9 µm²)", isok(o) && ov <= 1e-9, signif(ov, 3))
ck("U4 contatto: distanza in (0, 3e-5] (fessura 2 x rientro)", isok(o) && dd > 0 && dd <= 3e-5, signif(dd, 3))
ck("U4 dentro i propri territori, area = bersaglio", isok(o) && inside(o, x4) && all(rel(o$cell_df$nucleus_area, 95) <= 1e-6), msg(o))
# U5 phi = 1 e bersaglio = A_T: area = territorio rientrato, senza errore
o <- safe(DCG(x1$cen, x1$tess, CT(100)))
ck("U5 phi = 1, bersaglio = A_T: area = territorio rientrato", isok(o) && o$cell_df$nucleus_area_cap < 100 && o$cell_df$nucleus_area_cap > 100 - 4e-4 &&
   rel(o$cell_df$nucleus_area, o$cell_df$nucleus_area_cap) <= 1e-6, msg(o))
# U6 territorio non convesso (L)
xL <- mk(list(st_polygon(list(rbind(c(0, 0), c(20, 0), c(20, 10), c(10, 10), c(10, 20), c(0, 20), c(0, 0))))), 5, 5)
o <- safe(DCG(xL$cen, xL$tess, CT(250)))
ck("U6 territorio non convesso: valido, dentro, area = bersaglio", isok(o) && all(st_is_valid(o$cell_nuclei)) && inside(o, xL) && rel(o$cell_df$nucleus_area, 250) <= 1e-6, msg(o))
# U7 territorio con un buco
hole <- rbind(c(8, 8), c(8, 12), c(12, 12), c(12, 8), c(8, 8))
xH <- mk(list(st_polygon(list(rbind(c(0, 0), c(20, 0), c(20, 20), c(0, 20), c(0, 0)), hole))), 4, 4)
o <- safe(DCG(xH$cen, xH$tess, CT(300)))
ck("U7 territorio con buco: valido, dentro, area = bersaglio", isok(o) && all(st_is_valid(o$cell_nuclei)) && inside(o, xH) && rel(o$cell_df$nucleus_area, 300) <= 1e-6, msg(o))

# U8 Voronoi sintetico (300 punti in 50 x 50, seme 3; con il rientro calcolato sulla griglia 13 territori diventano vuoti,
# trappola BL-118b): nuclei grandi (rho 0.9, sigma 0.3), invarianti C3 / nessuna sovrapposizione / area / nessun vuoto
set.seed(3); P <- cbind(runif(300, 0, 50), runif(300, 0, 50))
box <- st_as_sfc(st_bbox(c(xmin = 0, ymin = 0, xmax = 50, ymax = 50)))
vv <- st_collection_extract(st_voronoi(st_sfc(st_multipoint(P)), box, point_order = TRUE), "POLYGON")
vv <- st_intersection(vv, box); stopifnot(length(vv) == 300, all(diag(st_intersects(st_as_sf(as.data.frame(P), coords = 1:2), vv, sparse = FALSE))))
x8 <- mk(as.list(vv), P[, 1], P[, 2])
ct8 <- data.frame(cell_type = "c", density = 300 / 2500 * 1e6, nucleus_to_eq_ratio = 0.9, nucleus_area_sdlog = 0.3)
o8 <- safe(DCG(x8$cen, x8$tess, ct8, random_seed = 11))
ck("U8 Voronoi sintetico: eseguito senza errori", isok(o8), msg(o8))
if (isok(o8)) {
  g8 <- as_geos_geometry(o8$cell_nuclei); T8 <- as_geos_geometry(x8$tess$cell_territories)
  ck("U8 C3 (geos): 300/300 dentro", all(geos_within(g8, T8)), sum(!geos_within(g8, T8)))
  ck("U8 C3 (sf, indipendente)", inside(o8, x8))
  tr <- geos_strtree(g8); q <- geos_strtree_query(tr, g8); pi_ <- rep(seq_along(q), lengths(q)); pj <- unlist(q); k <- pi_ < pj
  ov <- geos_area(geos_intersection(g8[pi_[k]], g8[pj[k]]))
  ck("U8 nessuna coppia sovrapposta (> 1e-9 µm²)", sum(ov > 1e-9) == 0, sprintf("max %.2e", max(c(0, ov))))
  ck("U8 area = A* (rel <= 1e-6), nessun vuoto", max(rel(o8$cell_df$nucleus_area, o8$cell_df$nucleus_area_cap)) <= 1e-6 && !any(st_is_empty(o8$cell_nuclei)),
     signif(max(rel(o8$cell_df$nucleus_area, o8$cell_df$nucleus_area_cap)), 3))
  ck("U8 il test esercita il ritaglio (> 50 ritagliati, >= 1 al tetto)", sum(o8$cell_df$nucleus_clipped) > 50 && sum(o8$cell_df$nucleus_capped) >= 1,
     sprintf("ritagliati %d, al tetto %d", sum(o8$cell_df$nucleus_clipped), sum(o8$cell_df$nucleus_capped)))
  ck("U16 iterazioni: massimo <= 40 (la bisezione del prototipo ne usava fino a 53)", max(o8$cell_df$nucleus_n_iter) <= 40, max(o8$cell_df$nucleus_n_iter))
  ck("U13 sf via WKB: st_area = nucleus_area (rel <= 1e-9)", max(rel(as.numeric(st_area(o8$cell_nuclei)), o8$cell_df$nucleus_area)) <= 1e-9)
  # U11 equivalenza con il prototipo S1.4c (kclip_nuclei, K2, phi = 1) sugli stessi bersagli
  Tg <- geos_set_precision(as_geos_geometry(x8$tess$cell_territories), S14C_GRID)
  pr <- kclip_nuclei(P[, 1], P[, 2], o8$cell_df$nucleus_area_target, Tg, "K2", a_t = x8$tess$territory_df$territory_area, phi = 1)
  gp <- geos_set_precision(pr$nuc, 0); a8 <- o8$cell_df$nucleus_area_cap
  sdA <- geos_area(geos_sym_difference(g8, gp)); dA <- abs(o8$cell_df$nucleus_area - pr$df$area)
  bound <- dA + sqrt(2) * S14C_GRID * pmax(geos_length(g8), geos_length(gp))   # rumore della griglia: strisce <= sqrt(2) g lungo il bordo
  ck("U11 = prototipo S1.4c: |dA| / A* <= 1e-6, differenza simmetrica entro il rumore della griglia, stesse celle al tetto",
     max(dA / a8) <= 1e-6 && all(sdA <= bound + 1e-12) && identical(pr$df$capped, o8$cell_df$nucleus_capped),
     sprintf("max|dA|/A* %.2e; sym/A* > 1e-6 in %d celle (forma pre-registrata di G-4a)", max(dA / a8), sum(sdA / a8 > 1e-6)))
  # U17 determinismo
  o8b <- DCG(x8$cen, x8$tess, ct8, random_seed = 11)
  ck("U17 determinismo: stesso seme, cell_df e WKB identici", identical(o8$cell_df, o8b$cell_df) &&
     identical(geos_write_wkb(as_geos_geometry(o8$cell_nuclei)), geos_write_wkb(as_geos_geometry(o8b$cell_nuclei))))
  # U12 flusso casuale: bersagli = rho^2 A_ref (A_T / A_ref) exp(sigma z), z = rnorm del seme; uguali al modo "scale"
  set.seed(11); z <- rnorm(300); aref <- 1e6 / ct8$density
  exp_t <- 0.9^2 * aref * (x8$tess$territory_df$territory_area / aref) * exp(0.3 * z)
  os <- safe(DCG(x8$cen, x8$tess, ct8, random_seed = 11, nucleus_model = "scale"))
  ck("U12 bersagli dal flusso dichiarato (z prima di u) e uguali al modo scale", max(rel(o8$cell_df$nucleus_area_target, exp_t)) <= 1e-12 && isok(os) &&
     identical(os$cell_df$nucleus_area_target, o8$cell_df$nucleus_area_target), signif(max(rel(o8$cell_df$nucleus_area_target, exp_t)), 3))
}
# U9 ellisse (aspetto 2.5, angolo fisso quasi lungo y) in un rettangolo stretto: ritagliata, area = bersaglio
xR <- mk(list(st_polygon(list(rbind(c(0, 0), c(6, 0), c(6, 20), c(0, 20), c(0, 0))))), 3, 10)
o <- safe(DCG(xR$cen, xR$tess, CT(80, nucleus_aspect = 2.5, nucleus_orient = "1.4707963")))
ck("U9 ellisse: ritagliata, area = bersaglio, dentro, semiassi in rapporto 2.5", isok(o) && o$cell_df$nucleus_clipped && rel(o$cell_df$nucleus_area, 80) <= 1e-6 &&
   inside(o, xR) && abs(o$cell_df$nucleus_major_r / o$cell_df$nucleus_minor_r - 2.5) < 1e-12, msg(o))
# U10 tipi multipli con phi diversi
ctm <- data.frame(cell_type = c("a", "b"), density = 1000, nucleus_radius = sqrt(80 / pi), nucleus_nc_max = c(0.5, 1))
o <- safe(DCG(transform(x4$cen, cell_type = c("a", "b")), x4$tess, ctm))
ck("U10 tipi multipli: a al tetto 0.5 x area rientrata, b all'area bersaglio", isok(o) && o$cell_df$nucleus_capped[1] && !o$cell_df$nucleus_capped[2] &&
   abs(o$cell_df$nucleus_area[1] - 0.5 * (10 - 2e-5)^2) < 1e-4 && rel(o$cell_df$nucleus_area[2], 80) <= 1e-6, msg(o))
# U13b tipo della colonna geometrica
o <- safe(DCG(x1$cen, x1$tess, CT(90)))
ck("U13b cell_nuclei sfc_POLYGON, area sf = nucleus_area", isok(o) && inherits(o$cell_nuclei, "sfc_POLYGON") && rel(as.numeric(st_area(o$cell_nuclei)), o$cell_df$nucleus_area) <= 1e-9, msg(o))
# U13c helper di conversione: un vettore geos con modello di precisione (geos_set_precision) va convertito con l'area giusta;
# la conversione diretta st_as_sfc() lo legge come fattore di scala (area 0). Sul percorso di derive_cell_geometry() i nuclei
# escono da geos_intersection_prec() senza modello di precisione: questo test e' l'unico che vede il mutante M6.
gq <- geos_set_precision(geos_make_polygon(c(0.1234567, 10.3, 10.3, 0.1234567), c(0.2, 0.2, 10.7, 10.7)), 1e-6)
sq_sf <- safe(E$.dcg_to_sfc(gq))
ck("U13c .dcg_to_sfc con modello di precisione: area = area geos", isok(sq_sf) && rel(as.numeric(st_area(sq_sf)), geos_area(gq)) <= 1e-9,
   if (isok(sq_sf)) as.numeric(st_area(sq_sf)) else msg(sq_sf))
# U14 controllo degli argomenti
er <- function(e) !isok(safe(e))
ck("U14a nucleus_nc_max = 0 rifiutato", er(DCG(x1$cen, x1$tess, CT(50, nucleus_nc_max = 0))))
ck("U14b nucleus_nc_max = 1.2 rifiutato", er(DCG(x1$cen, x1$tess, CT(50, nucleus_nc_max = 1.2))))
ck("U14c clip_inset_um = 0 rifiutato", er(DCG(x1$cen, x1$tess, CT(50), clip_inset_um = 0)))
ck("U14d area_rel_tol = 1e-2 rifiutato", er(DCG(x1$cen, x1$tess, CT(50), area_rel_tol = 1e-2)))
ck("U14e nucleus_model sconosciuto rifiutato", er(DCG(x1$cen, x1$tess, CT(50), nucleus_model = "foo")))
# U15 schema e centro
o <- safe(DCG(xL$cen, xL$tess, CT(250)))
cols <- c("cell_id", "region_id", "cluster_id", "cell_type", "x", "y", "territory_area", "eq_radius", "nucleus_shape", "nucleus_radius",
          "nucleus_major_r", "nucleus_minor_r", "nucleus_orient", "nucleus_cx", "nucleus_cy", "nucleus_area", "cytoplasm_area",
          "nucleus_area_target", "nucleus_area_cap", "nucleus_capped", "nucleus_clipped", "nucleus_scale", "nucleus_n_iter", "center_dist")
ck("U15a schema di cell_df (clip) e centro = generatore", isok(o) && identical(names(o$cell_df), cols) && o$cell_df$nucleus_cx == 5 && o$cell_df$nucleus_cy == 5 &&
   o$info$nucleus_model == "clip" && o$info$nucleus_center == "generator", msg(o))
o <- safe(DCG(xL$cen, xL$tess, CT(250), nucleus_center = "poi"))
ck("U15b clip con centro POI: dentro, area = bersaglio", isok(o) && inside(o, xL) && rel(o$cell_df$nucleus_area, 250) <= 1e-6 && o$cell_df$nucleus_cx != 5, msg(o))
nf <- sum(!unlist(res)); cat(sprintf("RESULT %d/%d PASS (mutante %s)\n", length(res) - nf, length(res), MUT))
quit(status = if (nf) 1 else 0)
