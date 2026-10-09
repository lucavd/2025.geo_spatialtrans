# results/S1.3/review/code/arbiter_check.R — revisione del codice S1.3 (RA-code), compito 4: arbitro hp_cell()
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.3/review/code/arbiter_check.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
args <- c("none", "/mnt/micron/geo_spatialtrans/S1.3"); commandArgs <- function(trailingOnly = TRUE) args
source("tools/S1.3_roi.R")                       # carica anche tools/S1.3_arbiter.R (hp_cell, arbitrate)
D <- "/mnt/micron/geo_spatialtrans/S1.3/review_code_data"; OUTD <- "results/S1.3/review/code"
f17 <- function(v) sprintf("%.17g", v)
rows <- list(); rec <- function(case, cell, a_hp, a_ref, ns_hp, ns_ref, K, bounded, note = "") {
  rows[[length(rows) + 1]] <<- data.frame(case = case, cell = cell, a_hp = f17(a_hp), a_ref = f17(a_ref),
    rel_err = f17(abs(a_hp - a_ref) / a_ref), ns_hp = ns_hp, ns_ref = ns_ref, K = K, bounded = bounded, note = note)
}
# cella esatta di riferimento: GEOS st_voronoi sull'insieme completo (envelope ampio)
geos_cell <- function(i, x, y) {
  env <- st_polygon(list(rbind(c(-1e5, -1e5), c(1e5, -1e5), c(1e5, 1e5), c(-1e5, 1e5), c(-1e5, -1e5))))
  v <- st_collection_extract(st_sfc(st_voronoi(st_multipoint(cbind(x, y)), envelope = env, point_order = TRUE)), "POLYGON")
  v[[i]]
}
nsd <- function(p) { m <- unclass(p)[[1]]; nrow(m) - 1L }
# 1) esagonale (a = 10, 15 x 15), cella centrale
a <- 10; g <- expand.grid(i = 0:14, j = 0:14); x <- g$i * a + (g$j %% 2) * a / 2; y <- g$j * a * sqrt(3) / 2
i <- which(g$i == 7 & g$j == 7); h <- hp_cell(i, x, y)
rec("esagonale a=10", i, h[["area"]], sqrt(3) / 2 * a^2, h[["nsides"]], 6L, h[["K"]], h[["bounded"]], "atteso analitico")
# 2) quadrato (passo 10, 15 x 15), cella centrale (vertici di grado 4: 4 punti cocircolari)
g <- expand.grid(i = 0:14, j = 0:14); x <- 10 * g$i; y <- 10 * g$j
i <- which(g$i == 7 & g$j == 7); h <- hp_cell(i, x, y)
rec("quadrato a=10", i, h[["area"]], 100, h[["nsides"]], 4L, h[["K"]], h[["bounded"]], sprintf("min_edge %.3g", h[["min_edge"]]))
# 2b) quadrato perturbato di 1e-9 (vertici quasi degeneri)
set.seed(3); xp <- x + runif(length(x), -1e-9, 1e-9); yp <- y + runif(length(y), -1e-9, 1e-9); h <- hp_cell(i, xp, yp); gc <- geos_cell(i, xp, yp)
rec("quadrato perturbato 1e-9", i, h[["area"]], as.numeric(st_area(st_sfc(gc))), h[["nsides"]], nsd(gc), h[["K"]], h[["bounded"]], sprintf("min_edge %.3g; ns_ref = vertici GEOS", h[["min_edge"]]))
# 3) celle sul bordo di un buco circolare (raggio 60) in un reticolo esagonale 400 x 400
a <- 10; g <- expand.grid(i = 0:40, j = 0:46); x <- g$i * a + (g$j %% 2) * a / 2; y <- g$j * a * sqrt(3) / 2
keep <- (x - 200)^2 + (y - 200)^2 > 60^2; x <- x[keep]; y <- y[keep]
rim <- which(sqrt((x - 200)^2 + (y - 200)^2) < 60 + a)
for (i in rim) { h <- hp_cell(i, x, y); gc <- geos_cell(i, x, y)
  rec("buco r=60 (esagonale)", i, h[["area"]], as.numeric(st_area(st_sfc(gc))), h[["nsides"]], nsd(gc), h[["K"]], h[["bounded"]]) }
# 4) controesempio costruito alla regola di arresto: anello di 32 punti (r = 1, varco di ±60° verso +x),
#    32 punti d'ombra (r = 1.2, verso -x), P = (1.5, 0) che taglia la punta, 2000 punti lontani (r 20-40)
th <- seq(60, 300, length.out = 32) * pi / 180; ts <- seq(120, 240, length.out = 32) * pi / 180
set.seed(7); rf <- runif(2000, 20, 40); af <- runif(2000, 0, 2 * pi)
x <- c(0, cos(th), 1.2 * cos(ts), 1.5, rf * cos(af)); y <- c(0, sin(th), 1.2 * sin(ts), 0, rf * sin(af))
h <- hp_cell(1L, x, y); gc <- geos_cell(1L, x, y); hfull <- hp_cell(1L, x, y, K0 = length(x) - 1L)
rec("controesempio arresto", 1L, h[["area"]], as.numeric(st_area(st_sfc(gc))), h[["nsides"]], nsd(gc), h[["K"]], h[["bounded"]], "ref = GEOS")
rec("controesempio arresto, K0 = n-1", 1L, hfull[["area"]], as.numeric(st_area(st_sfc(gc))), hfull[["nsides"]], nsd(gc), hfull[["K"]], hfull[["bounded"]], "ref = GEOS")
tab <- do.call(rbind, rows); write.csv(tab, file.path(OUTD, "arbiter_known_cases.csv"), row.names = FALSE)
print(tab[!grepl("buco", tab$case), ]); bb <- tab[grepl("buco", tab$case), ]
cat("buco: n =", nrow(bb), " max rel err =", max(as.numeric(bb$rel_err)), " lati diversi =", sum(bb$ns_hp != bb$ns_ref), "\n")
# 5) hp_cell su TUTTE le celle con tile dentro il quadrato del ROI, ROI reali A3 r1 e A5 r1 (tile del pacchetto come riferimento)
sys_rows <- list()
for (cs in list(c("A3", "r1", "cellpose_rgb"), c("A5", "r1", "cellpose_rgb"))) {
  w <- r3_roi_window(cs[1], cs[2], rois)
  gen <- as.data.frame(read_parquet(file.path(R3, "real", sprintf("%s_%s_%s_gen.parquet", cs[1], cs[2], cs[3]))))
  s <- s13_tessellate(gen$x, gen$y, w, backend = "geos"); tl <- s$tv$tiles
  inside <- which(lengths(st_within(tl, w$frame)) > 0)
  hh <- do.call(rbind, mclapply(inside, function(i) hp_cell(i, gen$x, gen$y), mc.cores = 30L))
  ag <- as.numeric(st_area(tl[inside])); ng <- vapply(tl[inside], nsd, 1L)
  sys_rows[[length(sys_rows) + 1]] <- data.frame(case = paste(cs, collapse = "_"), cell = inside, a_geos = f17(ag), a_hp = f17(hh[, "area"]),
    rel = f17(abs(ag - hh[, "area"]) / hh[, "area"]), ns_geos = ng, ns_hp = hh[, "nsides"], K = hh[, "K"], bounded = hh[, "bounded"], min_edge = f17(hh[, "min_edge"]))
}
st <- do.call(rbind, sys_rows); write.csv(st, file.path(D, "arbiter_systematic.csv"), row.names = FALSE)
for (cc in unique(st$case)) { z <- st[st$case == cc, ]
  cat(cc, ": celle", nrow(z), " rel > 1e-9:", sum(as.numeric(z$rel) > 1e-9), " max rel:", max(as.numeric(z$rel)), " lati diversi:", sum(z$ns_geos != z$ns_hp),
      " non limitate:", sum(z$bounded == 0), " K max:", max(z$K), "\n") }
# 6) punti delle righe nulle di S1.3_c10a_arbitration.csv (replica 1) per il ricalcolo in Python
ar <- read.csv("results/S1.3/S1.3_c10a_arbitration.csv"); nu <- unique(ar[ar$set == "null", c("archetype", "roi_id", "model")])
for (k in seq_len(nrow(nu))) {
  A <- nu$archetype[k]; roi <- nu$roi_id[k]; i <- which(rois$archetype == A & rois$roi_id == roi); a <- ARCH[ARCH$archetype == A, ]
  w <- r3_roi_window(A, roi, rois); n_obs <- readRDS(file.path(R3, "real", sprintf("%s_%s_%s_roi.rds", A, roi, a$primary)))$n
  p <- gen_points(i, nu$model[k], 1L, w, n_obs, a)
  write.csv(data.frame(x = f17(p$x), y = f17(p$y)), file.path(D, sprintf("arb_null_%s_%s_%s.csv", A, roi, nu$model[k])), row.names = FALSE, quote = FALSE)
}
cat("fatto\n")
