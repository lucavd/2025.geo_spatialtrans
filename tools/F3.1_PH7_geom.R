# tools/F3.1_PH7_geom.R - F3.1, indagine PH-7 (dichiarata il 2026-11-18): triangolazione fallita con la semina a compartimenti (BL-259).
# Copia di lavoro: .s27_cdt() e' sostituita SOLO in questa sessione R (il pacchetto non cambia). Fasi: (1) riproduzione con cattura della geometria che fallisce;
# (2) descrizione e prova delle correzioni candidate sulla sola geometria; (3) area A5_r3 rifatta con la correzione; (4) area di controllo A5_r1 (Q6) con e senza
# la correzione: conteggi identici. Uso: Rscript --vanilla tools/F3.1_PH7_geom.R <out_dir>
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
suppressPackageStartupMessages({ library(Matrix); library(sf); library(geos) })
for (f in sort(list.files("R", pattern = "^(04b|05b).*\\.R$", full.names = TRUE))) source(f)
out <- commandArgs(TRUE)[1]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
lay_dir <- "/mnt/micron/geo_spatialtrans/S1.6/layers"
all_tags <- sub("_U6\\.rds$", "", sort(list.files(lay_dir, pattern = "_U6\\.rds$")))
set.seed(20261114L, kind = "Mersenne-Twister"); seed_tab <- matrix(sample.int(.Machine$integer.max, length(all_tags) * 10), ncol = 10, dimnames = list(all_tags, NULL))
orig <- .s27_cdt
captured <- new.env()
.s27_cdt <- function(g, tol = 1e-4) {            # come l'originale, ma salva la geometria che fallisce prima di fermarsi
  r <- tryCatch(orig(g, tol), error = function(e) e)
  if (inherits(r, "error")) {
    k <- as.integer(sub(".*geometria ([0-9]+).*", "\\1", conditionMessage(r)))
    captured$g <- g[k]; captured$k <- k; captured$msg <- conditionMessage(r); stop(r)
  }
  r
}
L <- readRDS(file.path(lay_dir, "A5_r3_Q6.rds")); s <- seed_tab["A5_r3", 1L]
e1 <- tryCatch({ simulate_visium_hd_pieces(L, seed = s, frame_open = TRUE, workers = 1L, assemble = "full", keep_profiles = FALSE); "nessun errore" }, error = function(e) conditionMessage(e))
cat("(1) riproduzione:", e1, "\n")
if (is.null(captured$g)) { cat("geometria non catturata\n"); quit(status = 1) }
g <- captured$g; writeLines(geos_write_wkt(g), file.path(out, "PH7_geometry.wkt"))
xy <- wk::wk_coords(g); n <- nrow(xy)
ed <- sqrt(diff(xy$x)^2 + diff(xy$y)^2)[xy$ring_id[-1] == xy$ring_id[-n]]
cat("(2) geometria", captured$k, ": tipo", geos_type(g), "| valida", geos_is_valid(g), "|", geos_is_valid_detail(g)$reason, "| area", geos_area(g), "um2 | vertici", n,
    "| anelli", length(unique(xy$ring_id)), "| spigolo minimo", signif(min(ed), 3), "um | spigoli < 1e-3 um:", sum(ed < 1e-3), "\n")
try_fix <- function(name, h) {
  t <- tryCatch(geos_constrained_delaunay_triangles(h), error = function(e) NULL)
  a <- if (is.null(t)) NA else sum(geos_area(t))
  cat(sprintf("    %-38s triangolazione %-7s area triangoli %s (originale %.6g, scarto relativo %s)\n", name, if (is.null(t)) "FALLITA" else "ok",
              format(a, digits = 7), geos_area(g), if (is.na(a)) "-" else format(signif(a / geos_area(g) - 1, 3))))
  !is.null(t)
}
ok_a <- try_fix("make_valid", geos_make_valid(g))
ok_b <- try_fix("set_precision 1e-6", geos_set_precision(g, 1e-6))
ok_c <- try_fix("simplify 1e-3", geos_simplify_preserve_topology(g, 1e-3))
ok_d <- try_fix("make_valid + set_precision 1e-6", geos_set_precision(geos_make_valid(g), 1e-6))
ok_e <- try_fix("erosione 1e-6 um (dopo PH-7b)", geos_buffer(g, -1e-6))
fixer <- if (ok_e) function(h) geos::geos_buffer(h, -1e-6) else NULL   # scelta dopo la diagnosi PH-7b: un buco tocca il bordo esterno in un punto
if (is.null(fixer)) { cat("nessuna correzione candidata funziona\n"); quit(status = 0) }
n_fix <- 0L
.s27_cdt <- function(g, tol = 1e-4) {            # candidata: come l'originale, con un ultimo ripiego (fixer) dopo la semplificazione
  tri <- tryCatch(geos::geos_constrained_delaunay_triangles(g), error = function(e) NULL)
  if (!is.null(tri)) { attr(tri, "fallback") <- integer(0); return(tri) }
  outl <- vector("list", length(g)); fb <- integer(0)
  for (k in seq_along(g)) {
    t1 <- tryCatch(geos::geos_constrained_delaunay_triangles(g[k]), error = function(e) NULL)
    if (is.null(t1)) {
      t1 <- tryCatch(geos::geos_constrained_delaunay_triangles(geos::geos_simplify_preserve_topology(g[k], tol)), error = function(e) NULL)
      if (is.null(t1)) { t1 <- tryCatch(geos::geos_constrained_delaunay_triangles(fixer(g[k])), error = function(e) NULL); n_fix <<- n_fix + 1L }
      if (is.null(t1)) stop(sprintf(".s27_cdt(): triangolazione fallita anche dopo le correzioni (geometria %d).", k), call. = FALSE)
      fb <- c(fb, k)
    }
    outl[[k]] <- t1
  }
  tri <- do.call(c, outl); attr(tri, "fallback") <- fb; tri
}
o3 <- tryCatch(simulate_visium_hd_pieces(L, seed = s, frame_open = TRUE, workers = 1L, assemble = "full", keep_profiles = FALSE), error = function(e) conditionMessage(e))
cat("(3) A5_r3 con la correzione:", if (is.character(o3)) o3 else sprintf("ok, UMI totali %.0f, geometrie corrette con il nuovo ripiego %d", sum(o3$counts_2um), n_fix), "\n")
Lc <- readRDS(file.path(lay_dir, "A5_r1_Q6.rds")); sc <- seed_tab["A5_r1", 1L]
n_fix <- 0L; oc_fix <- simulate_visium_hd_pieces(Lc, seed = sc, frame_open = TRUE, workers = 1L, assemble = "full", keep_profiles = FALSE)
nf <- n_fix; .s27_cdt <- orig
oc_orig <- simulate_visium_hd_pieces(Lc, seed = sc, frame_open = TRUE, workers = 1L, assemble = "full", keep_profiles = FALSE)
cat("(4) controllo A5_r1 (Q6): conteggi a 8 um identici con e senza la correzione:", identical(oc_fix$counts_2um, oc_orig$counts_2um), "| nuovo ripiego usato", nf, "volte\n")
