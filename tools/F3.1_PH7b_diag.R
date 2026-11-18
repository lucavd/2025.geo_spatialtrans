# tools/F3.1_PH7b_diag.R - F3.1, PH-7 seguito: diagnosi della geometria 168 (A5_r3, Q6) salvata in PH7_geometry.wkt e correzioni ulteriori,
# provate dopo che le quattro pre-dichiarate sono fallite (non pre-registrate). Uso: Rscript --vanilla tools/F3.1_PH7b_diag.R <PH7_dir>
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu"); suppressPackageStartupMessages({ library(geos); library(wk) })
d <- commandArgs(TRUE)[1]; g <- geos_read_wkt(readLines(file.path(d, "PH7_geometry.wkt")))
cat("GEOS", as.character(geos_version()), "| valida:", geos_is_valid(g), "|", geos_is_valid_detail(g)$reason, "\n")
sh <- geos_ring_n(g, 1L); nh <- geos_num_interior_rings(g); cat("buchi:", nh, "\n")
for (i in seq_len(nh)) { h <- geos_ring_n(g, i + 1L); ii <- geos_intersection(sh, h)
  cat(sprintf("  buco %d: area %.3f um2, vertici %d, distanza dal bordo esterno %.4g um, contatto: %s\n", i, geos_area(geos_make_polygon(wk_coords(h)$x, wk_coords(h)$y)),
              geos_num_coordinates(h), geos_distance(sh, h), if (geos_is_empty(ii)) "nessuno" else geos_write_wkt(ii))) }
if (nh > 1) for (i in 1:(nh - 1)) for (j in (i + 1):nh) { dd <- geos_distance(geos_ring_n(g, i + 1L), geos_ring_n(g, j + 1L)); if (dd < 1e-3) cat(sprintf("  buchi %d e %d a distanza %.3g\n", i, j, dd)) }
xy <- wk_coords(g); key <- paste(round(xy$x, 9), round(xy$y, 9)); dup <- table(key); cat("vertici ripetuti (escluse chiusure degli anelli):", sum(dup > 1) - (nh + 1), "\n")
A0 <- geos_area(g)
try_fix <- function(name, h) {
  t <- tryCatch(geos_constrained_delaunay_triangles(h), error = function(e) conditionMessage(e))
  if (is.character(t)) { cat(sprintf("  %-34s FALLITA: %s\n", name, substr(t, 1, 120))); return(invisible(FALSE)) }
  a <- sum(geos_area(geos_geometry_n(t, seq_len(geos_num_geometries(t)))))
  cat(sprintf("  %-34s ok | area triangoli %.6f (originale %.6f, scarto relativo %.2e)\n", name, a, A0, a / A0 - 1)); invisible(TRUE)
}
cat("correzioni:\n")
try_fix("nessuna", g)
try_fix("buffer 0", geos_buffer(g, 0))
for (e in c(1e-6, 1e-4, 1e-3)) try_fix(sprintf("erosione %g um", e), geos_buffer(g, -e))
for (p in c(1e-4, 1e-3, 1e-2)) try_fix(sprintf("set_precision %g", p), geos_set_precision(g, p))
for (t in c(1e-2, 1e-1)) try_fix(sprintf("simplify %g", t), geos_simplify_preserve_topology(g, t))
try_fix("senza buchi (solo controllo)", geos_make_polygon(wk_coords(sh)$x, wk_coords(sh)$y))
grid_split <- function(h, step) { b <- geos_extent(h); xs <- seq(b$xmin, b$xmax + step, by = step); ys <- seq(b$ymin, b$ymax + step, by = step); pcs <- list()
  for (x in xs[-length(xs)]) for (y in ys[-length(ys)]) { p <- geos_intersection(h, geos_create_rectangle(x, y, x + step, y + step)); if (!geos_is_empty(p) && geos_area(p) > 0) pcs[[length(pcs) + 1]] <- p }
  pcs }
for (st in c(50, 20)) { pcs <- grid_split(g, st); ok <- TRUE; a <- 0
  for (p in pcs) { t <- tryCatch(geos_constrained_delaunay_triangles(p), error = function(e) NULL); if (is.null(t)) { ok <- FALSE; break }; a <- a + sum(geos_area(geos_geometry_n(t, seq_len(geos_num_geometries(t))))) }
  cat(sprintf("  taglio a griglia %g um (%d pezzi)    %s | area %.6f, scarto relativo %.2e\n", st, length(pcs), if (ok) "ok" else "FALLITA", a, a / A0 - 1)) }
