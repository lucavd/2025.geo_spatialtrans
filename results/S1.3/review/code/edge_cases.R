# results/S1.3/review/code/edge_cases.R — revisione del codice S1.3 (RA-code), compito 1: casi limite e tracciamento
# Uso (root del repo): R_LIBS_SITE=/nonexistent R_LIBS_USER=/nonexistent Rscript --vanilla results/S1.3/review/code/edge_cases.R
.libPaths("renv/library/linux-ubuntu-noble/R-4.6/x86_64-pc-linux-gnu")
args <- c("none", "/mnt/micron/geo_spatialtrans/S1.3"); commandArgs <- function(trailingOnly = TRUE) args
source("tools/S1.3_roi.R")
OUTD <- "results/S1.3/review/code"
out <- list(); rec <- function(case, metric, value) { out[[length(out) + 1]] <<- data.frame(case = case, metric = metric, value = as.character(value)); cat(sprintf("%-38s %-45s %s\n", case, metric, value)) }
sq <- function(x0, y0, x1, y1) st_polygon(list(rbind(c(x0, y0), c(x1, y0), c(x1, y1), c(x0, y1), c(x0, y0))))
mkreg <- function(g) list(region_df = data.frame(region_id = 1L, cluster_id = 1L), region_polygons = st_sfc(g))
mkcen <- function(x, y) data.frame(cell_id = seq_along(x), region_id = 1L, x = x, y = y)
TVs <- function(...) tryCatch(tessellate_voronoi(..., verbose = FALSE), error = function(e) e)
summ <- function(lab, tv, A) {
  if (inherits(tv, "error")) { rec(lab, "errore", conditionMessage(tv)); return(invisible()) }
  rec(lab, "aree", paste(sprintf("%.10g", tv$territory_df$territory_area), collapse = ";"))
  rec(lab, "rel |A - Σ|", sprintf("%.3g", abs(A - sum(tv$territory_df$territory_area)) / A))
  rec(lab, "tipi", paste(as.character(st_geometry_type(tv$cell_territories)), collapse = ";"))
  rec(lab, "frammenti/riass./agganci/multi/rip.", with(tv$info, sprintf("%d/%d/%d/%d/%d", n_fragments, n_fragments_reassigned, n_snapped, n_multipart, n_repaired)))
}
R1 <- sq(0, 0, 100, 60)
summ("3 punti collineari", TVs(mkcen(c(20, 50, 80), c(30, 30, 30)), mkreg(R1)), 6000)
summ("4 punti cocircolari (vertice grado 4)", TVs(mkcen(c(25, 75, 25, 75), c(15, 15, 45, 45)), mkreg(R1)), 6000)
summ("generatore sul bordo della regione", TVs(mkcen(c(0, 50), c(30, 30)), mkreg(R1)), 6000)
summ("generatore su un vertice della regione", TVs(mkcen(c(0, 50), c(0, 30)), mkreg(R1)), 6000)
isl <- st_multipolygon(list(unclass(sq(0, 0, 100, 60)), unclass(sq(120, 0, 130, 10))))
summ("MULTIPOLYGON con isola senza generatore", TVs(mkcen(c(30, 90), c(30, 30)), mkreg(isl)), 6100)
# orfano che tocca il ricevente solo in un punto: L di due quadrati che si toccano per un vertice
two <- st_multipolygon(list(unclass(sq(0, 0, 50, 50)), unclass(sq(50, 50, 100, 100))))
summ("regione: 2 quadrati che si toccano in 1 punto", TVs(mkcen(c(10, 40, 60), c(10, 45, 90)), mkreg(two)), 5000)
# smussatura su territori multiparte (A5 r1: 2 orfani senza confinante restano al generatore)
w <- r3_roi_window("A5", "r1", rois); gen <- as.data.frame(read_parquet(file.path(R3, "real", "A5_r1_cellpose_rgb_gen.parquet")))
cen <- mkcen(gen$x, gen$y); t0 <- tessellate_voronoi(cen, one_region(w), verbose = FALSE)
mp <- which(as.character(st_geometry_type(t0$cell_territories)) != "POLYGON")
t1 <- tessellate_voronoi(cen, one_region(w), corner_smoothing = 1/3, verbose = FALSE)
rec("A5 r1 multiparte, cs=1/3", "celle multiparte (cs=0)", paste(mp, collapse = ";"))
rec("A5 r1 multiparte, cs=1/3", "area cs=0", paste(sprintf("%.6g", t0$territory_df$territory_area[mp]), collapse = ";"))
rec("A5 r1 multiparte, cs=1/3", "area cs=1/3", paste(sprintf("%.6g", t1$territory_df$territory_area[mp]), collapse = ";"))
rec("A5 r1 multiparte, cs=1/3", "tipo cs=1/3", paste(as.character(st_geometry_type(t1$cell_territories[mp])), collapse = ";"))
rec("A5 r1 multiparte, cs=1/3", "n_smooth_pieces_dropped (tutte)", t1$info$n_smooth_pieces_dropped)
# st_make_valid nella smussatura (non contato in n_repaired): quanti anelli smussati sono invalidi, cs = 1/3, 2/3, 1
E <- tv_env()
for (lab in c("A1_r1", "A5_r1", "A6_r1")) {
  A <- sub("_.*", "", lab); roi <- sub(".*_", "", lab); a <- ARCH[ARCH$archetype == A, ]
  w <- r3_roi_window(A, roi, rois); gen <- as.data.frame(read_parquet(file.path(R3, "real", sprintf("%s_%s_%s_gen.parquet", A, roi, a$primary))))
  t0 <- tessellate_voronoi(mkcen(gen$x, gen$y), one_region(w), verbose = FALSE)
  for (ni in 1:3) {
    bad <- sum(vapply(t0$cell_territories, function(g) {
      sm <- if (inherits(g, "POLYGON")) st_polygon(lapply(unclass(g), E$.tv_chaikin_ring, n_iter = ni)) else
        st_multipolygon(lapply(unclass(g), function(p) lapply(p, E$.tv_chaikin_ring, n_iter = ni)))
      !isTRUE(st_is_valid(st_sfc(sm)))
    }, TRUE))
    ts <- tessellate_voronoi(mkcen(gen$x, gen$y), one_region(w), corner_smoothing = c(1/3, 2/3, 1)[ni], verbose = FALSE)
    rec(sprintf("%s smussatura n_iter=%d", lab, ni), "Chaikin invalidi -> st_make_valid non contato", sprintf("%d su %d; info$n_repaired = %d", bad, length(t0$cell_territories), ts$info$n_repaired))
  }
}
# tracciamento della riassegnazione su tutti i 40 ROI reali (GEOS): lunghezza condivisa, unione multiparte, aggancio riuscito
E <- tv_env(); src <- deparse(E$.tv_fragments)
i_u <- grep("geoms\\[\\[j\\]\\] <- u; bbm\\[j, \\] <- bb4\\(u\\)|geoms\\[\\[j\\]\\] <- u", src)[1]
src[i_u] <- paste0(src[i_u], "; TR[[length(TR) + 1]] <<- c(len = max(len), multi = as.numeric(inherits(u, 'MULTIPOLYGON') && length(unclass(u)) > 1), snapped = as.numeric(n_snap > SN0)); SN0 <<- n_snap")
f <- eval(parse(text = src)); environment(f) <- E; E$.tv_fragments <- f
jobs <- do.call(rbind, lapply(seq_len(nrow(rois)), function(i) { a <- ARCH[ARCH$archetype == rois$archetype[i], ]
  data.frame(A = rois$archetype[i], roi = rois$roi_id[i], m = c(a$primary, if (!is.na(a$sens)) a$sens), stringsAsFactors = FALSE) }))
trs <- mclapply(seq_len(nrow(jobs)), function(k) {
  TR <<- list(); SN0 <<- 0L
  w <- r3_roi_window(jobs$A[k], jobs$roi[k], rois); gen <- as.data.frame(read_parquet(file.path(R3, "real", sprintf("%s_%s_%s_gen.parquet", jobs$A[k], jobs$roi[k], jobs$m[k]))))
  tv <- E$tessellate_voronoi(mkcen(gen$x, gen$y), one_region(w), verbose = FALSE)
  tr <- if (length(TR)) as.data.frame(do.call(rbind, TR)) else data.frame(len = numeric(0), multi = numeric(0), snapped = numeric(0))
  data.frame(A = jobs$A[k], roi = jobs$roi[k], m = jobs$m[k], n_frag = tv$info$n_fragments, n_reass = tv$info$n_fragments_reassigned,
             n_snapped = tv$info$n_snapped, n_multipart = tv$info$n_multipart, n_union_multi_after = sum(tr$multi == 1),
             n_snap_ok = sum(tr$snapped == 1 & tr$multi == 0), n_snap_failed = sum(tr$snapped == 1 & tr$multi == 1),
             n_len_lt_1e6 = sum(tr$len < 1e-6), n_len_lt_1e6_multi = sum(tr$len < 1e-6 & tr$multi == 1), min_len = if (nrow(tr)) min(tr$len) else NA)
}, mc.cores = 30L, mc.preschedule = FALSE)
bad <- vapply(trs, function(z) !is.data.frame(z), TRUE); if (any(bad)) print(trs[bad])
TRt <- do.call(rbind, trs[!bad]); write.csv(TRt, file.path(OUTD, "fragments_trace_40roi.csv"), row.names = FALSE)
print(colSums(TRt[, c("n_frag", "n_reass", "n_snapped", "n_multipart", "n_union_multi_after", "n_snap_ok", "n_snap_failed", "n_len_lt_1e6", "n_len_lt_1e6_multi")]))
write.csv(do.call(rbind, out), file.path(OUTD, "edge_cases.csv"), row.names = FALSE)
